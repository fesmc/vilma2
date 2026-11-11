module vilma_tensor_sh
   !! Tensor spherical-harmonic dyadic transforms (Martinec 2000, Appendix B), the
   !! machinery rung 6 (laterally-varying viscosity) needs — GENERAL order (mmax≥0).
   !!
   !! The Maxwell memory τ and strain ε are second-order symmetric tensors on the
   !! sphere, expanded in the tensor spherical harmonics Z^λ: the spheroidal
   !! λ∈{1,2,5,6} and, once viscosity varies laterally, the toroidal λ∈{3,4}.
   !! Coefficient arrays carry either TLAM_SPH = 4 channels (spheroidal only, the
   !! 1-D case) or TLAM = 6, and the transforms dispatch on which: local channel
   !! 1..4 = λ 1,2,5,6 as always, 5 = λ3, 6 = λ4 appended, so the spheroidal
   !! subset is always the leading four. For lateral viscosity the update τ⁺=(1−M)τ−2μM·ε is
   !! pointwise in PHYSICAL space, so the tensor is reconstructed on the grid via its
   !! six dyadic components (eqs 90/91, B10/B11) and projected back. The six physical
   !! components rr, rθ, rφ, θθ, θφ, φφ relate to the coefficients by
   !!   rr        = Σ T¹ Y                          (Z¹)
   !!   rθ, rφ    = Σ T² (E, F)                      (Z²,  E=∂_θY, F=(1/sinθ)∂_φY)
   !!   θθ        = Σ [−l(l+1)T⁵ Y + T⁶ G]          (Z⁵ trace + Z⁶)
   !!   φφ        = Σ [−l(l+1)T⁵ Y − T⁶ G]
   !!   θφ        = Σ 4 T⁶ H                          (Z⁶,  G,H = B11 second derivatives)
   !! and the toroidal channels add
   !!   rθ, rφ   += Σ T³ (−F, E)                     (Z³ = e_r×∇₁ companion of Z²)
   !!   θθ, φφ   += Σ ∓T⁴ H,   θφ += Σ T⁴ G         (Z⁴ = Z⁶ turned by 45°, halved)
   !! Z³ = sym(e_r ⊗ e_r×∇₁Y) and Z⁴ = sym(∇₁(e_r×∇₁Y)), so the toroidal
   !! displacement u = W(r) e_r×∇₁Y has strain (W′−W/r) Z³ + (W/r) Z⁴ —
   !! the λ=3,4 rows of Martinec eq 87, re-derived rather than transcribed.
   !!
   !! Synthesis is EXACT via grid identities — no recurrence, no re-analysis. The
   !! spin-2 G,H (the only pieces SHTns has no routine for) come from scalar + vector
   !! synths plus algebraic grid factors, using ∂_φ = "multiply coeffs by im" (exact
   !! on the known input) and ∂_θθ via ∇₁²Y = −l(l+1)Y:
   !!   Sg ≡ Σ T⁶ G = ∇₁²f − 2cotθ·g_θ − 2(1/sinθ)·∂_φ g_φ
   !!   Sh ≡ Σ T⁶ H = (1/sinθ)·∂_φ g_θ − cotθ·g_φ
   !! with f=synth(T⁶), (g_θ,g_φ)=sph_synth(T⁶), ∂_φ(·)=sph_synth(im·T⁶).
   !!
   !! Analysis: channels 1,2,5 (and 3) invert through SHTns's own scalar/vector
   !! analyses (a synth/analysis pair) with the −l(l+1) factor for the trace; the
   !! spin-2 channels 6 and 4 use the adjoint of their synthesis (`Sg*`,`Sh*`) — the
   !! same ops with synth↔analysis swapped — each normalised by a per-degree factor
   !! calibrated once at init. Z⁴⊥Z⁶ (B12) keeps the two projections diagonal only
   !! as far as the grid quadrature is exact; test_tensor_sh measures the leak.
   !! Validated by the round trip and the physical ∫τ:ε double-dot vs the B13 norms.
   use vilma_precision, only: wp
   use vilma_sht,       only: sht_grid, sht_free_cfg, sht_grid_lmidx, sht_grid_clone_cfg, sht_grid_synthesis, sht_grid_sph_synthesis, sht_grid_analysis, sht_grid_sph_analysis, &
                           sht_grid_sphtor_synthesis, sht_grid_sphtor_analysis
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated
   !$ use omp_lib
   implicit none
   private

   public :: tensor_sh
   public :: tensor_sh_init, tensor_sh_synth, tensor_sh_analysis, tensor_sh_thread_cfg, tensor_sh_destroy
   integer, parameter, public :: TLAM_SPH = 4      ! λ = 1,2,5,6 → local 1..4
   integer, parameter, public :: TLAM     = 6      ! + λ = 3,4   → local 5,6
   ! Dyadic-field plane indices (the third dimension of the dyad array).
   integer, parameter, public :: DY_RR = 1, DY_RT = 2, DY_RP = 3, &
                                 DY_TT = 4, DY_TP = 5, DY_PP = 6

   ! Scratch of one transform: grid fields and coefficient vectors. A synthesis
   ! uses the grid slots G_* below, an analysis the first five; they never run at
   ! the same time on one thread.
   integer, parameter :: NGRID = 11, NCOEF = 3
   integer, parameter :: G_TR = 1, G_TT = 2, G_PP = 3, G_TP = 4, &   ! tensor_sh_synth
                         G_SG = 5, G_SH = 6, &                       ! spin2_synth(_tor)
                         G_GT = 7, G_GP = 8, G_GTF = 9, G_GPF = 10, G_LAP = 11   ! spin2_fields
   integer, parameter :: G_VT = 1, G_VP = 2, &                       ! tensor_sh_analysis
                         G_D = 3, G_AT = 4, G_AP = 5                 ! spin2_adjoint(_tor)

   type :: tensor_sh_thread
      !! What one OpenMP thread transforms with: its own SHTns config (a single
      !! config is NOT safe for concurrent calls) and its own scratch. The scratch
      !! is on the heap, allocated by the thread at its first transform (so the
      !! pages are local to it): as automatic arrays the ~13 grid fields of a
      !! transform overflowed the default thread stack from lmax ~ 128 on.
      type(c_ptr) :: cfg = c_null_ptr
      real(wp),    allocatable :: g(:,:,:)   !! (nphi, nlat, NGRID) grid scratch
      complex(wp), allocatable :: c(:,:)     !! (nlm, NCOEF) coefficient scratch
   end type tensor_sh_thread

   type :: tensor_sh
      integer :: lmax = 0, nlm = 0, nphi = 0, nlat = 0
      integer,  allocatable :: ldeg(:)   !! (nlm) degree l of each coefficient
      integer,  allocatable :: mord(:)   !! (nlm) order m of each coefficient
      real(wp), allocatable :: llp1(:)   !! (nlm) l(l+1)
      real(wp), allocatable :: cott(:)   !! (nlat) cotθ at the Gauss latitudes
      real(wp), allocatable :: invsin(:) !! (nlat) 1/sinθ
      real(wp), allocatable :: n6(:)     !! (0:lmax) spin-2 channel norm S₆*S₆ (calibrated)
      real(wp), allocatable :: n4(:)     !! (0:lmax) toroidal spin-2 norm S₄*S₄ (calibrated)
      ! One tensor_sh_thread per OpenMP thread, so the element loop in the memory
      ! advance can run the dyadic transforms concurrently. Built serially at init.
      ! A pointer, so that a transform can use its thread's scratch through an
      ! intent(in) tensor_sh; each thread touches its own element only.
      type(tensor_sh_thread), pointer :: pool(:) => null()   !! (npool)
      integer :: npool = 0
   end type tensor_sh

contains

   subroutine tensor_sh_init(self, sht)
      !! Set the transformer up for the grid sht. Idempotent: a transformer that
      !! already serves this grid and thread count is left as it is (its configs
      !! and calibration are functions of those alone), and any other one is
      !! destroyed first, so no SHTns config is ever dropped without being freed.
      type(tensor_sh), intent(inout) :: self
      type(sht_grid),   intent(in)    :: sht
      complex(wp), allocatable :: c6(:), craw(:)
      real(wp),    allocatable :: tt(:,:), pp(:,:), tp(:,:)
      integer :: l, m, lm, i, npool

      npool = 1
      !$ npool = omp_get_max_threads()
      if (associated(self%pool)) then
         if (self%lmax == sht%lmax .and. self%nlm == sht%nlm .and. self%nphi == sht%nphi &
             .and. self%nlat == sht%nlat .and. self%npool == npool) return
      end if
      call tensor_sh_destroy(self)

      self%lmax = sht%lmax;  self%nlm = sht%nlm
      self%nphi = sht%nphi;  self%nlat = sht%nlat
      allocate(self%ldeg(self%nlm), self%mord(self%nlm), self%llp1(self%nlm))
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, sht%lmax
            lm = sht_grid_lmidx(sht, l, m)
            self%ldeg(lm) = l;  self%mord(lm) = m
            self%llp1(lm) = real(l,wp)*real(l+1,wp)
         end do
      end do
      allocate(self%cott(self%nlat), self%invsin(self%nlat))
      do i = 1, self%nlat
         self%cott(i)   = cos(sht%colat(i))/sin(sht%colat(i))
         self%invsin(i) = 1.0_wp/sin(sht%colat(i))
      end do

      ! Per-thread configs (serial creation — FFTW planning is not thread-safe).
      ! Before the calibration, which transforms with this thread's scratch.
      self%npool = npool
      allocate(self%pool(self%npool))
      do i = 1, self%npool
         self%pool(i)%cfg = sht_grid_clone_cfg(sht)
      end do

      ! Calibrate the spin-2 channel norms per degree: send a unit Z⁶ (Z⁴)
      ! coefficient through synthesis then its adjoint; the diagonal response is
      ! n6(l) (n4(l)). The norm depends on the degree only, so one order per degree
      ! is enough: the highest the grid carries, m = l where it can (sectoral: both
      ! G and H present), and m = 0 (H ≡ 0, only G contributes) for a degree that
      ! has no other — every degree of an axisymmetric grid (mmax = 0), and the
      ! degrees below mres of a grid with an order stride.
      allocate(self%n6(0:self%lmax));  self%n6 = 0.0_wp
      allocate(self%n4(0:self%lmax));  self%n4 = 0.0_wp
      allocate(c6(self%nlm), craw(self%nlm))
      allocate(tt(self%nphi,self%nlat), pp(self%nphi,self%nlat), tp(self%nphi,self%nlat))
      do l = 2, self%lmax
         ! Round DOWN to a multiple of mres: shtns_lmidx assumes it, and with
         ! mres > 1 the bare min() would name an order the config does not carry.
         m = (min(l, sht%mmax*sht%mres)/sht%mres)*sht%mres
         lm = sht_grid_lmidx(sht, l, m)
         c6 = (0.0_wp, 0.0_wp);  c6(lm) = (1.0_wp, 0.0_wp)
         call spin2_synth(self, sht, c6, tt, pp, tp)      ! Sg→±(tt,pp), 4Sh→tp
         call spin2_adjoint(self, sht, tt, pp, tp, craw)  ! raw S₆* (unnormalised)
         self%n6(l) = real(craw(lm), wp)
         call spin2_synth_tor(self, sht, c6, tt, pp, tp)  ! the same coefficient as Z⁴
         call spin2_adjoint_tor(self, sht, tt, pp, tp, craw)
         self%n4(l) = real(craw(lm), wp)
      end do
      deallocate(c6, craw, tt, pp, tp)
   end subroutine tensor_sh_init

   function thread_slot(self) result(t)
      !! The calling OpenMP thread's element of the pool (1-based index = thread
      !! id + 1; serial / non-OpenMP builds always get pool(1)), with its scratch
      !! allocated.
      type(tensor_sh), intent(in) :: self
      type(tensor_sh_thread), pointer :: t
      integer :: tid
      tid = 0
      !$ tid = omp_get_thread_num()
      ! The pool is sized at init from omp_get_max_threads(). A host that raises
      ! the thread count afterwards, or that calls us from inside its OWN parallel
      ! region (where nested parallelism is off, so every host thread reports
      ! tid = 0 and would share pool(1)), breaks the one-config-per-thread
      ! contract this exists to keep. Fail loudly rather than transform on a
      ! garbage handle or silently share one config.
      if (tid < 0 .or. tid >= self%npool) error stop &
         'vilma_tensor_sh: thread id outside the config pool — the thread &
         &count changed after tensor_sh_init, or this was called from a nested &
         &parallel region'
      t => self%pool(tid+1)
      if (.not. allocated(t%g)) then
         allocate(t%g(self%nphi, self%nlat, NGRID), t%c(self%nlm, NCOEF))
         t%g = 0.0_wp;  t%c = (0.0_wp, 0.0_wp)         ! first touch, by this thread
      end if
   end function thread_slot

   function tensor_sh_thread_cfg(self) result(cfg)
      !! The calling OpenMP thread's private SHTns config (see thread_slot).
      type(tensor_sh), intent(in) :: self
      type(c_ptr) :: cfg
      type(tensor_sh_thread), pointer :: t
      t => thread_slot(self)
      cfg = t%cfg
   end function tensor_sh_thread_cfg

   ! --- synthesis -------------------------------------------------------------

   subroutine tensor_sh_synth(self, sht, c, dyad, cfg)
      !! Tensor-harmonic coefficients c(λ, nlm) → six dyadic grid fields, with
      !! size(c,1) = TLAM_SPH (spheroidal) or TLAM (plus toroidal λ=3,4).
      !! Pass `cfg` (a thread_cfg handle) to transform on a thread-local config.
      type(tensor_sh), intent(in)  :: self
      type(sht_grid),   intent(in)  :: sht
      complex(wp),      intent(in)  :: c(:,:)        !! (TLAM_SPH or TLAM, nlm)
      real(wp),         intent(out) :: dyad(:,:,:)   !! (nphi, nlat, 6)
      type(c_ptr), intent(in), optional :: cfg
      type(tensor_sh_thread), pointer :: w
      logical     :: tor
      tor = toroidal_channels(c)
      w => thread_slot(self)
      associate (tr => w%g(:,:,G_TR), tt => w%g(:,:,G_TT), pp => w%g(:,:,G_PP), &
                 tp => w%g(:,:,G_TP), scaled => w%c(:,1))
      ! rr (Z¹) and rθ,rφ (Z², + Z³ when toroidal)
      call sht_grid_synthesis(sht, c(1,:), dyad(:,:,DY_RR), cfg)
      if (tor) then
         call sht_grid_sphtor_synthesis(sht, c(2,:), c(5,:), dyad(:,:,DY_RT), dyad(:,:,DY_RP), cfg)
      else
         call sht_grid_sph_synthesis(sht, c(2,:), dyad(:,:,DY_RT), dyad(:,:,DY_RP), cfg)
      end if
      ! trace from Z⁵:  −l(l+1) Y
      scaled = -self%llp1*c(3,:)
      call sht_grid_synthesis(sht, scaled, tr, cfg)
      ! spin-2 from Z⁶
      call spin2_synth(self, sht, c(4,:), tt, pp, tp, cfg)   ! tt=Sg, pp=−Sg, tp=4Sh
      dyad(:,:,DY_TT) = tr + tt
      dyad(:,:,DY_PP) = tr - tt
      dyad(:,:,DY_TP) = tp
      if (tor) then                                           ! + Z⁴: (−Sh, +Sh, Sg)
         call spin2_synth_tor(self, sht, c(6,:), tt, pp, tp, cfg)
         dyad(:,:,DY_TT) = dyad(:,:,DY_TT) + tt
         dyad(:,:,DY_PP) = dyad(:,:,DY_PP) + pp
         dyad(:,:,DY_TP) = dyad(:,:,DY_TP) + tp
      end if
      end associate
   end subroutine tensor_sh_synth

   logical function toroidal_channels(c) result(tor)
      !! Whether a coefficient block carries the toroidal channels. Anything but
      !! the two supported widths is a caller error, and a silent one: a block of
      !! 5 would drop Z⁴ and a block of 3 would read past Z⁵.
      complex(wp), intent(in) :: c(:,:)
      select case (size(c,1))
      case (TLAM_SPH);  tor = .false.
      case (TLAM);      tor = .true.
      case default
         error stop 'vilma_tensor_sh: coefficient block must have TLAM_SPH or TLAM channels'
      end select
   end function toroidal_channels

   subroutine spin2_synth(self, sht, c6, tt, pp, tp, cfg)
      !! Z⁶ contribution: tt=Sg=Σc6·G, pp=−Sg, tp=4Sh=4Σc6·H, via the exact grid
      !! identities (no recurrence, no re-analysis).
      type(tensor_sh), intent(in)  :: self
      type(sht_grid),   intent(in)  :: sht
      complex(wp),      intent(in)  :: c6(:)
      real(wp),         intent(out) :: tt(:,:), pp(:,:), tp(:,:)
      type(c_ptr), intent(in), optional :: cfg
      type(tensor_sh_thread), pointer :: w
      w => thread_slot(self)
      associate (sg => w%g(:,:,G_SG), sh => w%g(:,:,G_SH))
      call spin2_fields(self, sht, c6, sg, sh, cfg)
      tt = sg;  pp = -sg;  tp = 4.0_wp*sh
      end associate
   end subroutine spin2_synth

   subroutine spin2_synth_tor(self, sht, c4, tt, pp, tp, cfg)
      !! Z⁴ contribution, Z⁴ = G e_θφ − H (e_θθ − e_φφ): tt=−Sh, pp=+Sh, tp=Sg,
      !! with Sg, Sh the same grid fields as Z⁶ takes, of the coefficients c4.
      type(tensor_sh), intent(in)  :: self
      type(sht_grid),   intent(in)  :: sht
      complex(wp),      intent(in)  :: c4(:)
      real(wp),         intent(out) :: tt(:,:), pp(:,:), tp(:,:)
      type(c_ptr), intent(in), optional :: cfg
      type(tensor_sh_thread), pointer :: w
      w => thread_slot(self)
      associate (sg => w%g(:,:,G_SG), sh => w%g(:,:,G_SH))
      call spin2_fields(self, sht, c4, sg, sh, cfg)
      tt = -sh;  pp = sh;  tp = sg
      end associate
   end subroutine spin2_synth_tor

   subroutine spin2_fields(self, sht, c6, sg, sh, cfg)
      !! The two spin-2 grid fields of coefficients c6: Sg = Σc6·G, Sh = Σc6·H.
      type(tensor_sh), intent(in)  :: self
      type(sht_grid),   intent(in)  :: sht
      complex(wp),      intent(in)  :: c6(:)
      real(wp),         intent(out) :: sg(:,:), sh(:,:)
      type(c_ptr), intent(in), optional :: cfg
      type(tensor_sh_thread), pointer :: w
      integer :: i
      w => thread_slot(self)
      associate (gt => w%g(:,:,G_GT), gp => w%g(:,:,G_GP), gtf => w%g(:,:,G_GTF), &
                 gpf => w%g(:,:,G_GPF), lap => w%g(:,:,G_LAP), imc => w%c(:,2), nlc => w%c(:,3))
      imc = cmplx(0.0_wp, real(self%mord,wp), wp)*c6            ! im·c6  (= ∂_φ on coeffs)
      nlc = -self%llp1*c6                                       ! −l(l+1)·c6
      ! f itself is NOT synthesized: only ∇₁²f enters Sg, and that is the
      ! synthesis of −l(l+1)·c6 below. This routine is called once per
      ! tensor_sh_synth, i.e. two to three times per 3-D element per step, so the
      ! transform it does not do is worth the comment saying why.
      call sht_grid_sph_synthesis(sht, c6, gt, gp, cfg)                   ! g_θ, g_φ
      call sht_grid_sph_synthesis(sht, imc, gtf, gpf, cfg)                ! ∂_φ g_θ, ∂_φ g_φ
      call sht_grid_synthesis(sht, nlc, lap, cfg)                         ! ∇₁²f = −l(l+1)f
      ! Sg = ∇₁²f − 2cotθ g_θ − 2(1/sinθ)∂_φ g_φ ;  Sh = (1/sinθ)∂_φ g_θ − cotθ g_φ
      do i = 1, self%nlat
         sg(:,i) = lap(:,i) - 2.0_wp*(gt(:,i)*self%cott(i)) - 2.0_wp*(gpf(:,i)*self%invsin(i))
         sh(:,i) = gtf(:,i)*self%invsin(i) - gp(:,i)*self%cott(i)
      end do
      end associate
   end subroutine spin2_fields

   ! --- analysis --------------------------------------------------------------

   subroutine tensor_sh_analysis(self, sht, dyad, c, cfg)
      !! Six dyadic grid fields → tensor-harmonic coefficients, spheroidal only
      !! or with the toroidal channels according to size(c,1) (see synth).
      !! Pass `cfg` (a thread_cfg handle) to transform on a thread-local config.
      type(tensor_sh), intent(in)    :: self
      type(sht_grid),   intent(in)    :: sht
      real(wp),         intent(inout) :: dyad(:,:,:)   !! (nphi,nlat,6); SHTns overwrites
      complex(wp),      intent(out)   :: c(:,:)        !! (TLAM_SPH or TLAM, nlm)
      type(c_ptr), intent(in), optional :: cfg
      type(tensor_sh_thread), pointer :: w
      integer     :: lm
      logical     :: tor
      tor = toroidal_channels(c)
      w => thread_slot(self)
      associate (vt => w%g(:,:,G_VT), vp => w%g(:,:,G_VP), craw => w%c(:,1))
      ! rr (Z¹): scalar analysis (inverse of synth)
      call sht_grid_analysis(sht, dyad(:,:,DY_RR), c(1,:), cfg)
      ! rθ,rφ (Z², Z³): vector analysis (inverse of the sph/sphtor synth)
      vt = dyad(:,:,DY_RT);  vp = dyad(:,:,DY_RP)
      if (tor) then
         call sht_grid_sphtor_analysis(sht, vt, vp, c(2,:), c(5,:), cfg)
      else
         call sht_grid_sph_analysis(sht, vt, vp, c(2,:), cfg)
      end if
      ! trace (Z⁵): analysis(θθ+φφ) = −2 l(l+1) T⁵
      vt = dyad(:,:,DY_TT) + dyad(:,:,DY_PP)
      call sht_grid_analysis(sht, vt, craw, cfg)
      do lm = 1, self%nlm
         if (self%llp1(lm) > 0.0_wp) then
            c(3,lm) = craw(lm)/(-2.0_wp*self%llp1(lm))
         else
            c(3,lm) = (0.0_wp, 0.0_wp)
         end if
      end do
      ! spin-2 (Z⁶): adjoint of spin2_synth, normalised by the calibrated per-degree n6
      call spin2_adjoint(self, sht, dyad(:,:,DY_TT), dyad(:,:,DY_PP), dyad(:,:,DY_TP), craw, cfg)
      do lm = 1, self%nlm
         if (self%n6(self%ldeg(lm)) /= 0.0_wp) then
            c(4,lm) = craw(lm)/self%n6(self%ldeg(lm))
         else
            c(4,lm) = (0.0_wp, 0.0_wp)
         end if
      end do
      if (tor) then
         ! toroidal spin-2 (Z⁴): its own adjoint, normalised by the calibrated n4
         call spin2_adjoint_tor(self, sht, dyad(:,:,DY_TT), dyad(:,:,DY_PP), dyad(:,:,DY_TP), craw, cfg)
         do lm = 1, self%nlm
            if (self%n4(self%ldeg(lm)) /= 0.0_wp) then
               c(6,lm) = craw(lm)/self%n4(self%ldeg(lm))
            else
               c(6,lm) = (0.0_wp, 0.0_wp)
            end if
         end do
      end if
      end associate
   end subroutine tensor_sh_analysis

   subroutine spin2_adjoint(self, sht, dtt, dpp, dtp, craw, cfg)
      !! Unnormalised adjoint of spin2_synth: craw = S_g*(θθ−φφ) + S_h*(θφ), where the
      !! *-operators swap each forward op (synth↔analysis, grid-multiply self-adjoint,
      !! im → −im). The ∫dΩ adjoint of synth IS analysis (= Wᵀ·synth), so this returns
      !! ∫(Z⁶-basis):(reconstructed tensor) — the projection numerator. dtt/dpp/dtp are
      !! overwritten by the SHTns analyses.
      type(tensor_sh), intent(in)    :: self
      type(sht_grid),   intent(in)    :: sht
      real(wp),         intent(inout) :: dtt(:,:), dpp(:,:), dtp(:,:)
      complex(wp),      intent(out)   :: craw(:)
      type(c_ptr), intent(in), optional :: cfg
      type(tensor_sh_thread), pointer :: w
      w => thread_slot(self)
      associate (D => w%g(:,:,G_D), vt => w%g(:,:,G_AT), vp => w%g(:,:,G_AP), &
                 q => w%c(:,2), s => w%c(:,3))
      ! The ∫dΩ adjoint of the spheroidal vector synth (SHsph_to_spat) is l(l+1)·
      ! spat_to_SHsphtor (its INVERSE differs from its adjoint by the spheroidal norm
      ! l(l+1)); the scalar synth's adjoint is plain analysis (orthonormal). So every
      ! vector-analysis result below is scaled by llp1 to be the true adjoint.
      D = dtt - dpp                                   ! θθ−φφ feeds S_g*
      ! S_g*(D) = −l(l+1)·analysis(D) − 2·sphAnal(cotθ·D,0).S + 2 im·sphAnal(0,(1/sinθ)·D).S
      ! (the analysis overwrites its input, so it takes a copy of D)
      vt = D
      call sht_grid_analysis(sht, vt, q, cfg);   craw = -self%llp1*q
      call byprof(D, self%cott, vt);  vp = 0.0_wp
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg);   craw = craw - 2.0_wp*self%llp1*s
      vt = 0.0_wp;  call byprof(D, self%invsin, vp)
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg)
      craw = craw + 2.0_wp*cmplx(0.0_wp, real(self%mord,wp), wp)*self%llp1*s
      ! θφ contributes ∫T:Z⁶|_θφ = a_θφ·4H·(e_θφ:e_θφ=½) = 2 a_θφ H ⇒ 2·S_h*(θφ), with
      ! S_h*(D) = −im·llp1·sphAnal((1/sinθ)D,0).S − llp1·sphAnal(0,cotθ·D).S.
      call byprof(dtp, self%invsin, vt);  vp = 0.0_wp
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg)
      craw = craw - 2.0_wp*cmplx(0.0_wp, real(self%mord,wp), wp)*self%llp1*s
      vt = 0.0_wp;  call byprof(dtp, self%cott, vp)
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg)
      craw = craw - 2.0_wp*self%llp1*s
      end associate
   end subroutine spin2_adjoint

   subroutine spin2_adjoint_tor(self, sht, dtt, dpp, dtp, craw, cfg)
      !! Unnormalised adjoint of spin2_synth_tor, the Z⁴ projection numerator
      !! ∫(Z⁴-basis):(tensor). With Z⁴ = G e_θφ − H(e_θθ − e_φφ) and e_θφ:e_θφ = ½,
      !!   ∫τ:Z⁴ = ½·S_g*(θφ) − S_h*(θθ−φφ),
      !! the Z⁶ projection with the two spin-2 inputs exchanged. S_g*, S_h* are the
      !! operators written out in spin2_adjoint; the inputs are read, not written.
      type(tensor_sh), intent(in)    :: self
      type(sht_grid),   intent(in)    :: sht
      real(wp),         intent(in)    :: dtt(:,:), dpp(:,:), dtp(:,:)
      complex(wp),      intent(out)   :: craw(:)
      type(c_ptr), intent(in), optional :: cfg
      type(tensor_sh_thread), pointer :: w
      w => thread_slot(self)
      associate (D => w%g(:,:,G_D), vt => w%g(:,:,G_AT), vp => w%g(:,:,G_AP), &
                 q => w%c(:,2), s => w%c(:,3))
      ! ½·S_g*(θφ) = ½·[−llp1·analysis(D) − 2·llp1·sphAnal(cotθ D,0) + 2 im·llp1·sphAnal(0,D/sinθ)]
      ! (im = i·m below; the analysis overwrites its input, so it takes a copy)
      vt = dtp
      call sht_grid_analysis(sht, vt, q, cfg);   craw = -0.5_wp*self%llp1*q
      call byprof(dtp, self%cott, vt);  vp = 0.0_wp
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg);   craw = craw - self%llp1*s
      vt = 0.0_wp;  call byprof(dtp, self%invsin, vp)
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg)
      craw = craw + cmplx(0.0_wp, real(self%mord,wp), wp)*self%llp1*s
      ! −S_h*(θθ−φφ) = +im·llp1·sphAnal(D/sinθ,0) + llp1·sphAnal(0,cotθ D)
      D = dtt - dpp
      call byprof(D, self%invsin, vt);  vp = 0.0_wp
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg)
      craw = craw + cmplx(0.0_wp, real(self%mord,wp), wp)*self%llp1*s
      vt = 0.0_wp;  call byprof(D, self%cott, vp)
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg);   craw = craw + self%llp1*s
      end associate
   end subroutine spin2_adjoint_tor

   ! --- helpers ---------------------------------------------------------------

   pure subroutine byprof(field, prof, out)
      !! out = each Gauss-latitude column of a (nphi,nlat) field times prof(nlat).
      real(wp), intent(in)  :: field(:,:), prof(:)
      real(wp), intent(out) :: out(:,:)
      integer  :: i
      do i = 1, size(field,2)
         out(:,i) = field(:,i)*prof(i)
      end do
   end subroutine byprof

   subroutine tensor_sh_destroy(self)
      type(tensor_sh), intent(inout) :: self
      integer :: i
      if (allocated(self%ldeg))   deallocate(self%ldeg, self%mord, self%llp1)
      if (allocated(self%cott))   deallocate(self%cott, self%invsin)
      if (allocated(self%n6))     deallocate(self%n6)
      if (allocated(self%n4))     deallocate(self%n4)
      if (associated(self%pool)) then
         do i = 1, self%npool
            if (c_associated(self%pool(i)%cfg)) call sht_free_cfg(self%pool(i)%cfg)
         end do
         deallocate(self%pool)          ! with each thread's scratch
      end if
      nullify(self%pool)
      self%npool = 0;  self%lmax = 0;  self%nlm = 0
   end subroutine tensor_sh_destroy

end module vilma_tensor_sh
