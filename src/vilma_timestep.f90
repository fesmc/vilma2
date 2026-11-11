module vilma_timestep
   !! Adaptive time-stepping strategies for the viscoelastic + sea-level model.
   !! Isolated from vilma_coupling so the model can carry more than one stepping strategy
   !! (fixed sub-steps, adaptive, …) behind a common driver.
   !!
   !! `adaptive_stepper` advances the VE+SLE model across a coupling interval [t0,t1]
   !! with the ice load LINEARLY INTERPOLATED between its endpoints. It carries TWO
   !! sub-stepping strategies, chosen by the memory scheme:
   !!
   !!  • EXPLICIT (forward-Euler, the default): the 1st-order memory has no embedded
   !!    error signal, so Δt is set a priori by the Maxwell STABILITY ceiling
   !!    Δt ≤ cfl/max(μ/η) and the interval is divided into equal sub-steps (this is
   !!    how VILMA1 marches its explicit memory). The ceiling is exact — max(μ/η) is taken
   !!    over every element and, with lateral viscosity, every grid point — so the scheme
   !!    is stable by construction; a sub-step whose memory ∞-norm comes back non-finite
   !!    (a corrupt input or state) is rolled back and halved, and the run stops when
   !!    halving does not cure it.
   !!
   !!  • IMPLICIT (trapezoidal, SCHEME_TRAP): Δt is chosen by a step-doubling local-error
   !!    estimate on the memory (the §3c adaptive-dt controller) + the SLE↔memory
   !!    co-convergence (§3c 3b): each candidate step is taken once at Δt and once as two
   !!    Δt/2 sub-steps; the coarse/fine difference is the local error (∝ Δt^{p+1}). A
   !!    rejected step rolls back and retries smaller; an accepted step keeps the fine
   !!    state and grows Δt.
   !!
   !! Δt enters the response only through Mk = (μ/η)Δt, so changing it (resp%set_dt) is a
   !! cheap rescale — no operator re-factorization (the band LU is Δt-independent).
   use vilma_precision,    only: wp, is_finite
   use vilma_constants,    only: sec_per_year
   use vilma_sht,          only: sht_grid
   use vilma_response,     only: response_coarse_fine_error, response_stash_coarse, response_prime_sigma, response_restore_state, response_set_dt, response_save_state, response_memory_norm, response_max_rate, response
   use vilma_sle,          only: sle_solve, sle_solver, sle_result, sle_host
   use vilma_rotation,     only: rotation_state, rotation_save_state, rotation_restore_state
   use vilma_viscoelastic, only: scheme_order, scheme_is_implicit
   implicit none
   private

   public :: adaptive_stepper
   public :: stepper_advance

   type :: adaptive_stepper
      !! Embedded step-doubling controller. Tolerances are on the memory ∞-norm (the
      !! surface observable is algebraic in the memory τ, so this also bounds the
      !! observable's local error). Defaults are placeholders — set per problem.
      real(wp) :: rtol      = 1.0e-4_wp   !! relative local-error tolerance (memory norm)
      real(wp) :: atol      = 1.0e-3_wp   !! absolute local-error floor (memory units)
      real(wp) :: safety    = 0.9_wp      !! step-size safety factor
      real(wp) :: grow_max  = 5.0_wp      !! max Δt growth per accepted step
      real(wp) :: shrink_min = 0.2_wp     !! min Δt shrink per step (reject or grow)
      real(wp) :: dt_min    = 0.0_wp      !! Δt floor (0 = none); a step at the floor is
                                          !! accepted even if over tolerance (no infinite
                                          !! subdivision) and tallied in n_floor. For the
                                          !! explicit scheme it OVERRIDES the stability
                                          !! ceiling when larger (dt_min_overrides)
      real(wp) :: dt_max    = huge(1.0_wp)!! Δt ceiling
      !! Explicit (forward-Euler) sub-stepping. The 1st-order memory carries no
      !! embedded local-error signal, so Δt is set a priori by the Maxwell stability
      !! ceiling Δt ≤ cfl/max(μ/η): the interval is split into equal sub-steps with
      !! Maxwell number M = μΔt/η ≤ cfl. cfl = 1 matches the rotation channel's
      !! convention (M ≤ 1, a factor-2 margin below the |1−M| = 1 stability bound).
      real(wp) :: cfl       = 1.0_wp      !! explicit sub-step Maxwell-number ceiling
      !! Set (and a warning printed, once) when dt_min exceeded the explicit stability
      !! ceiling cfl/max(μ/η) and so set the sub-step: the run is then outside the
      !! Maxwell number the user asked for, and unstable beyond M = 2.
      logical  :: dt_min_overrides = .false.
      real(wp) :: dt_try    = 0.0_wp      !! next-step Δt suggestion (carried across
                                          !! intervals; 0 ⇒ first guess = whole interval)
      integer  :: n_accept  = 0           !! accepted steps (cumulative)
      integer  :: n_reject  = 0           !! rejected attempts (cumulative)
      integer  :: n_floor   = 0           !! steps accepted at dt_min over tolerance
      integer  :: n_solve   = 0           !! SLE solves issued (the dominant cost unit)
      real(wp) :: worst_mass_resid = 0.0_wp  !! worst SLE mass residual over the LAST advance()
      !! Convergence of the SLE solves of the LAST advance() (sle_result): how many there
      !! were, how many ended at a cap of the inner (water load), coastline or memory
      !! iteration, the worst inner residual max|ΔS|/max|S| (sle%tol bounds it when
      !! converged) and the most cells still changing ocean/land state at the end of a solve.
      integer  :: n_sle       = 0
      integer  :: n_cap_inner = 0
      integer  :: n_cap_coast = 0
      integer  :: n_cap_mem   = 0
      real(wp) :: worst_resid = 0.0_wp
      integer  :: worst_coast_flip = 0
      !! PROFILE: wall-clock [s] of the rollback bookkeeping — response_save_state /
      !! _restore_state (a full copy of the (NLAM,ne,nk) memory arrays) and the
      !! response_memory_norm scans that feed the guards. This is paid on EVERY
      !! sub-step, outside sle_solve, so without a timer it lands in the driver's
      !! residual bucket and reads as SLE cost.
      real(wp) :: t_guard = 0.0_wp
   end type adaptive_stepper

contains

   subroutine stepper_advance(self, sht, resp, sle, topo0, ice0, ice1, ice_ref, &
                              t0, t1, rsl, C, rot, sigma_out, u_out, N_out, dphi_out, &
                              host_ref, host0, host1, dIg_out, fill_out)
      !! Advance the VE+SLE model from t0 to t1 with the ice load interpolated linearly
      !! ice(t) = ice0 + (t−t0)/(t1−t0)·(ice1−ice0), absolute; the SLE load is
      !! d_ice(t) = ice(t) − ice_ref (total change from the reference state). Δt is
      !! chosen adaptively; on return resp holds the end-of-interval memory at t1 and
      !! rsl/C hold the converged sea-level fields there.
      type(adaptive_stepper),  intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      type(response),       intent(inout) :: resp
      type(sle_solver),         intent(inout) :: sle
      real(wp),                 intent(in)    :: topo0(:,:)
      real(wp),                 intent(in)    :: ice0(:,:), ice1(:,:)  !! abs. ice at t0,t1
      real(wp),                 intent(in)    :: ice_ref(:,:)          !! reference ice
      real(wp),                 intent(in)    :: t0, t1
      real(wp),                 intent(inout) :: rsl(:,:)
      real(wp),                 intent(out)   :: C(:,:)
      !! rot (optional): rotational feedback (vilma_rotation), stepped with the solid
      !! Earth: every SLE solve iterates it to the fixed point at the end of its own
      !! sub-step, and it is saved / restored wherever the response is.
      type(rotation_state), optional, intent(inout) :: rot
      !! sigma_out (optional): the SLE's converged spectral surface mass load [kg m⁻²]
      !! at t1 — the SAME load the response saw, including the subgrid sloping-coast
      !! term. Returned so the caller drives the end-of-interval rotational feedback
      !! from it rather than re-deriving the load (vilma_coupling). The final SLE solve of
      !! the run lands on the accepted state at t1, so its load is the t1 load.
      complex(wp),    optional, intent(out)   :: sigma_out(:)
      !! u_out, N_out, dphi_out (optional): the converged solid uplift u and geoid
      !! rise N (spectral, nlm) and the uniform sea-surface shift Δφ at t1, from the
      !! same final SLE solve as rsl, so rsl = N − u + Δφ (on the grid, to round-off).
      complex(wp),    optional, intent(out)   :: u_out(:), N_out(:)
      real(wp),       optional, intent(out)   :: dphi_out
      !! host_ref, host0, host1 (optional, all or none): host-resolved terms (vilma_sle,
      !! sle_host) at the reference and at t0 and t1. They ramp linearly over the interval,
      !! field by field and with the same formula as the ice, so every sub-step solves
      !! with host(t) = host0 + (t−t0)/(t1−t0)·(host1−host0). The ice then is the own
      !! part's (1−w). The σ_n priming uses host0; the dt = 0 seed uses host1 with ice0
      !! (like the ice, see KNOWN ISSUE below). Without them the stepper is the host-free one.
      type(sle_host), optional, intent(in)    :: host_ref, host0, host1
      !! dIg_out (optional): the grounded-ice increment ΔI_g [m] of the final SLE solve (t1).
      real(wp),       optional, intent(out)   :: dIg_out(:,:)
      !! fill_out (optional): the flooded-area fill of the final SLE solve (t1) integrated over
      !! the unit sphere [m sr]: (1) the whole fill, (2) the own part's share (sle_result).
      real(wp),       optional, intent(out)   :: fill_out(2)

      type(sle_result)      :: res
      type(sle_host)        :: host_now
      real(wp), allocatable :: rsl_n(:,:), ice_now(:,:), dice_now(:,:), dIg_last(:,:)
      complex(wp), allocatable :: sig0(:), sig_last(:), u_last(:), N_last(:)
      real(wp) :: dphi_last, fill_last(2)
      real(wp) :: t, dt, err_inf, tau_inf, errsc, fac, ricfac, expo, span
      real(wp) :: rate, dt_stab
      integer  :: p, np, nl, n_sub
      logical  :: at_floor, accept, hmode
      integer(kind=8) :: gc0, gc1, grate          ! PROFILE: see %t_guard

      span = t1 - t0
      hmode = present(host_ref)
      if (hmode .neqv. (present(host0) .and. present(host1))) &
         error stop 'stepper_advance: host_ref, host0 and host1 go together'
      if (hmode) then
         allocate(host_now%w(sht%nphi,sht%nlat), host_now%I_g(sht%nphi,sht%nlat), &
                  host_now%C(sht%nphi,sht%nlat), host_now%F(sht%nphi,sht%nlat))
      end if
      ! Capture the converged load on every SLE solve; after the run sig_last holds
      ! the final (t1) solve's load. Always allocated so solve_at can fill it cheaply.
      allocate(sig_last(sht%nlm), u_last(sht%nlm), N_last(sht%nlm), dIg_last(sht%nphi,sht%nlat))
      dIg_last = 0.0_wp
      dphi_last = 0.0_wp;  fill_last = 0.0_wp
      self%worst_mass_resid = 0.0_wp             ! per-interval diagnostics (reset each advance)
      self%n_sle = 0;  self%n_cap_inner = 0;  self%n_cap_coast = 0;  self%n_cap_mem = 0
      self%worst_resid = 0.0_wp;  self%worst_coast_flip = 0

      ! Zero-length interval: no time or memory advance, but still refresh the diagnostic
      ! fields (rsl, C) so they are consistent with the CURRENT memory + ice load. A
      ! report-only SLE solve converges the load against the frozen entering memory and
      ! returns without advancing it. The drive uses solid_earth_update(...,dt=0) to seed
      ! the entering ice and populate rsl/z_bed/C from a seeded or restored memory state;
      ! a cross-resolution restart restores only the spectral memory (not the 2-D rsl), so
      ! without this the caller's rsl would stay zero on the first output sample.
      !
      ! KNOWN ISSUE -- the diagnostic is built from ice0, the ice the model is
      ! already carrying, NOT ice1, the slice the caller just handed in. The caller
      ! (solid_earth_update) overwrites gg%h_ice with the new ice the moment this
      ! returns, so after a dt=0 call rsl/C describe the PREVIOUS ice while h_ice
      ! describes the new one, and update_bsl then combines the new h_ice with the
      ! old C. Three consequences:
      !   * solid_earth_update(se, X, 0.0) is not idempotent -- calling it twice
      !     with the same X gives different rsl/C.
      !   * the two backends differ at t0. vilma_drive seeds with dt=0 and
      !     solid_earth_init leaves h_ice = h_ice_eq, so the native seed reports
      !     the REFERENCE coastline while the VILMA1 path (which sets h_ice before
      !     its own diagnostics) reports the LGM one. They differ over every cell
      !     where LGM grounded ice sits on sub-sea-level bed -- Hudson Bay, the
      !     Barents shelf, West Antarctica -- in both C_ocean and the t0 bsl.
      !   * it affects only the t0 diagnostic sample, not the transient physics.
      ! Switching to ice1 is the obviously right expression and was tried; it makes
      ! solid_earth_spinup stop relaxing (relax_hold's first pass then measures
      ! ~4e-7 m of bed motion instead of 136 m and exits "converged" immediately,
      ! leaving 29 m of subsidence where the equilibrium is 239 m). That interaction
      ! is not understood, and initialization I1 -- which every production
      ! deglaciation depends on -- runs through this path, so the one-line change is
      ! NOT safe to make on its own. Fixing it means understanding why a seeded rsl
      ! suppresses the subsequent viscous advance; test_spinup is the reproducer.
      if (span <= 0.0_wp) then
         allocate(ice_now(sht%nphi, sht%nlat), dice_now(sht%nphi, sht%nlat))
         ice_now  = ice0
         dice_now = ice0 - ice_ref
         call solve(dice_now, ice_now, sig_last, 0.0_wp, host1, report_only=.true.)
         call tally()
         call hand_back()
         return
      end if
      if (.not. scheme_is_implicit(resp%scheme)) then
         ! --- explicit (forward-Euler) memory: a-priori stability sub-stepping -------
         ! The 1st-order memory carries no embedded error signal, so step-doubling is
         ! meaningless. Instead the interval is divided into equal sub-steps sized by
         ! the Maxwell stability ceiling Δt ≤ cfl/max(μ/η) (M = μΔt/η ≤ cfl), as VILMA1
         ! advances its explicit memory (small fixed Δt). max(μ/η) is exact
         ! (response_max_rate), so within the ceiling the memory cannot run away and no
         ! growth test is made; the one check is that the memory norm stays finite.
         np = sht%nphi;  nl = sht%nlat
         allocate(rsl_n(np,nl), ice_now(np,nl), dice_now(np,nl))
         rate = response_max_rate(resp)
         if (rate > 0.0_wp) then
            dt_stab = self%cfl/rate                  ! M = μΔt/η ≤ cfl
         else
            dt_stab = span                           ! purely elastic: one step suffices
         end if
         dt_stab = min(dt_stab, self%dt_max)
         if (self%dt_min > dt_stab) then
            ! the floor wins over stability: say so once, with the Maxwell number it implies
            if (.not. self%dt_min_overrides) write(*,'(a,es10.3,a,es10.3,a,f0.2,a)') &
               ' WARNING vilma_timestep: dt_min = ', self%dt_min/sec_per_year, &
               ' yr exceeds the explicit sub-step limit ', dt_stab/sec_per_year, &
               ' yr and overrides it: Maxwell number M = ', rate*self%dt_min, &
               ' (cfl asked for less; forward Euler is unstable for M > 2)'
            self%dt_min_overrides = .true.
            dt_stab = self%dt_min
         end if
         n_sub = max(1, ceiling(span/dt_stab - 1.0e-9_wp))
         dt    = span/real(n_sub, wp)                 ! nominal (equal) sub-step
         t = t0
         do
            if (t >= t1 - 1.0e-9_wp*span) exit
            dt   = min(dt, t1 - t)
            rsl_n = rsl
            call system_clock(gc0, grate)
            call response_save_state(resp)
            if (present(rot)) call rotation_save_state(rot)
            call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
            call response_set_dt(resp, dt)
            call solve_at(t + dt, dt)                 ! advances memory by dt
            call system_clock(gc0, grate)
            err_inf = response_memory_norm(resp)
            call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
            accept  = is_finite(err_inf)             ! .false. for NaN / Inf
            if (accept) then
               t = t + dt;  self%n_accept = self%n_accept + 1
               dt = min(span/real(n_sub, wp), 2.0_wp*dt)   ! recover toward nominal
            else
               call system_clock(gc0, grate)
               call response_restore_state(resp);  rsl = rsl_n
               if (present(rot)) call rotation_restore_state(rot)
               call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
               self%n_reject = self%n_reject + 1
               dt = 0.5_wp*dt
               if (dt <= 1.0e-12_wp*span) error stop &
                  'vilma_timestep: the memory is not finite at any explicit sub-step (corrupt input or state)'
            end if
         end do
         call hand_back()
         return
      end if

      np = sht%nphi;  nl = sht%nlat
      allocate(rsl_n(np,nl), ice_now(np,nl), dice_now(np,nl))

      ! Seed the trapezoidal start-of-step load σ_0 once, if not already tracked: at the
      ! interval start the response carries its entering memory, so a report-only SLE
      ! solve (no memory/time advance) gives the load consistent with it — at t=0 that
      ! is the elastic-consistent load. Without this the first step falls back to the
      ! σ_{n+1} proxy for σ_n (O(Δt) for the relaxing SLE load → 2nd order on that step).
      if (.not. resp%sigma_primed) then
         allocate(sig0(sht%nlm))
         ice_now = ice0;  dice_now = ice0 - ice_ref
         call solve(dice_now, ice_now, sig0, 0.0_wp, host0, report_only=.true.)
         call response_prime_sigma(resp, sig0)
      end if

      p      = scheme_order(resp%scheme)          ! global order (2 for trapezoidal)
      ricfac = real(2**p - 1, wp)                 ! Richardson factor for the estimate
      expo   = 1.0_wp/real(p + 1, wp)             ! step-size exponent (local err ∝ Δt^{p+1})

      t = t0
      if (self%dt_try <= 0.0_wp) self%dt_try = span
      do
         if (t >= t1 - 1.0e-9_wp*span) exit
         dt = min(self%dt_try, t1 - t)
         at_floor = (self%dt_min > 0.0_wp .and. dt <= self%dt_min*(1.0_wp + 1.0e-9_wp))

         ! --- field step-doubling around [t, t+dt] -------------------------------
         rsl_n = rsl
         call system_clock(gc0, grate)
         call response_save_state(resp)                   ! buffer A = τ_n
         if (present(rot)) call rotation_save_state(rot)
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         call response_set_dt(resp, dt)
         call solve_at(t + dt, dt)                 ! coarse: one Δt
         call system_clock(gc0, grate)
         call response_stash_coarse(resp)                  ! buffer B = τ_coarse
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         call system_clock(gc0, grate)
         call response_restore_state(resp);  rsl = rsl_n   ! back to τ_n (and its rsl seed)
         if (present(rot)) call rotation_restore_state(rot) ! the rotation follows the fine pair
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         call response_set_dt(resp, 0.5_wp*dt)
         call solve_at(t + 0.5_wp*dt, 0.5_wp*dt)   ! fine sub-step 1
         call solve_at(t + dt, 0.5_wp*dt)          ! fine sub-step 2 → τ_fine
         call system_clock(gc0, grate)
         call response_coarse_fine_error(resp, err_inf, tau_inf)
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         call response_set_dt(resp, dt)                       ! leave resp%dt at the step size
         errsc = (err_inf/ricfac) / (self%atol + self%rtol*tau_inf)
         ! Two guards. Without them a non-finite estimate makes fac NaN and dt_try NaN,
         ! so t never advances; and with the default dt_min = 0 a persistently rejected
         ! step shrinks geometrically toward zero forever. Both are an infinite hang
         ! rather than an error.
         if (.not. is_finite(errsc)) error stop &
            'vilma_timestep: implicit local-error estimate went non-finite &
            &(viscosity too stiff, or the state has already diverged)'

         ! --- accept / reject ----------------------------------------------------
         accept = (errsc <= 1.0_wp) .or. at_floor
         if (accept) then
            t = t + dt;  self%n_accept = self%n_accept + 1
            if (at_floor .and. errsc > 1.0_wp) self%n_floor = self%n_floor + 1
         else
            call system_clock(gc0, grate)
            call response_restore_state(resp);  rsl = rsl_n
            if (present(rot)) call rotation_restore_state(rot)
            call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
            self%n_reject = self%n_reject + 1
         end if

         ! --- step-size update ---------------------------------------------------
         if (errsc <= 0.0_wp) then
            fac = self%grow_max
         else
            fac = self%safety * errsc**(-expo)
         end if
         fac = min(self%grow_max, max(self%shrink_min, fac))
         self%dt_try = min(self%dt_max, max(self%dt_min, dt*fac))
         if (self%dt_try <= 1.0e-12_wp*span) error stop &
            'vilma_timestep: implicit adaptive step collapsed (tolerance unreachable)'
      end do

      ! sig_last now holds the final accepted step's fine sub-step at t1.
      call hand_back()

   contains

      subroutine hand_back()
         !! Return the optional outputs of the final SLE solve (the one at t1).
         if (present(sigma_out)) sigma_out = sig_last
         if (present(u_out))     u_out     = u_last
         if (present(N_out))     N_out     = N_last
         if (present(dphi_out))  dphi_out  = dphi_last
         if (present(dIg_out))   dIg_out   = dIg_last
         if (present(fill_out))  fill_out  = fill_last
      end subroutine hand_back

      subroutine tally()
         !! Record the solve just made (res): the fields handed back from the final one,
         !! and the interval's convergence record.
         self%worst_mass_resid = max(self%worst_mass_resid, res%mass_resid)
         u_last = res%u_lm;  N_last = res%N_lm;  dphi_last = res%esl
         fill_last = [res%fill_int, res%fill_own]
         self%n_sle = self%n_sle + 1
         if (res%inner_capped) self%n_cap_inner = self%n_cap_inner + 1
         if (res%coast_capped) self%n_cap_coast = self%n_cap_coast + 1
         if (res%mem_capped)   self%n_cap_mem   = self%n_cap_mem + 1
         self%worst_resid      = max(self%worst_resid, res%resid_rel)
         self%worst_coast_flip = max(self%worst_coast_flip, res%n_coast_flip)
      end subroutine tally

      subroutine solve(dice, ice, sigma, h, hn, report_only)
         !! The one place sle_solve is called: with the host terms hn (and the reference)
         !! when the stepper has them, else exactly the host-free call. Solves into the
         !! stepper's rsl, C and res; sigma is the spectral load, h the rotation step.
         real(wp),    intent(in)  :: dice(:,:), ice(:,:)
         complex(wp), intent(out) :: sigma(:)
         real(wp),    intent(in)  :: h
         type(sle_host), optional, intent(in) :: hn
         logical,        optional, intent(in) :: report_only
         if (hmode) then
            call sle_solve(sle, sht, resp, dice, ice, topo0, rsl, C, res, &
                           report_only=report_only, sigma_lm=sigma, rot=rot, rot_dt=h, &
                           host=hn, host_ref=host_ref, dIg_out=dIg_last)
         else
            call sle_solve(sle, sht, resp, dice, ice, topo0, rsl, C, res, &
                           report_only=report_only, sigma_lm=sigma, rot=rot, rot_dt=h, &
                           dIg_out=dIg_last)
         end if
      end subroutine solve

      subroutine solve_at(t_eval, h)
         !! One SLE solve (co-converged memory advance) with the ice load interpolated
         !! at t_eval. Advances resp by the current resp%dt (set by the caller) and the
         !! rotation by h, the same sub-step (passed: resp%dt is not set for elastic).
         real(wp), intent(in) :: t_eval, h
         real(wp) :: frac
         frac = (t_eval - t0)/span
         ice_now  = ice0 + frac*(ice1 - ice0)
         dice_now = ice_now - ice_ref
         self%n_solve = self%n_solve + 1
         if (hmode) then
            host_now%w   = host0%w   + frac*(host1%w   - host0%w)
            host_now%I_g = host0%I_g + frac*(host1%I_g - host0%I_g)
            host_now%C   = host0%C   + frac*(host1%C   - host0%C)
            host_now%F   = host0%F   + frac*(host1%F   - host0%F)
            call solve(dice_now, ice_now, sig_last, h, host_now)
         else
            call solve(dice_now, ice_now, sig_last, h)
         end if
         call tally()
      end subroutine solve_at

   end subroutine stepper_advance

end module vilma_timestep
