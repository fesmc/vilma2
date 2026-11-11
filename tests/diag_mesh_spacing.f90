program diag_mesh_spacing
   !! Radial-mesh resolution study (not a pass/fail test): what the depth-based
   !! element size of radial_mesh_build buys, and how far the default sizes are
   !! from converged. Usage:
   !!
   !!   diag_mesh_spacing.x [lvz3d]
   !!
   !! (A) element counts of the default structures under the earlier rule (size
   !!     from the LAYER's mid-depth) and the present one (size from the
   !!     ELEMENT's depth).
   !! (B) M3-L70-V01 elastic and fluid loading Love numbers against the
   !!     normal-mode table, for mesh_dr_scale = 2, 1, 0.5, 0.25.
   !! (C) a 100 km ice disc (0 -> 100 m over 100 yr, then held; the load of
   !!     test_benchmark_lvz) on M3-L70-V01 and on two coarse user layerings,
   !!     each on the earlier rule's mesh, the present one at several scales,
   !!     and a 0.25x reference:
   !!       - "to670": one Maxwell layer surface -> 670 km, lower mantle, core;
   !!       - "toCMB": one Maxwell layer surface -> CMB, core; alone, and with a
   !!         global 100 km low-viscosity layer (1e19 Pa s, 100-200 km depth)
   !!         given at the NODES and averaged onto elements, as a 3-D viscosity
   !!         file is. Its reference carries the layer as true layers.
   !!     and on PREM with its earlier (220) and present (219 element) mesh.
   !! (D) with the argument "lvz3d": the Weerdesteijn (2023) cylinder LVZ of
   !!     test_benchmark_lvz at mesh_dr_scale = 2, 1, 0.5 (minutes, threaded).
   use vilma_precision,       only: wp
   use vilma_constants,       only: kyr
   use vilma_params,          only: vilma_param_class
   use vilma_earth_structure, only: earth_model, earth_layer, build_M3L70V01, build_earth, &
                                    RHEOL_MAXWELL, RHEOL_FLUID
   use vilma_radial_fe,       only: radial_operator_solve, radial_operator_assemble, radial_mesh_build, &
                                    radial_mesh, radial_operator, loading_love, radial_fe_finalize
   use vilma_response,        only: response, response_init_ve, response_apply, response_begin_step, &
                                    response_commit_step, response_enable_lateral_visc, response_destroy
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_eval_point, sht_grid_lmidx
   implicit none

   character(*), parameter :: REF = 'data/benchmarks/love_M3-L70-V01/mod_M3-L70-V01'
   integer,  parameter :: NMAX = 256, LMAX = 512, NLATF = 3
   real(wp), parameter :: km = 1.0e3_wp, YR = 0.001_wp*kyr, A = 6371.0_wp*km
   real(wp), parameter :: DT = 2.5_wp*YR, R_ICE = 100.0_wp*km, RHO_ICE = 931.0_wp
   real(wp), parameter :: H_MAX = 100.0_wp, T_RAMP = 100.0_wp*YR
   real(wp), parameter :: DR_DEF(3) = [5.0_wp, 10.0_wp, 40.0_wp]*km, DEPTH_DEF(2) = [70.0_wp, 670.0_wp]*km
   real(wp), parameter :: NO_DEPTH(0) = [real(wp) ::]
   real(wp), parameter :: LVZ_PERT = -2.0_wp
   integer,  parameter :: LVZ_NONE = 0, LVZ_NODES = 1, LVZ_CYL = 2
   real(wp), parameter :: PREM_TOP(11) = [6371.0_wp, 6291.0_wp, 6151.0_wp, 5971.0_wp, &
        5771.0_wp, 5701.0_wp, 5600.0_wp, 4943.0_wp, 4287.0_wp, 3630.0_wp, 3480.0_wp]*km

   type(sht_grid)    :: sht
   type(earth_model) :: e
   real(wp) :: he(NMAX), le(NMAX), ke(NMAX), hf(NMAX), lf(NMAX), kf(NMAX)
   real(wp) :: uref(2), u(2)
   real(wp), allocatable :: disc_lm(:)
   character(16) :: arg
   logical  :: okread
   integer  :: is
   real(wp), parameter :: SCALES(4) = [2.0_wp, 1.0_wp, 0.5_wp, 0.25_wp]

   call get_command_argument(1, arg)
   call sht_grid_init(sht, LMAX, nlat=NLATF*LMAX, nphi=2, mmax=0)
   allocate(disc_lm(sht%nlm))
   call cap_coeffs(R_ICE/A, disc_lm)

   if (trim(arg) == 'lvz3d') then
      write(*,'(a)') ' (D) Weerdesteijn 2023 cylinder LVZ on M3-L70-V01 (ref -0.75 m, -1.23 m at 200 yr)'
      write(*,'(a)') '     scale    ne   U_homog[m]   U_lvz[m]'
      do is = 1, 3
         e = build_M3L70V01();  e%mesh_dr = SCALES(is)*DR_DEF;  e%mesh_depth = DEPTH_DEF
         call disc_run(e, LVZ_NONE, 70.0_wp*km, 170.0_wp*km, 200.0_wp*YR, uref)
         call disc_run(e, LVZ_CYL,  70.0_wp*km, 170.0_wp*km, 200.0_wp*YR, u)
         write(*,'(f10.2,i6,2f12.4)') SCALES(is), ne_of(e), uref(1), u(1)
      end do
      call radial_fe_finalize()
      stop
   end if

   ! --- (A) element counts ------------------------------------------------------
   write(*,'(a)') ' (A) elements: rule by layer mid-depth -> rule by element depth'
   e = build_M3L70V01()
   write(*,'(a,i5,a,i5)') '     M3-L70-V01      ', ne_layer_rule(e), ' ->', ne_of(e)
   e = stack(PREM_TOP, [(1.0e21_wp, is = 1, 10)])
   write(*,'(a,i5,a,i5)') '     PREM, 11 layers ', ne_layer_rule(e), ' ->', ne_of(e)
   e = stack([A, A - 670.0_wp*km, 3480.0_wp*km], [1.0e21_wp, 2.0e21_wp])
   write(*,'(a,i5,a,i5)') '     to670           ', ne_layer_rule(e), ' ->', ne_of(e)
   e = stack([A, 3480.0_wp*km], [1.0e21_wp])
   write(*,'(a,i5,a,i5)') '     toCMB           ', ne_layer_rule(e), ' ->', ne_of(e)

   ! --- (B) Love numbers vs the normal-mode table -------------------------------
   call read_ref(REF, he, le, ke, hf, lf, kf, okread)
   if (.not. okread) error stop 'diag_mesh_spacing: cannot read the Love table (run from the repo root)'
   write(*,'(/,a)') ' (B) M3-L70-V01 Love numbers vs the normal-mode table: max |rel. error| [%]'
   write(*,'(a)') '     scale    ne | elastic h, l, k: j<=48          j<=256             | fluid h, l, k: j<=8'
   do is = 1, 4
      call love_errors(SCALES(is))
   end do

   ! --- (C) disc load -----------------------------------------------------------
   write(*,'(/,a)') ' (C) 100 km disc, U at the centre [m] at 200 yr and 1000 yr; error against the reference [%]'
   write(*,'(a)') '     case / mesh                          ne   U(200)    U(1000)   err200  err1000'

   write(*,'(a)') '   M3-L70-V01'
   e = build_M3L70V01()
   call ref_run(e, LVZ_NONE)
   do is = 1, 3
      call scaled_run(e, 'element depth', SCALES(is), LVZ_NONE)
   end do

   write(*,'(a)') '   to670: Maxwell 0-670 km in one layer'
   e = stack([A, A - 670.0_wp*km, 3480.0_wp*km], [1.0e21_wp, 2.0e21_wp])
   call ref_run(e, LVZ_NONE)
   call table_run(e, 'layer mid-depth (10 km | 40 km)', [10.0_wp, 40.0_wp]*km, [670.0_wp*km], LVZ_NONE)
   do is = 1, 3
      call scaled_run(e, 'element depth', SCALES(is), LVZ_NONE)
   end do

   write(*,'(a)') '   toCMB: Maxwell 0-2891 km in one layer'
   e = stack([A, 3480.0_wp*km], [1.0e21_wp])
   call ref_run(e, LVZ_NONE)
   call table_run(e, 'layer mid-depth (40 km)', [40.0_wp*km], NO_DEPTH, LVZ_NONE)
   do is = 1, 3
      call scaled_run(e, 'element depth', SCALES(is), LVZ_NONE)
   end do

   write(*,'(a)') '   toCMB + global low-viscosity layer 100-200 km (1e19 Pa s) given at the nodes'
   e = stack([A, A - 100.0_wp*km, A - 200.0_wp*km, 3480.0_wp*km], [1.0e21_wp, 1.0e19_wp, 1.0e21_wp])
   call ref_run(e, LVZ_NONE)                     ! the layer as true layers, 0.25x
   e = stack([A, 3480.0_wp*km], [1.0e21_wp])
   call table_run(e, 'layer mid-depth (40 km)', [40.0_wp*km], NO_DEPTH, LVZ_NODES)
   do is = 1, 4
      call scaled_run(e, 'element depth', SCALES(is), LVZ_NODES)
   end do

   write(*,'(a)') '   PREM, 11 layers (Clemens 1-D viscosity)'
   e = prem()
   call ref_run(e, LVZ_NONE)
   call table_run(e, 'earlier mesh (zones at 80, 670 km)', DR_DEF, [80.0_wp, 670.0_wp]*km, LVZ_NONE)
   do is = 1, 3
      call scaled_run(e, 'element depth', SCALES(is), LVZ_NONE)
   end do
   call radial_fe_finalize()

contains

   function stack(r_top, eta) result(em)
      !! Maxwell layers with the given tops (surface first) and viscosities over a
      !! fluid core; one mantle material throughout, so the layering changes
      !! nothing but the mesh and the viscosity.
      real(wp), intent(in) :: r_top(:), eta(:)
      type(earth_model) :: em
      integer :: k, n
      n = size(r_top)
      em%name = "stack";  em%r_earth = r_top(1);  em%r_core = r_top(n)
      allocate(em%layers(n))
      do k = 1, n - 1
         em%layers(k) = earth_layer(r_top(k+1), r_top(k), 4450.0_wp, 7.0363e10_wp, eta(k), RHEOL_MAXWELL)
      end do
      em%layers(n) = earth_layer(0.0_wp, r_top(n), 10750.0_wp, 0.0_wp, 0.0_wp, RHEOL_FLUID)
   end function stack

   function prem() result(em)
      !! The PREM block of input/vilma_defaults.nml.
      type(earth_model) :: em
      type(vilma_param_class) :: p
      p%earth = "PREM";  p%n_layer = 11
      p%r_top(1:11) = PREM_TOP
      p%r_bot(1:10) = PREM_TOP(2:11);  p%r_bot(11) = 0.0_wp
      p%eta(1:11)   = [1.0e40_wp, 4.0e20_wp, 4.0e20_wp, 4.0e20_wp, 4.0e20_wp, &
                       1.0e22_wp, 1.0e22_wp, 1.0e22_wp, 1.0e22_wp, 1.0e22_wp, 0.0_wp]
      p%rheology(1:11) = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]
      em = build_earth(p)
   end function prem

   integer function ne_of(em) result(ne)
      type(earth_model), intent(in) :: em
      type(radial_mesh) :: m
      call radial_mesh_build(m, em)
      ne = m%ne
   end function ne_of

   integer function ne_layer_rule(em) result(ne)
      !! Element count of the earlier rule: uniform elements per layer, sized by
      !! the target at the layer's mid-depth.
      type(earth_model), intent(in) :: em
      integer  :: k
      real(wp) :: depth_mid, dr
      ne = 0
      do k = 1, size(em%layers)
         depth_mid = em%r_earth - 0.5_wp*(em%layers(k)%r_bot + em%layers(k)%r_top)
         dr = DR_DEF(3)
         if (depth_mid <= DEPTH_DEF(2)) dr = DR_DEF(2)
         if (depth_mid <= DEPTH_DEF(1)) dr = DR_DEF(1)
         ne = ne + max(1, ceiling((em%layers(k)%r_top - em%layers(k)%r_bot)/dr))
      end do
   end function ne_layer_rule

   subroutine love_errors(scale)
      real(wp), intent(in) :: scale
      type(earth_model)     :: em, ef
      type(radial_mesh)     :: m
      type(radial_operator) :: op
      real(wp) :: uu, vv, ff, h, l, k, ee(3,2), ef8(3)
      integer  :: j
      em = build_M3L70V01();  em%mesh_dr = scale*DR_DEF;  em%mesh_depth = DEPTH_DEF
      ef = em
      ef%layers(2:4)%mu = 0.0_wp;  ef%layers(2:4)%rheology = RHEOL_FLUID
      call radial_mesh_build(m, em)
      ee = 0.0_wp;  ef8 = 0.0_wp
      do j = 2, NMAX
         call radial_operator_assemble(op, em, m, j)
         call radial_operator_solve(op, 1.0_wp, uu, vv, ff)
         call loading_love(em, j, 1.0_wp, uu, vv, ff, h, l, k)
         if (j <= 48) ee(:,1) = max(ee(:,1), [reldiff(h,he(j)), reldiff(l,le(j)), reldiff(k,ke(j))])
         ee(:,2) = max(ee(:,2), [reldiff(h,he(j)), reldiff(l,le(j)), reldiff(k,ke(j))])
         if (j <= 8) then
            call radial_operator_assemble(op, ef, m, j)
            call radial_operator_solve(op, 1.0_wp, uu, vv, ff)
            call loading_love(ef, j, 1.0_wp, uu, vv, ff, h, l, k)
            ef8 = max(ef8, [reldiff(h,hf(j)), reldiff(l,lf(j)), reldiff(k,kf(j))])
         end if
      end do
      write(*,'(f10.2,i6,a,3f8.4,2x,3f8.4,a,3f8.4)') scale, m%ne, ' |', 100.0_wp*ee(:,1), &
           100.0_wp*ee(:,2), ' |', 100.0_wp*ef8
   end subroutine love_errors

   subroutine ref_run(em, lvz)
      !! The 0.25x mesh of the element-depth rule: sets uref.
      type(earth_model), intent(in) :: em
      integer,           intent(in) :: lvz
      type(earth_model) :: er
      er = em;  er%mesh_dr = 0.25_wp*DR_DEF;  er%mesh_depth = DEPTH_DEF
      call disc_run(er, lvz, 100.0_wp*km, 200.0_wp*km, 1000.0_wp*YR, uref)
      write(*,'(5x,a36,i5,2f10.4)') 'reference (element depth, 0.25x)    ', ne_of(er), uref
   end subroutine ref_run

   subroutine scaled_run(em, label, scale, lvz)
      type(earth_model), intent(in) :: em
      character(*),      intent(in) :: label
      real(wp),          intent(in) :: scale
      integer,           intent(in) :: lvz
      character(36) :: lab
      write(lab,'(a,a,f4.2,a)') label, ', ', scale, 'x'
      call table_run(em, lab, scale*DR_DEF, DEPTH_DEF, lvz)
   end subroutine scaled_run

   subroutine table_run(em, label, dr, depth, lvz)
      type(earth_model), intent(in) :: em
      character(*),      intent(in) :: label
      real(wp),          intent(in) :: dr(:), depth(:)
      integer,           intent(in) :: lvz
      type(earth_model) :: er
      character(36) :: lab
      real(wp) :: uu(2)
      er = em;  er%mesh_dr = dr;  er%mesh_depth = depth
      call disc_run(er, lvz, 100.0_wp*km, 200.0_wp*km, 1000.0_wp*YR, uu)
      lab = label
      write(*,'(5x,a36,i5,2f10.4,2f9.3)') lab, ne_of(er), uu, 100.0_wp*(uu - uref)/abs(uref)
   end subroutine table_run

   subroutine disc_run(em, lvz, lvz_top, lvz_bot, t_end, uc)
      !! The disc load of test_benchmark_lvz on earth `em`; uc = U at the load
      !! centre at 200 yr and at t_end. lvz: LVZ_NODES, a global layer between the
      !! depths lvz_top and lvz_bot set at the nodes (half weight on a node at
      !! either bound) and averaged onto elements; LVZ_CYL, the benchmark's
      !! cylinder under the load, set by element midpoint.
      type(earth_model), intent(in)  :: em
      integer,           intent(in)  :: lvz
      real(wp),          intent(in)  :: lvz_top, lvz_bot, t_end
      real(wp),          intent(out) :: uc(2)
      type(response) :: ve
      complex(wp), allocatable :: slm(:), ulm(:), nlm(:)
      real(wp),    allocatable :: pert(:,:,:), pnode(:)
      real(wp) :: t, depth, uval
      integer  :: i, ie, j, nstep

      call response_init_ve(ve, em, sht, DT)
      allocate(slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))
      if (lvz /= LVZ_NONE) then
         allocate(pert(sht%nphi, sht%nlat, ve%ne), pnode(ve%nr));  pert = 0.0_wp
         do i = 1, ve%nr
            depth = em%r_earth - ve%r(i)
            pnode(i) = 0.0_wp
            if (depth > lvz_top - 1.0_wp .and. depth < lvz_bot + 1.0_wp) pnode(i) = 0.5_wp*LVZ_PERT
            if (depth > lvz_top + 1.0_wp .and. depth < lvz_bot - 1.0_wp) pnode(i) = LVZ_PERT
         end do
         do ie = 1, ve%ne
            if (lvz == LVZ_NODES) then
               pert(:,:,ie) = 0.5_wp*(pnode(ie) + pnode(ie+1))
            else
               depth = em%r_earth - 0.5_wp*(ve%r(ie) + ve%r(ie+1))
               if (depth > lvz_top .and. depth < lvz_bot) then
                  do j = 1, sht%nlat
                     if (sht%colat(j) <= R_ICE/em%r_earth) pert(:,j,ie) = LVZ_PERT
                  end do
               end if
            end if
         end do
         call response_enable_lateral_visc(ve, sht, pert)
      end if
      nstep = nint(t_end/DT)
      do i = 1, nstep
         t = real(i,wp)*DT
         slm = RHO_ICE*min(t/T_RAMP, 1.0_wp)*H_MAX*disc_lm
         call response_begin_step(ve, sht)
         call response_apply(ve, sht, slm, ulm, nlm)
         call response_commit_step(ve, sht, slm)
         if (i == nint(200.0_wp*YR/DT) .or. i == nstep) then
            call sht_grid_eval_point(sht, ulm, 0.0_wp, 0.0_wp, uval)
            if (i == nint(200.0_wp*YR/DT)) uc(1) = uval
            if (i == nstep) uc(2) = uval
         end if
      end do
      call response_destroy(ve)
   end subroutine disc_run

   subroutine cap_coeffs(tc, lm_coeffs)
      !! Exact m=0 coefficients (SHTns orthonormal) of the unit cap indicator, as
      !! in test_benchmark_lvz.
      real(wp),    intent(in)  :: tc
      complex(wp), intent(out) :: lm_coeffs(:)
      real(wp) :: c, p(0:LMAX+1)
      real(wp), parameter :: sqpi = 1.7724538509055160_wp
      integer  :: l
      c = cos(tc)
      p(0) = 1.0_wp;  p(1) = c
      do l = 1, LMAX
         p(l+1) = (real(2*l+1,wp)*c*p(l) - real(l,wp)*p(l-1))/real(l+1,wp)
      end do
      lm_coeffs = (0.0_wp, 0.0_wp)
      lm_coeffs(sht_grid_lmidx(sht,0,0)) = cmplx(sqpi*(1.0_wp - c), 0.0_wp, wp)
      do l = 1, LMAX
         lm_coeffs(sht_grid_lmidx(sht,l,0)) = cmplx( &
            sqpi/sqrt(real(2*l+1,wp))*(p(l-1) - p(l+1)), 0.0_wp, wp)
      end do
   end subroutine cap_coeffs

   pure real(wp) function reldiff(a, b) result(d)
      real(wp), intent(in) :: a, b
      d = 0.0_wp
      if (abs(b) > 0.0_wp) d = abs(a-b)/abs(b)
   end function reldiff

   subroutine read_ref(fname, he, le, ke, hf, lf, kf, okread)
      !! The normal-mode table, as test_benchmark_love reads it.
      character(len=*), intent(in)  :: fname
      real(wp),         intent(out) :: he(:), le(:), ke(:), hf(:), lf(:), kf(:)
      logical,          intent(out) :: okread
      integer  :: u, i, ni, nm, mm, ios
      real(wp) :: kv, hv, lv, dum
      okread = .false.
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do i = 1, 5
         read(u,*,iostat=ios)
         if (ios /= 0) then;  close(u);  return;  end if
      end do
      do
         read(u,*,iostat=ios) ni, nm, kv, hv, lv
         if (ios /= 0) exit
         if (ni < 1 .or. ni > size(he)) exit
         ke(ni) = kv;  he(ni) = hv;  le(ni) = lv
         do mm = 1, nm
            read(u,*,iostat=ios) dum
            if (ios /= 0) then;  close(u);  return;  end if
         end do
         read(u,*,iostat=ios) ni, nm, kv, hv, lv
         if (ios /= 0) then;  close(u);  return;  end if
         kf(ni) = kv;  hf(ni) = hv;  lf(ni) = lv
      end do
      close(u)
      okread = .true.
   end subroutine read_ref

end program diag_mesh_spacing
