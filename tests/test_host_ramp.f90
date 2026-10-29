program test_host_ramp
   !! Host-resolved terms in the adaptive stepper (stepper_advance with host_ref, host0,
   !! host1): they ramp linearly over the interval like the ice.
   !!
   !! Checks:
   !!   (1) elastic response, one step: the step is one solve at t1, so it equals
   !!       sle_solve with host1 and ice1;
   !!   (2) dt = 0 seed: equals a report-only sle_solve with host1 and ice0, and not the
   !!       one with host0;
   !!   (3) ramp: one advance over four sub-steps (forward Euler, Maxwell sphere) equals four
   !!       advances over the quarter intervals with the ice and host terms interpolated to
   !!       their ends, and differs from holding the host terms at host0;
   !!   (4) trapezoidal scheme: the σ_n priming is the report-only solve with host0 (and
   !!       ice0): a run that primes itself equals one primed by hand with host0, and
   !!       differs from one primed with host1;
   !!   (5) dIg_out of the final solve is returned.
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, sec_per_year
   use vilma_earth_structure, only: earth_model, earth_layer, build_M3L70V01, RHEOL_MAXWELL
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response_destroy, response_prime_sigma, response_memory_norm, response, &
                                    response_init_elastic, response_init_ve
   use vilma_viscoelastic,    only: SCHEME_TRAP
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_sle,             only: sle_solve, sle_solver, sle_result, sle_host
   use vilma_timestep,        only: stepper_advance, adaptive_stepper
   implicit none

   integer,  parameter :: LMAX = 16
   real(wp), parameter :: DEG = pi/180.0_wp, yr = sec_per_year
   real(wp), parameter :: km = 1.0e3_wp, a_e = 6371.0_wp*km, rho_e = 5511.0_wp, mu_e = 1.4e11_wp, eta_e = 1.0e21_wp
   type(sht_grid)      :: sht
   type(earth_model)   :: e_ve, e_el
   real(wp), allocatable :: topo0(:,:), ice1(:,:), zero(:,:)
   type(sle_host)      :: hr, h0, h1
   logical             :: ok
   integer             :: np, nl

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX+2, nphi=4*LMAX)
   np = sht%nphi;  nl = sht%nlat
   e_el = build_M3L70V01()
   e_ve%name = "maxwell";  e_ve%r_earth = a_e;  e_ve%r_core = 0.0_wp
   allocate(e_ve%layers(1))
   e_ve%layers(1) = earth_layer(0.0_wp, a_e, rho_e, mu_e, eta_e, RHEOL_MAXWELL)
   allocate(topo0(np,nl), ice1(np,nl), zero(np,nl))
   zero = 0.0_wp
   call make_fields()

   call check_elastic_step()
   call check_dt0()
   call check_ramp()
   call check_priming()

   write(*,'(a)') ''
   call sht_grid_destroy(sht)
   if (ok) then
      write(*,'(a)') ' PASS: host terms ramp with the ice, seed at host1 and prime at host0'
      call radial_fe_finalize()
   else
      write(*,'(a)') ' FAIL: host terms in the stepper'
      call radial_fe_finalize()
      error stop 1
   end if

contains

   subroutine make_fields()
      !! Bed: land cap (colat < 50, +500 m), shelf ring (50..62, -100 m), bare ring (62..70,
      !! -60 m), deep ocean. Own ice at t1: 2 km cap, 300 m on the shelf; at t0 none. Host
      !! region: w = 1 (lon < 90), 1/2 (lon < 180) in colat 20..70; its ice at t1 is 3 km on
      !! the cap and 200 m beyond, its coast at colat 54; at t0 half the weight, a quarter
      !! of the ice, its coast at 58 (the reference's).
      integer  :: i, j
      real(wp) :: th, lon
      call alloc(hr, .false.);  call alloc(h0, .true.);  call alloc(h1, .true.)
      do j = 1, nl
         th = sht%colat(j)/DEG
         do i = 1, np
            lon = sht%lon(i)/DEG
            if (th < 50.0_wp) then
               topo0(i,j) = 500.0_wp
            else if (th < 62.0_wp) then
               topo0(i,j) = -100.0_wp
            else if (th < 70.0_wp) then
               topo0(i,j) = -60.0_wp
            else
               topo0(i,j) = -4000.0_wp
            end if
            ice1(i,j) = 0.0_wp
            if (th < 30.0_wp) ice1(i,j) = 2000.0_wp
            if (th >= 50.0_wp .and. th < 62.0_wp) ice1(i,j) = 300.0_wp
            if (th >= 20.0_wp .and. th < 70.0_wp .and. lon < 180.0_wp) then
               hr%w(i,j) = merge(1.0_wp, 0.5_wp, lon < 90.0_wp)
               call terms(h1, i, j, hr%w(i,j),        merge(3000.0_wp, 200.0_wp, th < 30.0_wp), th < 54.0_wp, 0.0_wp)
               call terms(h0, i, j, 0.5_wp*hr%w(i,j), merge(750.0_wp,   50.0_wp,  th < 30.0_wp), th < 58.0_wp, 0.0_wp)
               ! the reference: w, ice-free, coast at 58
               hr%C(i,j) = hr%w(i,j)*merge(0.0_wp, 1.0_wp, th < 58.0_wp)
               ! the t1 fill: from the reference coast to the t1 one, over its own bed
               h1%F(i,j) = merge(200.0_wp, -50.0_wp, th < 54.0_wp)*(h1%C(i,j) - hr%C(i,j))
            end if
         end do
      end do
   end subroutine make_fields

   subroutine alloc(h, with_f)
      type(sle_host), intent(inout) :: h
      logical,        intent(in)    :: with_f
      allocate(h%w(np,nl), h%I_g(np,nl), h%C(np,nl), source=0.0_wp)
      if (with_f) allocate(h%F(np,nl), source=0.0_wp)
   end subroutine alloc

   subroutine terms(h, i, j, w, H_ice, land, F)
      type(sle_host), intent(inout) :: h
      integer,        intent(in)    :: i, j
      real(wp),       intent(in)    :: w, H_ice, F
      logical,        intent(in)    :: land
      h%w(i,j) = w
      h%I_g(i,j) = merge(w*H_ice, 0.0_wp, land)
      h%C(i,j) = merge(0.0_wp, w, land)
      h%F(i,j) = F
   end subroutine terms

   subroutine lerp(ha, hb, f, hc)
      type(sle_host), intent(in)    :: ha, hb
      real(wp),       intent(in)    :: f
      type(sle_host), intent(inout) :: hc
      if (.not. allocated(hc%w)) call alloc(hc, .true.)
      hc%w   = ha%w   + f*(hb%w   - ha%w)
      hc%I_g = ha%I_g + f*(hb%I_g - ha%I_g)
      hc%C   = ha%C   + f*(hb%C   - ha%C)
      hc%F   = ha%F   + f*(hb%F   - ha%F)
   end subroutine lerp

   subroutine setup_ve(resp, sle, st, trap)
      type(response),        intent(out) :: resp
      type(sle_solver),      intent(out) :: sle
      type(adaptive_stepper), intent(out) :: st
      logical,               intent(in)  :: trap
      call response_init_ve(resp, e_ve, sht, 100.0_wp*yr)
      if (trap) then
         resp%scheme = SCHEME_TRAP;  resp%couple_tol = 1.0e-9_wp
      end if
      sle%n_outer = 12;  sle%n_inner = 80;  sle%tol = 1.0e-12_wp;  sle%warm_start = .true.
      st%dt_max = 200.0_wp*yr                      ! quarter of the 800 yr interval
      st%dt_min = 25.0_wp*yr                       ! the migrating coast makes the step-doubling
                                                   ! estimate noisy; accept at the floor
   end subroutine setup_ve

   subroutine check_elastic_step()
      type(response)        :: resp
      type(sle_solver)      :: s, s2
      type(adaptive_stepper) :: st
      type(sle_result)      :: res
      real(wp), allocatable :: rsl(:,:), C(:,:), rsl2(:,:), C2(:,:), dig(:,:), dig2(:,:)
      allocate(rsl(np,nl), C(np,nl), rsl2(np,nl), C2(np,nl), dig(np,nl), dig2(np,nl))
      call response_init_elastic(resp, e_el, lmax=LMAX)
      rsl = 0.0_wp;  rsl2 = 0.0_wp
      call stepper_advance(st, sht, resp, s, topo0, zero, ice1, zero, 0.0_wp, 100.0_wp*yr, rsl, C, &
                           host_ref=hr, host0=h0, host1=h1, dIg_out=dig)
      call sle_solve(s2, sht, resp, ice1, ice1, topo0, rsl2, C2, res, host=h1, host_ref=hr, dIg_out=dig2)
      write(*,'(a)') ' (1) elastic, one step: stepper == sle_solve at host1, ice1'
      write(*,'(a,es10.2,a,es10.2,a,es10.2)') '     |drsl| =', maxval(abs(rsl - rsl2)), ' m   |dC| =', &
           maxval(abs(C - C2)), '   |dIg| =', maxval(abs(dig - dig2))
      if (maxval(abs(rsl - rsl2)) > 1.0e-9_wp .or. maxval(abs(C - C2)) > 1.0e-12_wp &
          .or. maxval(abs(dig - dig2)) > 1.0e-9_wp) then
         write(*,'(a)') '     FAIL: stepper step differs from the direct solve';  ok = .false.
      end if
      if (maxval(abs(dig)) < 1.0_wp) then
         write(*,'(a)') '     FAIL: dIg_out empty';  ok = .false.
      end if
      call response_destroy(resp)
   end subroutine check_elastic_step

   subroutine check_dt0()
      type(response)        :: resp
      type(sle_solver)      :: s, s2, s3
      type(adaptive_stepper) :: st
      type(sle_result)      :: res
      real(wp), allocatable :: rsl(:,:), C(:,:), rsl1(:,:), rsl0(:,:), C1(:,:)
      real(wp), allocatable :: icep(:,:)
      allocate(rsl(np,nl), C(np,nl), rsl1(np,nl), rsl0(np,nl), C1(np,nl), icep(np,nl))
      icep = 0.4_wp*ice1                           ! the ice the model carries (ice0)
      call response_init_elastic(resp, e_el, lmax=LMAX)
      rsl = 0.0_wp;  rsl1 = 0.0_wp;  rsl0 = 0.0_wp
      call stepper_advance(st, sht, resp, s, topo0, icep, ice1, zero, 0.0_wp, 0.0_wp, rsl, C, &
                           host_ref=hr, host0=h0, host1=h1)
      call sle_solve(s2, sht, resp, icep, icep, topo0, rsl1, C1, res, report_only=.true., host=h1, host_ref=hr)
      call sle_solve(s3, sht, resp, icep, icep, topo0, rsl0, C1, res, report_only=.true., host=h0, host_ref=hr)
      write(*,'(a)') ' (2) dt = 0: seed with host1 and ice0'
      write(*,'(a,es10.2,a,es10.2)') '     |drsl - solve(host1)| =', maxval(abs(rsl - rsl1)), &
           ' m   |drsl - solve(host0)| =', maxval(abs(rsl - rsl0))
      if (maxval(abs(rsl - rsl1)) > 1.0e-9_wp) then
         write(*,'(a)') '     FAIL: dt = 0 is not the host1 solve';  ok = .false.
      end if
      if (maxval(abs(rsl - rsl0)) < 1.0e-3_wp) then
         write(*,'(a)') '     FAIL: host1 and host0 solves do not differ, the check is blind';  ok = .false.
      end if
      call response_destroy(resp)
   end subroutine check_dt0

   subroutine check_ramp()
      type(response)        :: resp_a, resp_b, resp_c
      type(sle_solver)      :: s_a, s_b, s_c
      type(adaptive_stepper) :: st_a, st_b, st_c
      type(sle_host)        :: ha, hb
      real(wp), allocatable :: rsl_a(:,:), C_a(:,:), dig_a(:,:), rsl_b(:,:), C_b(:,:), dig_b(:,:)
      real(wp), allocatable :: rsl_c(:,:), C_c(:,:), ia(:,:), ib(:,:)
      real(wp) :: dphi_a, dphi_b, t, span, da, dc, dm
      integer  :: k
      allocate(rsl_a(np,nl), C_a(np,nl), dig_a(np,nl), rsl_b(np,nl), C_b(np,nl), dig_b(np,nl), &
               rsl_c(np,nl), C_c(np,nl), ia(np,nl), ib(np,nl))
      span = 800.0_wp*yr
      ! a: one advance, four sub-steps
      call setup_ve(resp_a, s_a, st_a, .false.)
      rsl_a = 0.0_wp
      call stepper_advance(st_a, sht, resp_a, s_a, topo0, zero, ice1, zero, 0.0_wp, span, rsl_a, C_a, &
                           dphi_out=dphi_a, host_ref=hr, host0=h0, host1=h1, dIg_out=dig_a)
      ! b: four advances over the quarter intervals, endpoints interpolated by hand
      call setup_ve(resp_b, s_b, st_b, .false.)
      rsl_b = 0.0_wp
      do k = 1, 4
         call lerp(h0, h1, real(k-1,wp)/4.0_wp, ha);  call lerp(h0, h1, real(k,wp)/4.0_wp, hb)
         ia = real(k-1,wp)/4.0_wp*ice1;  ib = real(k,wp)/4.0_wp*ice1
         t = real(k-1,wp)*0.25_wp*span
         call stepper_advance(st_b, sht, resp_b, s_b, topo0, ia, ib, zero, t, t + 0.25_wp*span, rsl_b, C_b, &
                              dphi_out=dphi_b, host_ref=hr, host0=ha, host1=hb, dIg_out=dig_b)
      end do
      ! c: the host terms held at host0 instead of ramping
      call setup_ve(resp_c, s_c, st_c, .false.)
      rsl_c = 0.0_wp
      call stepper_advance(st_c, sht, resp_c, s_c, topo0, zero, ice1, zero, 0.0_wp, span, rsl_c, C_c, &
                           host_ref=hr, host0=h0, host1=h0)
      da = maxval(abs(rsl_a - rsl_b));  dc = maxval(abs(rsl_a - rsl_c))
      dm = abs(response_memory_norm(resp_a) - response_memory_norm(resp_b))
      write(*,'(a)') ' (3) ramp: one advance of four sub-steps vs four advances at interpolated ends'
      write(*,'(a,i0,a,i0,a)') '     sub-steps ', st_a%n_accept, ' and ', st_b%n_accept, ' (expect 4 and 4)'
      write(*,'(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') '     |drsl| =', da, ' m   |ddphi| =', abs(dphi_a - dphi_b), &
           '   |dIg| =', maxval(abs(dig_a - dig_b)), '   |dmem| =', dm
      write(*,'(a,es10.2,a)') '     host terms held at host0 instead: |drsl| =', dc, ' m'
      if (st_a%n_accept /= 4 .or. st_b%n_accept /= 4) then
         write(*,'(a)') '     FAIL: sub-step count';  ok = .false.
      end if
      if (da > 1.0e-9_wp .or. abs(dphi_a - dphi_b) > 1.0e-9_wp .or. maxval(abs(dig_a - dig_b)) > 1.0e-9_wp &
          .or. maxval(abs(C_a - C_b)) > 1.0e-12_wp) then
         write(*,'(a)') '     FAIL: ramp is not the interpolation of the host terms';  ok = .false.
      end if
      if (dc < 1.0e-3_wp) then
         write(*,'(a)') '     FAIL: ramped and held host terms agree, the ramp is not tested';  ok = .false.
      end if
      call response_destroy(resp_a);  call response_destroy(resp_b);  call response_destroy(resp_c)
   end subroutine check_ramp

   subroutine check_priming()
      type(response)        :: resp_a, resp_b, resp_c
      type(sle_solver)      :: s_a, s_b, s_c
      type(adaptive_stepper) :: st_a, st_b, st_c
      type(sle_result)      :: res
      real(wp), allocatable :: rsl_a(:,:), C_a(:,:), rsl_b(:,:), C_b(:,:), rsl_c(:,:), C_c(:,:), rp(:,:), Cp(:,:)
      complex(wp), allocatable :: sig0(:)
      real(wp) :: span, d0, d1
      allocate(rsl_a(np,nl), C_a(np,nl), rsl_b(np,nl), C_b(np,nl), rsl_c(np,nl), C_c(np,nl), &
               rp(np,nl), Cp(np,nl), sig0(sht%nlm))
      span = 400.0_wp*yr
      ! a: the stepper primes σ_n itself
      call setup_ve(resp_a, s_a, st_a, .true.)
      rsl_a = 0.0_wp
      call stepper_advance(st_a, sht, resp_a, s_a, topo0, 0.3_wp*ice1, ice1, zero, 0.0_wp, span, rsl_a, C_a, &
                           host_ref=hr, host0=h0, host1=h1)
      ! b: primed by hand with host0 and ice0
      call setup_ve(resp_b, s_b, st_b, .true.)
      rsl_b = 0.0_wp
      call sle_solve(s_b, sht, resp_b, 0.3_wp*ice1, 0.3_wp*ice1, topo0, rsl_b, C_b, res, report_only=.true., &
                     sigma_lm=sig0, host=h0, host_ref=hr)
      call response_prime_sigma(resp_b, sig0)
      call stepper_advance(st_b, sht, resp_b, s_b, topo0, 0.3_wp*ice1, ice1, zero, 0.0_wp, span, rsl_b, C_b, &
                           host_ref=hr, host0=h0, host1=h1)
      ! c: primed by hand with host1
      call setup_ve(resp_c, s_c, st_c, .true.)
      rsl_c = 0.0_wp
      call sle_solve(s_c, sht, resp_c, 0.3_wp*ice1, 0.3_wp*ice1, topo0, rsl_c, C_c, res, report_only=.true., &
                     sigma_lm=sig0, host=h1, host_ref=hr)
      call response_prime_sigma(resp_c, sig0)
      call stepper_advance(st_c, sht, resp_c, s_c, topo0, 0.3_wp*ice1, ice1, zero, 0.0_wp, span, rsl_c, C_c, &
                           host_ref=hr, host0=h0, host1=h1)
      d0 = maxval(abs(rsl_a - rsl_b));  d1 = maxval(abs(rsl_a - rsl_c))
      write(*,'(a)') ' (4) trapezoidal: σ_n priming at host0'
      write(*,'(a,es10.2,a,es10.2,a)') '     |drsl - primed by hand at host0| =', d0, &
           ' m   |drsl - primed at host1| =', d1, ' m'
      if (d0 > 1.0e-9_wp) then
         write(*,'(a)') '     FAIL: the stepper does not prime with host0';  ok = .false.
      end if
      if (d1 < 1.0e-6_wp) then
         write(*,'(a)') '     FAIL: priming at host0 and host1 agree, the check is blind';  ok = .false.
      end if
      call response_destroy(resp_a);  call response_destroy(resp_b);  call response_destroy(resp_c)
   end subroutine check_priming

end program test_host_ramp
