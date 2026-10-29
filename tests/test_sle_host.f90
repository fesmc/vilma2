program test_sle_host
   !! Host-resolved terms in the sea-level solve (sle_solve with host, host_ref): a host
   !! region of weight w, where the host hands over its own grounded-ice column I_g,
   !! ocean fraction C and flooded-area fill F, and the solve applies its own flotation
   !! and coastline rule on the remaining 1 - w. Self-gravitating elastic response.
   !!
   !! Checks:
   !!   (1) w = 0, zero terms: the host formulas reduce to the host-free solve (rsl, C,
   !!       dphi, dIg_out to round-off; the coupling dispatches around this case, so
   !!       bit-identity is not asked here);
   !!   (2) w = 1 everywhere with the host terms of a fixed coastline == the host-free
   !!       fixed-ocean solve; own ice is NaN there and must not leak; the outer loop
   !!       ends after one pass;
   !!   (3) w = 1/2 with the host terms of half the own state, no coastline motion ==
   !!       the host-free solve of the whole state;
   !!   (4) w varying (0, 1/2, 1 in sectors), the host's coastline offset from the own
   !!       one, a grounded marine margin that floods and unfloods the own part:
   !!       the ocean mass balance  ∫C·rsl = -(rho_i/rho_w) ∫ΔI_g + ∫fill  to round-off,
   !!       formed here from the returned fields, the own coastline converged, and the
   !!       fields finite.
   use vilma_precision,       only: wp
   use vilma_constants,       only: rho_ice, rho_water, pi
   use vilma_earth_structure, only: earth_model, build_M3L70V01
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response_destroy, response, response_init_elastic
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_surface_integral
   use vilma_sle,             only: sle_solve, sle_solver, sle_result, sle_host, ocean_function
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none

   integer,  parameter :: LMAX = 16
   real(wp), parameter :: DEG = pi/180.0_wp
   type(sht_grid)      :: sht
   type(earth_model)   :: e
   type(response)      :: resp
   logical             :: ok
   integer             :: np, nl

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX+2, nphi=4*LMAX)
   np = sht%nphi;  nl = sht%nlat
   e = build_M3L70V01()
   call response_init_elastic(resp, e, lmax=LMAX)

   call check_zero_weight(ok)
   call check_full_weight(ok)
   call check_half_weight(ok)
   call check_mass_balance(ok)

   write(*,'(a)') ''
   call response_destroy(resp)
   call sht_grid_destroy(sht)
   if (ok) then
      write(*,'(a)') ' PASS: host-resolved terms reduce to the host-free solve and conserve ocean mass'
      call radial_fe_finalize()
   else
      write(*,'(a)') ' FAIL: host-resolved terms in the sea-level solve'
      call radial_fe_finalize()
      error stop 1
   end if

contains

   subroutine solver(sle)
      !! Converge deeply, so that an iteration count that differs by one between two
      !! mathematically equal formulations cannot show up in the comparisons.
      type(sle_solver), intent(out) :: sle
      sle%n_outer = 12;  sle%n_inner = 80;  sle%tol = 1.0e-13_wp
   end subroutine solver

   subroutine make_fields(topo0, ice, shelf)
      !! Reference bed: land cap (colat < 50, +500 m), shelf ring (50..62, -100 m), bare
      !! shallow ring (62..70, -60 m), deep ocean (-4000 m); shelf = .false. makes both
      !! rings deep ocean, so that no coast can move. Ice: 2 km on the cap (colat < 30)
      !! and 300 m on the shelf ring (grounded: 931*300 > 1000*100). The sea level falls
      !! by about 100 m under this load, so the bare ring emerges in the coastline loop.
      real(wp), intent(out) :: topo0(:,:), ice(:,:)
      logical,  intent(in)  :: shelf
      integer  :: i, j
      real(wp) :: th
      do j = 1, nl
         th = sht%colat(j)/DEG
         do i = 1, np
            if (th < 50.0_wp) then
               topo0(i,j) = 500.0_wp
            else if (th < 62.0_wp .and. shelf) then
               topo0(i,j) = -100.0_wp
            else if (th < 70.0_wp .and. shelf) then
               topo0(i,j) = -60.0_wp
            else
               topo0(i,j) = -4000.0_wp
            end if
            ice(i,j) = 0.0_wp
            if (th < 30.0_wp) ice(i,j) = 2000.0_wp
            if (th >= 50.0_wp .and. th < 62.0_wp .and. shelf) ice(i,j) = 300.0_wp
         end do
      end do
   end subroutine make_fields

   subroutine alloc_host(h, with_f)
      type(sle_host), intent(inout) :: h
      logical,        intent(in)    :: with_f
      allocate(h%w(np,nl), h%I_g(np,nl), h%C(np,nl), source=0.0_wp)
      if (with_f) allocate(h%F(np,nl), source=0.0_wp)
   end subroutine alloc_host

   subroutine check_zero_weight(ok)
      logical, intent(inout) :: ok
      type(sle_solver) :: s0, s1
      type(sle_result) :: r0, r1
      type(sle_host)   :: h, hr
      real(wp), allocatable :: topo0(:,:), ice(:,:), rsl0(:,:), C0(:,:), rsl1(:,:), C1(:,:), dig(:,:)
      real(wp) :: dr, dd
      allocate(topo0(np,nl), ice(np,nl), rsl0(np,nl), C0(np,nl), rsl1(np,nl), C1(np,nl), dig(np,nl))
      call make_fields(topo0, ice, .true.)
      call solver(s0);  call solver(s1)
      call alloc_host(h, .true.);  call alloc_host(hr, .false.)
      call sle_solve(s0, sht, resp, ice, ice, topo0, rsl0, C0, r0)
      call sle_solve(s1, sht, resp, ice, ice, topo0, rsl1, C1, r1, host=h, host_ref=hr, dIg_out=dig)
      dr = maxval(abs(rsl1 - rsl0));  dd = abs(r1%esl - r0%esl)
      write(*,'(a)') ' (1) w = 0, zero terms vs host-free'
      write(*,'(a,es10.2,a,es10.2,a,es10.2,a,i0)') '     |drsl| =', dr, ' m   |ddphi| =', dd, &
           ' m   |dC| =', maxval(abs(C1 - C0)), '   outer ', r1%n_outer_done
      if (dr > 1.0e-9_wp .or. dd > 1.0e-9_wp .or. any(C1 /= C0)) then
         write(*,'(a)') '     FAIL: w = 0 differs from the host-free solve';  ok = .false.
      end if
      if (r1%n_outer_done /= r0%n_outer_done .or. r1%n_coast_flip /= 0) then
         write(*,'(a)') '     FAIL: own coastline loop differs from the host-free one';  ok = .false.
      end if
      if (r0%n_outer_done < 2) then
         write(*,'(a)') '     FAIL: scenario does not move the coastline';  ok = .false.
      end if
   end subroutine check_zero_weight

   subroutine check_full_weight(ok)
      logical, intent(inout) :: ok
      type(sle_solver) :: s0, s1
      type(sle_result) :: r0, r1
      type(sle_host)   :: h, hr
      real(wp), allocatable :: topo0(:,:), ice(:,:), Cref(:,:), rsl0(:,:), C0(:,:), rsl1(:,:), C1(:,:), dig(:,:)
      real(wp), allocatable :: junk(:,:)
      real(wp) :: dr, dd
      allocate(topo0(np,nl), ice(np,nl), Cref(np,nl), rsl0(np,nl), C0(np,nl), rsl1(np,nl), C1(np,nl), &
               dig(np,nl), junk(np,nl))
      call make_fields(topo0, ice, .true.)
      call ocean_function(topo0, ice*0.0_wp, Cref, rho_ice, rho_water)   ! reference coast: ice-free
      ! host-free reference: the coastline held at the reference one
      call solver(s0);  s0%fixed_ocean = .true.
      call sle_solve(s0, sht, resp, ice, ice, topo0, rsl0, C0, r0)
      ! the same state through the host: w = 1, the host's columns at the held coastline
      call solver(s1)
      call alloc_host(h, .true.);  call alloc_host(hr, .false.)
      h%w = 1.0_wp;  hr%w = 1.0_wp
      h%C = Cref;  hr%C = Cref
      h%I_g = ice*(1.0_wp - Cref);  hr%I_g = 0.0_wp
      junk = ieee_value(1.0_wp, ieee_quiet_nan)          ! own ice and topography where w = 1
      call sle_solve(s1, sht, resp, junk, junk, topo0, rsl1, C1, r1, host=h, host_ref=hr, dIg_out=dig)
      dr = maxval(abs(rsl1 - rsl0));  dd = abs(r1%esl - r0%esl)
      write(*,'(a)') ' (2) w = 1, host terms of a held coastline vs host-free fixed-ocean'
      write(*,'(a,es10.2,a,es10.2,a,es10.2,a,i0)') '     |drsl| =', dr, ' m   |ddphi| =', dd, &
           ' m   |dC| =', maxval(abs(C1 - C0)), '   outer ', r1%n_outer_done
      if (dr > 1.0e-9_wp .or. dd > 1.0e-9_wp .or. any(C1 /= C0)) then
         write(*,'(a)') '     FAIL: w = 1 differs from the host-free fixed-ocean solve';  ok = .false.
      end if
      if (r1%n_outer_done /= 1) then
         write(*,'(a)') '     FAIL: w = 1 should end after one coastline pass';  ok = .false.
      end if
      if (any(dig /= h%I_g - hr%I_g)) then
         write(*,'(a)') '     FAIL: dIg_out is not I_g - I_g0 where w = 1';  ok = .false.
      end if
      if (any(rsl1 /= rsl1)) then
         write(*,'(a)') '     FAIL: own-part NaN leaked into rsl';  ok = .false.
      end if
   end subroutine check_full_weight

   subroutine check_half_weight(ok)
      logical, intent(inout) :: ok
      type(sle_solver) :: s0, s1
      type(sle_result) :: r0, r1
      type(sle_host)   :: h, hr
      real(wp), allocatable :: topo0(:,:), ice(:,:), Cref(:,:), rsl0(:,:), C0(:,:), rsl1(:,:), C1(:,:)
      real(wp) :: dr, dd
      allocate(topo0(np,nl), ice(np,nl), Cref(np,nl), rsl0(np,nl), C0(np,nl), rsl1(np,nl), C1(np,nl))
      call make_fields(topo0, ice, .false.)               ! no shelf: the coast cannot move
      call ocean_function(topo0, ice*0.0_wp, Cref, rho_ice, rho_water)
      call solver(s0)
      call sle_solve(s0, sht, resp, ice, ice, topo0, rsl0, C0, r0)
      call solver(s1)
      call alloc_host(h, .true.);  call alloc_host(hr, .false.)
      h%w = 0.5_wp;  hr%w = 0.5_wp
      h%C = 0.5_wp*Cref;  hr%C = 0.5_wp*Cref
      h%I_g = 0.5_wp*ice*(1.0_wp - Cref);  hr%I_g = 0.0_wp
      call sle_solve(s1, sht, resp, ice, ice, topo0, rsl1, C1, r1, host=h, host_ref=hr)
      dr = maxval(abs(rsl1 - rsl0));  dd = abs(r1%esl - r0%esl)
      write(*,'(a)') ' (3) w = 1/2, host terms of half the own state vs host-free'
      write(*,'(a,es10.2,a,es10.2,a,es10.2)') '     |drsl| =', dr, ' m   |ddphi| =', dd, &
           ' m   |dC| =', maxval(abs(C1 - C0))
      if (r0%n_outer_done /= 1) then
         write(*,'(a)') '     FAIL: scenario moved the coastline';  ok = .false.
      end if
      if (dr > 1.0e-9_wp .or. dd > 1.0e-9_wp .or. maxval(abs(C1 - C0)) > 1.0e-14_wp) then
         write(*,'(a)') '     FAIL: split state differs from the whole one';  ok = .false.
      end if
   end subroutine check_half_weight

   subroutine check_mass_balance(ok)
      logical, intent(inout) :: ok
      type(sle_solver) :: s
      type(sle_result) :: r
      type(sle_host)   :: h, hr
      real(wp), allocatable :: topo0(:,:), ice(:,:), dice(:,:), rsl(:,:), C(:,:), dig(:,:)
      real(wp), allocatable :: Cown(:,:), C0own(:,:), fill(:,:), Cref(:,:), fg(:,:), fg0(:,:), a(:,:)
      real(wp) :: th, lon, lhs, rhs, melt, fillint, scale
      integer  :: i, j
      allocate(topo0(np,nl), ice(np,nl), dice(np,nl), rsl(np,nl), C(np,nl), dig(np,nl), &
               Cown(np,nl), C0own(np,nl), fill(np,nl), Cref(np,nl), fg(np,nl), fg0(np,nl), a(np,nl))
      call make_fields(topo0, ice, .true.)
      call alloc_host(h, .true.);  call alloc_host(hr, .false.)
      do j = 1, nl
         th = sht%colat(j)/DEG
         do i = 1, np
            lon = sht%lon(i)/DEG
            if (th >= 20.0_wp .and. th < 70.0_wp) then
               if (lon < 90.0_wp) then
                  h%w(i,j) = 1.0_wp
               else if (lon < 180.0_wp) then
                  h%w(i,j) = 0.5_wp
               end if
            end if
            ! the host sees land to colat 54 now (the own part's coast is at 50 and 62) and,
            ! at the reference, to colat 58
            fg(i,j)  = merge(1.0_wp, 0.0_wp, th < 54.0_wp)
            fg0(i,j) = merge(1.0_wp, 0.0_wp, th < 58.0_wp)
         end do
      end do
      hr%w = h%w
      ! host terms: I_g = w H f_g with its own ice (200 m beyond colat 30 on its land, 1.5x
      ! the own cap), C = w (1 - f_g), F = w zeta0 (C - C_ref) with its own bed (+200 m
      ! land, -50 m shelf)
      where (h%w > 0.0_wp)
         h%I_g  = h%w*fg*merge(3000.0_wp, 200.0_wp, sht_colat_grid() < 30.0_wp)
         h%C    = h%w*(1.0_wp - fg)
         hr%C   = hr%w*(1.0_wp - fg0)
         h%F    = h%w*merge(200.0_wp, -50.0_wp, fg > 0.5_wp)*((1.0_wp - fg) - (1.0_wp - fg0))
      end where
      ! own part: the cap and shelf ice above, reference ice-free (d_ice = ice)
      dice = ice
      call solver(s)
      call sle_solve(s, sht, resp, dice, ice, topo0, rsl, C, r, host=h, host_ref=hr, dIg_out=dig)

      ! independent balance from the returned fields
      a = 1.0_wp - h%w
      call ocean_function(topo0, ice - dice, C0own, rho_ice, rho_water)
      where (a > 0.0_wp)
         Cown = (C - h%C)/a
      elsewhere
         Cown = 0.0_wp
      end where
      where (a > 0.0_wp)
         fill = h%F + a*topo0*(Cown - C0own)
      elsewhere
         fill = h%F
      end where
      lhs     = sht_grid_surface_integral(sht, C*rsl)
      melt    = -(rho_ice/rho_water)*sht_grid_surface_integral(sht, dig)
      fillint = sht_grid_surface_integral(sht, fill)
      rhs     = melt + fillint
      scale   = max(abs(melt), abs(fillint), 1.0e-30_wp)
      write(*,'(a)') ' (4) varying w, host coast offset, own marine margin floods'
      write(*,'(a,es11.3,a,es11.3,a,es11.3,a,es10.2)') '     int C rsl =', lhs, '   melt =', melt, &
           '   fill =', fillint, '   lhs-rhs =', lhs - rhs
      write(*,'(a,es10.2,a,es10.2,a,i0,a,i0)') '     balance error / scale =', abs(lhs - rhs)/scale, &
           '   mass_resid =', r%mass_resid, '   outer ', r%n_outer_done, '   flips ', r%n_coast_flip
      if (abs(lhs - rhs)/scale > 1.0e-10_wp) then
         write(*,'(a)') '     FAIL: ocean mass balance';  ok = .false.
      end if
      if (r%mass_resid > 1.0e-10_wp) then
         write(*,'(a)') '     FAIL: solver mass residual';  ok = .false.
      end if
      if (r%n_coast_flip /= 0 .or. r%n_outer_done < 2) then
         write(*,'(a)') '     FAIL: own coastline did not move and settle';  ok = .false.
      end if
      if (abs(fillint) < 1.0e-3_wp*abs(melt)) then
         write(*,'(a)') '     FAIL: fill is negligible, the scenario does not test it';  ok = .false.
      end if
      if (any(rsl /= rsl) .or. any(C < 0.0_wp) .or. any(C > 1.0_wp + 1.0e-12_wp)) then
         write(*,'(a)') '     FAIL: rsl not finite or C outside [0,1]';  ok = .false.
      end if
   end subroutine check_mass_balance

   function sht_colat_grid() result(th)
      real(wp) :: th(np,nl)
      integer  :: i, j
      do j = 1, nl
         do i = 1, np
            th(i,j) = sht%colat(j)/DEG
         end do
      end do
   end function sht_colat_grid

end program test_sle_host
