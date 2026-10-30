program test_bsl
   !! The two sea-level diagnostics of the host API, on a deglaciation whose coastline
   !! migrates (grounded marine ice on a shelf floods when it melts):
   !!
   !!   se%bsl     = ∫C·rsl dA / ∫C dA   (the ocean mean of C·rsl: the barystatic rise)
   !!   se%bsl_ive = −(ρ_i/ρ_w) ∫ΔI_g dA / ∫C dA   (ice-volume equivalent)
   !!
   !! Checks (after every real update, Gauss-native host):
   !!   (1) bsl equals the ocean mean of C·rsl formed here from se%gg%C, se%gg%rsl;
   !!   (2) the SLE's mass balance ties the two: ∫C·rsl dA = −(ρ_i/ρ_w)∫ΔI_g dA + ∫F dA
   !!       with the flooded-area fill F = z_bed_eq·(C − C0), so
   !!       bsl − bsl_ive = ∫F dA / ∫C dA, to the SLE's mass residual;
   !!   (4) the exported fill integral se%fill_int [m³] is ∫F dA, and with it the balance
   !!       closes to round-off: ∫C·rsl dA − fill_int = bsl_ive·∫C dA (the melt term); without
   !!       host terms the own part's share is all of it;
   !!   (3) the two differ while the shelf is grounded (C /= C0, so the fill is not
   !!       zero), so the rename is not a no-op: the water that fills newly submerged
   !!       area, or is displaced from it, does not move the barystatic level.
   use vilma_precision, only: wp
   use vilma_constants, only: pi
   use vilma_params,    only: vilma_param_class
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_surface_integral
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   implicit none

   integer, parameter :: LMAX = 16, NSTEP = 6
   real(wp), parameter :: TOL = 1.0e-8_wp               ! [m]
   type(sht_grid), target  :: sht
   type(vilma_param_class) :: p
   type(solid_earth)       :: se
   real(wp), allocatable   :: z_eq(:,:), h_eq(:,:), h_on(:,:), h_off(:,:)
   real(wp) :: c_int, bsl_ref, fill, err1, err2, worst1, worst2, max_gap
   real(wp) :: scl, a2, lhs, melt, err3, err4, worst3, worst4
   integer  :: i, j, step
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX+2, nphi=4*LMAX)
   allocate(z_eq(sht%nphi,sht%nlat), h_eq(sht%nphi,sht%nlat), &
            h_on(sht%nphi,sht%nlat), h_off(sht%nphi,sht%nlat))

   ! Reference: polar land cap (colat < 50, +500 m), shelf ring (50 <= colat < 62,
   ! -100 m), deep ocean (-4000 m); no reference ice. Load: 2 km on the cap
   ! (colat < 30) and 300 m of grounded marine ice on the shelf (931*300 > 1000*100).
   do j = 1, sht%nlat
      do i = 1, sht%nphi
         if (sht%colat(j)*180.0_wp/pi < 50.0_wp) then
            z_eq(i,j) = 500.0_wp
         else if (sht%colat(j)*180.0_wp/pi < 62.0_wp) then
            z_eq(i,j) = -100.0_wp
         else
            z_eq(i,j) = -4000.0_wp
         end if
         h_on(i,j) = 0.0_wp
         if (sht%colat(j)*180.0_wp/pi < 30.0_wp) h_on(i,j) = 2000.0_wp
         if (sht%colat(j)*180.0_wp/pi >= 50.0_wp .and. sht%colat(j)*180.0_wp/pi < 62.0_wp) h_on(i,j) = 300.0_wp
      end do
   end do
   h_eq = 0.0_wp;  h_off = 0.0_wp

   p%lmax = LMAX;  p%nlat = 2*LMAX+2;  p%nphi = 4*LMAX
   se%par = p
   call solid_earth_init(se, z_eq, h_eq)

   write(*,'(a,es8.1,a)') ' bsl: ocean mean of C*rsl, and its tie to bsl_ive; tolerance ', TOL, ' m'
   write(*,'(a)') '       step      bsl [m]   bsl_ive [m]   (bsl-ive) [m]   |(1)| [m]   |(2)| [m]'
   worst1 = 0.0_wp;  worst2 = 0.0_wp;  max_gap = 0.0_wp;  worst3 = 0.0_wp;  worst4 = 0.0_wp
   a2 = se%earth%r_earth**2
   do step = 1, 2*NSTEP                      ! build the load up, then take it away
      if (step <= NSTEP) then
         call solid_earth_update(se, h_on, 1.0e3_wp)
      else
         call solid_earth_update(se, h_off, 1.0e3_wp)
      end if
      c_int   = sht_grid_surface_integral(sht, se%gg%C)
      bsl_ref = sht_grid_surface_integral(sht, se%gg%C*se%gg%rsl)/c_int
      fill    = sht_grid_surface_integral(sht, z_eq*(se%gg%C - se%gg%C0))/c_int
      err1 = abs(se%bsl - bsl_ref)
      err2 = abs((se%bsl - se%bsl_ive) - fill)
      worst1 = max(worst1, err1);  worst2 = max(worst2, err2)
      ! (4) the exported fill integral [m³] and the balance with it
      fill = a2*sht_grid_surface_integral(sht, z_eq*(se%gg%C - se%gg%C0))
      lhs  = a2*sht_grid_surface_integral(sht, se%gg%C*se%gg%rsl)
      melt = a2*c_int*se%bsl_ive
      scl  = a2*sht_grid_surface_integral(sht, abs(se%gg%C*se%gg%rsl)) + abs(fill) + abs(melt) + 1.0e-30_wp
      err3 = abs(se%fill_int - fill)/scl
      err4 = abs(lhs - se%fill_int - melt)/scl
      if (se%fill_int_own /= se%fill_int) err3 = 1.0_wp
      worst3 = max(worst3, err3);  worst4 = max(worst4, err4)
      max_gap = max(max_gap, abs(se%bsl_ive - se%bsl))
      write(*,'(i9,4f14.6,2es12.2)') step, se%bsl, se%bsl_ive, se%bsl - se%bsl_ive, fill, err1, err2
   end do

   if (worst1 > TOL) then
      write(*,'(a)') '   FAIL: bsl is not the ocean mean of C*rsl'; ok = .false.
   end if
   if (worst2 > TOL) then
      write(*,'(a)') '   FAIL: bsl - bsl_ive does not match the flooded-area fill'; ok = .false.
   end if
   write(*,'(a,es10.2,a,es10.2)') ' fill_int vs recomputed (rel to scale):', worst3, '   balance with fill_int:', worst4
   if (worst3 > 1.0e-10_wp) then
      write(*,'(a)') '   FAIL: fill_int is not the fill integral'; ok = .false.
   end if
   if (worst4 > 1.0e-10_wp) then
      write(*,'(a)') '   FAIL: the ocean mass balance does not close with fill_int'; ok = .false.
   end if
   if (max_gap < 1.0e-3_wp) then
      write(*,'(a)') '   FAIL: bsl and bsl_ive never differ (no coastline change to tell them apart)'; ok = .false.
   end if

   call solid_earth_finalize(se)
   call sht_grid_destroy(sht)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: bsl is the ocean mean of C*rsl; bsl_ive is the ice-volume equivalent'
   else
      write(*,'(a)') ' FAIL: bsl / bsl_ive'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

end program test_bsl
