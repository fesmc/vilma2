program test_dz_split
   !! The split of rsl into its bed and sea-surface parts at the host API: after every
   !! solid_earth_update (and the dt = 0 seed)
   !!
   !!     rsl = dz_ss - dz_bed ,    dz_bed = u ,    dz_ss = N + Δφ ,
   !!
   !! to round-off, on the Gauss grid (se%gg) and on the host grid (se%rsl, se%dz_*),
   !! for the explicit (fe) and implicit (trap) stepper and with and without a host
   !! remap. The load is a grounded land cap plus grounded marine ice on a shelf,
   !! built up and then removed, so the coastline migrates and Δφ is non-trivial.
   !! Also checks the parts are not trivially zero (the identity would hold for 0 − 0)
   !! and that the reference state has dz_bed = dz_ss = 0.
   use vilma_precision, only: wp
   use vilma_constants, only: pi
   use vilma_params,    only: vilma_param_class
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   use coords,          only: grid_class, grid_init
   implicit none

   integer, parameter :: LMAX = 16, NSTEP = 6
   integer, parameter :: NLON = 96, NLAT = 48           ! host grid for the remap case
   real(wp), parameter :: TOL = 1.0e-13_wp              ! relative to max|rsl|
   type(sht_grid), target :: sht
   type(grid_class)       :: hgrid
   real(wp), allocatable  :: lon(:), lat(:)
   logical :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX+2, nphi=4*LMAX)

   write(*,'(a,es8.1,a)') ' dz split: rsl = dz_ss - dz_bed, tolerance ', TOL, ' x max|rsl|'

   ! Gauss-native host (passthrough)
   call run_case("fe,   Gauss-native ", "fe",   .false.)
   call run_case("trap, Gauss-native ", "trap", .false.)
   ! lon-lat host grid (conservative in, bilinear out)
   allocate(lon(NLON), lat(NLAT))
   call build_host_grid()
   call run_case("fe,   lon-lat host ", "fe",   .true.)

   call sht_grid_destroy(sht)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: rsl = dz_ss - dz_bed to round-off on the Gauss and host grids'
   else
      write(*,'(a)') ' FAIL: dz_bed / dz_ss split'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine run_case(label, scheme, use_host)
      character(len=*), intent(in) :: label, scheme
      logical,          intent(in) :: use_host
      type(vilma_param_class) :: p
      type(solid_earth)       :: se
      real(wp), allocatable   :: z_eq(:,:), h_eq(:,:), h_on(:,:), h_off(:,:), lo(:), la(:)
      real(wp) :: err_g, err_h, scale, worst_g, worst_h, max_u, max_n, max_dphi
      integer  :: nx, ny, step

      p%lmax = LMAX;  p%nlat = 2*LMAX+2;  p%nphi = 4*LMAX
      p%scheme = scheme
      if (use_host) then
         nx = NLON;  ny = NLAT
         allocate(lo(nx), la(ny));  lo = lon;  la = lat
      else
         nx = sht%nphi;  ny = sht%nlat
         allocate(lo(nx), la(ny))
         lo = sht%lon*180.0_wp/pi
         la = 90.0_wp - sht%colat*180.0_wp/pi
      end if
      allocate(z_eq(nx,ny), h_eq(nx,ny), h_on(nx,ny), h_off(nx,ny))
      call make_fields(lo, la, z_eq, h_eq, h_on)
      h_off = h_eq

      se%par = p
      if (use_host) then
         call solid_earth_init(se, z_eq, h_eq, grid=hgrid)
      else
         call solid_earth_init(se, z_eq, h_eq)
      end if

      ! reference: nothing has moved
      if (maxval(abs(se%dz_bed)) /= 0.0_wp .or. maxval(abs(se%dz_ss)) /= 0.0_wp .or. &
          maxval(abs(se%gg%dz_bed)) /= 0.0_wp .or. maxval(abs(se%gg%dz_ss)) /= 0.0_wp) then
         write(*,'(a)') '   FAIL: dz_bed / dz_ss not zero at the reference'; ok = .false.
      end if

      worst_g = 0.0_wp;  worst_h = 0.0_wp;  max_u = 0.0_wp;  max_n = 0.0_wp;  max_dphi = 0.0_wp
      do step = 0, 2*NSTEP
         ! step 0: the dt = 0 report-only seed; then build the load up, then take it away
         if (step == 0) then
            call solid_earth_update(se, h_on, 0.0_wp)
         else if (step <= NSTEP) then
            call solid_earth_update(se, h_on, 1.0e3_wp)
         else
            call solid_earth_update(se, h_off, 1.0e3_wp)
         end if
         scale = max(maxval(abs(se%gg%rsl)), 1.0_wp)
         err_g = maxval(abs(se%gg%rsl - (se%gg%dz_ss - se%gg%dz_bed)))/scale
         scale = max(maxval(abs(se%rsl)), 1.0_wp)
         err_h = maxval(abs(se%rsl - (se%dz_ss - se%dz_bed)))/scale
         worst_g = max(worst_g, err_g);  worst_h = max(worst_h, err_h)
         max_u = max(max_u, maxval(abs(se%gg%dz_bed)))
         max_n = max(max_n, maxval(abs(se%gg%dz_ss)))
         max_dphi = max(max_dphi, abs(se%dphi))
      end do

      write(*,'(2a,a,es9.2,a,es9.2,a,f8.2,a,f8.2,a,f7.2)') '   ', label, &
           ' worst |rsl-(dz_ss-dz_bed)|/max|rsl|: Gauss ', worst_g, '  host ', worst_h, &
           '   max|u|=', max_u, ' max|N+dphi|=', max_n, ' max|dphi|=', max_dphi
      if (worst_g > TOL .or. worst_h > TOL) then
         write(*,'(a)') '   FAIL: rsl /= dz_ss - dz_bed beyond round-off'; ok = .false.
      end if
      if (max_u < 1.0e-3_wp .or. max_n < 1.0e-3_wp .or. max_dphi < 1.0e-3_wp) then
         write(*,'(a)') '   FAIL: dz_bed, dz_ss or dphi trivially zero'; ok = .false.
      end if
      call solid_earth_finalize(se)

   end subroutine run_case

   subroutine make_fields(lo, la, z_eq, h_eq, h_on)
      !! Reference: polar land cap (lat > 40, +500 m), a shelf ring (28 < lat <= 40,
      !! -100 m) and deep ocean (-4000 m); no reference ice. Load: 2 km on the cap
      !! (lat > 60), 300 m of grounded marine ice on the shelf (931*300 > 1000*100).
      real(wp), intent(in)  :: lo(:), la(:)
      real(wp), intent(out) :: z_eq(:,:), h_eq(:,:), h_on(:,:)
      integer :: i, j
      do j = 1, size(la)
         do i = 1, size(lo)
            if (la(j) > 40.0_wp) then
               z_eq(i,j) = 500.0_wp
            else if (la(j) > 28.0_wp) then
               z_eq(i,j) = -100.0_wp
            else
               z_eq(i,j) = -4000.0_wp
            end if
            h_on(i,j) = 0.0_wp
            if (la(j) > 60.0_wp) h_on(i,j) = 2000.0_wp
            if (la(j) > 28.0_wp .and. la(j) <= 40.0_wp) h_on(i,j) = 300.0_wp
         end do
      end do
      h_eq = 0.0_wp
   end subroutine make_fields

   subroutine build_host_grid()
      integer :: i, j
      do i = 1, NLON;  lon(i) = -180.0_wp + (real(i,wp)-0.5_wp)*(360.0_wp/NLON);  end do
      do j = 1, NLAT;  lat(j) =  -90.0_wp + (real(j,wp)-0.5_wp)*(180.0_wp/NLAT);  end do
      call grid_init(hgrid, name="host", mtype="latlon", units="degrees", planet="Spherical Earth", &
                     lon180=.true., x=lon, y=lat)
   end subroutine build_host_grid

end program test_dz_split
