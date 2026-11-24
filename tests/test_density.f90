program test_density
   !! The densities are &vilma parameters: rho_ice, rho_water (melt water, the volume a melt
   !! adds to the ocean) and rho_seawater (flotation, the ocean load; vilma_sle), by default
   !! 917, 1000, 1028 (the GIA benchmarks' 931 and one water are set by the benchmark tests).
   !! Checks:
   !!
   !!   (1) vilma_par_load: the defaults file and the in-code defaults give 917 / 1000 / 1028,
   !!       a run config overrides them;
   !!   (2) the densities reach the whole model: the sea-level equation's volume balance
   !!       (bsl - bsl_ive = fill / ∫C, which uses the SLE's own ρ_i/ρ_w) holds against
   !!       bsl_ive formed here with ρ_i/ρ_w, and the ocean-volume residual of every solve
   !!       stays at round-off;
   !!   (3) flotation is the sea water's: the reference ocean function of a 100 m ice column
   !!       on 90 m of water is afloat at the defaults (917 x 100 < 1028 x 90 kg m^-2) and
   !!       a host's (910, 1000, 1028), grounded with one water of 1000 (91700 > 90000).
   use vilma_precision, only: wp
   use vilma_constants, only: pi
   use vilma_params,    only: vilma_param_class, vilma_par_load
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_surface_integral
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   implicit none

   integer, parameter :: LMAX = 16, NSTEP = 4
   real(wp), parameter :: TOL = 1.0e-8_wp               ! [m]
   real(wp), parameter :: RHO_ICE_DEF = 917.0_wp, RHO_WATER_DEF = 1000.0_wp, RHO_SEAWATER_DEF = 1028.0_wp   ! [kg m^-3]
   type(sht_grid), target :: sht
   real(wp), allocatable  :: z_eq(:,:), h_eq(:,:), h_on(:,:)
   integer :: ish, jsh
   logical :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX+2, nphi=4*LMAX)
   allocate(z_eq(sht%nphi,sht%nlat), h_eq(sht%nphi,sht%nlat), h_on(sht%nphi,sht%nlat))
   call make_fields()
   ish = 1;  jsh = nearest_row(55.0_wp)          ! a shelf cell

   call check_load()
   call run_case("default   ", RHO_ICE_DEF, RHO_WATER_DEF, RHO_SEAWATER_DEF, .true., 1.0_wp)
   call run_case("host      ", 910.0_wp, 1000.0_wp, 1028.0_wp, .false., 1.0_wp)
   call run_case("one water ", RHO_ICE_DEF, 1000.0_wp, 1000.0_wp, .false., 0.0_wp)

   call sht_grid_destroy(sht)
   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: densities are parameters; defaults are 917, 1000, 1028'
   else
      write(*,'(a)') ' FAIL: density parameters'
      error stop 1
   end if

contains

   subroutine check_load()
      !! (1) the in-code defaults, the nml route: defaults file, then a run config overriding both.
      type(vilma_param_class) :: p, p0
      character(len=*), parameter :: fn = "test_density.nml"
      integer :: u
      open(newunit=u, file=fn, status="replace", action="write")
      write(u,'(a)') "&vilma"
      write(u,'(a)') "  lmax = 16"
      write(u,'(a)') "  rho_ice = 910.0"
      write(u,'(a)') "  rho_water = 999.0"
      write(u,'(a)') "  rho_seawater = 1025.0"
      write(u,'(a)') "/"
      close(u)
      write(*,'(a,3f9.1)') ' density: in-code defaults rho_ice, rho_water, rho_seawater   =', p0%rho_ice, p0%rho_water, p0%rho_seawater
      if (p0%rho_ice /= RHO_ICE_DEF .or. p0%rho_water /= RHO_WATER_DEF .or. p0%rho_seawater /= RHO_SEAWATER_DEF) then
         write(*,'(a)') '   FAIL: in-code defaults are not 917 / 1000 / 1028'; ok = .false.
      end if
      call vilma_par_load(p, "input/vilma_defaults.nml")
      write(*,'(a,3f9.1)') ' density: defaults file gives rho_ice, rho_water, rho_seawater =', p%rho_ice, p%rho_water, p%rho_seawater
      if (p%rho_ice /= RHO_ICE_DEF .or. p%rho_water /= RHO_WATER_DEF .or. p%rho_seawater /= RHO_SEAWATER_DEF) then
         write(*,'(a)') '   FAIL: defaults file is not 917 / 1000 / 1028'; ok = .false.
      end if
      call vilma_par_load(p, fn, defaults_file="input/vilma_defaults.nml")
      write(*,'(a,3f9.1)') ' density: run config gives rho_ice, rho_water, rho_seawater =', p%rho_ice, p%rho_water, p%rho_seawater
      if (p%rho_ice /= 910.0_wp .or. p%rho_water /= 999.0_wp .or. p%rho_seawater /= 1025.0_wp) then
         write(*,'(a)') '   FAIL: run config did not override the densities'; ok = .false.
      end if
      call execute_command_line("rm -f "//fn)
   end subroutine check_load

   subroutine run_case(label, rho_i, rho_w, rho_sw, is_default, c0_expect)
      character(len=*), intent(in) :: label
      real(wp),         intent(in) :: rho_i, rho_w, rho_sw
      logical,          intent(in) :: is_default      !! leave the densities unset
      real(wp),         intent(in) :: c0_expect       !! the shelf's reference ocean function (1 afloat, 0 grounded)
      type(vilma_param_class) :: p
      type(solid_earth)       :: se
      real(wp) :: c_int, ive_ref, fill, worst_tie, worst_ive, worst_mass, c0_shelf
      integer  :: step

      p%lmax = LMAX;  p%nlat = 2*LMAX+2;  p%nphi = 4*LMAX
      if (.not. is_default) then;  p%rho_ice = rho_i;  p%rho_water = rho_w;  p%rho_seawater = rho_sw;  end if
      se%par = p
      call solid_earth_init(se, z_eq, h_eq)
      if (se%sle%rho_ice /= rho_i .or. se%sle%rho_water /= rho_w .or. se%sle%rho_seawater /= rho_sw) then
         write(*,'(a)') '   FAIL: densities did not reach the sea-level solver'; ok = .false.
      end if

      c0_shelf = se%gg%C0(ish,jsh)
      worst_tie = 0.0_wp;  worst_ive = 0.0_wp;  worst_mass = 0.0_wp
      do step = 1, NSTEP
         call solid_earth_update(se, h_on, 1.0e3_wp)
         c_int   = sht_grid_surface_integral(sht, se%gg%C)
         ive_ref = -(rho_i/rho_w) * sht_grid_surface_integral(sht, &
                      se%gg%h_ice*(1.0_wp - se%gg%C) - se%gg%h_ice_eq*(1.0_wp - se%gg%C0))/c_int
         fill    = sht_grid_surface_integral(sht, z_eq*(se%gg%C - se%gg%C0))/c_int
         worst_ive  = max(worst_ive, abs(se%bsl_ive - ive_ref))
         worst_tie  = max(worst_tie, abs((se%bsl - se%bsl_ive) - fill))
         worst_mass = max(worst_mass, se%worst_mass_resid)
      end do
      write(*,'(2a,f7.1,a,f7.1,a,f7.1,a,es9.2,a,es9.2,a,es9.2,a,f4.1)') ' density: ', label, rho_i, ' /', rho_w, ' /', rho_sw, &
           '  |bsl_ive err|', worst_ive, '  |mass tie|', worst_tie, '  mass resid', worst_mass, &
           '  C0(shelf)=', c0_shelf
      if (worst_ive > TOL .or. worst_tie > TOL) then
         write(*,'(a)') '   FAIL: sea-level budget does not use the parameter densities'; ok = .false.
      end if
      if (worst_mass > 1.0e-10_wp) then
         write(*,'(a)') '   FAIL: ocean volume not conserved'; ok = .false.
      end if
      ! (3) flotation by the sea water: afloat (C = 1) at 1028, grounded (C = 0) with one water of 1000
      if (c0_shelf /= c0_expect) then
         write(*,'(a,f4.1)') '   FAIL: the shelf ice''s reference ocean function should be ', c0_expect; ok = .false.
      end if
      call solid_earth_finalize(se)
   end subroutine run_case

   subroutine make_fields()
      !! Reference: polar land cap (colat < 50, +500 m), a shelf at -90 m (50 <= colat
      !! < 62), deep ocean elsewhere. Reference ice: 100 m on the shelf (grounded only with
      !! one water of 1000). Load: that, plus 2 km on the cap (colat < 30).
      integer :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            if (thd < 50.0_wp) then
               z_eq(i,j) = 500.0_wp
            else if (thd < 62.0_wp) then
               z_eq(i,j) = -90.0_wp
            else
               z_eq(i,j) = -4000.0_wp
            end if
            h_on(i,j) = 0.0_wp
            if (thd < 30.0_wp) h_on(i,j) = 2000.0_wp
            h_eq(i,j) = 0.0_wp
            if (thd >= 50.0_wp .and. thd < 62.0_wp) then
               h_eq(i,j) = 100.0_wp
               h_on(i,j) = 100.0_wp
            end if
         end do
      end do
   end subroutine make_fields

   integer function nearest_row(colat_deg) result(jbest)
      real(wp), intent(in) :: colat_deg
      integer  :: j
      real(wp) :: d, dbest
      jbest = 1;  dbest = huge(1.0_wp)
      do j = 1, sht%nlat
         d = abs(sht%colat(j)*180.0_wp/pi - colat_deg)
         if (d < dbest) then;  dbest = d;  jbest = j;  end if
      end do
   end function nearest_row

end program test_density
