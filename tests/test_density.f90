program test_density
   !! Ice and water densities are &vilma parameters (rho_ice, rho_water); the GIA
   !! benchmark values (931, 1000) are the defaults. Checks:
   !!
   !!   (1) vilma_par_load: the defaults file gives 931 / 1000, a run config overrides them;
   !!   (2) a host pair (910, 1028) reaches the whole model: the sea-level equation's mass
   !!       balance (bsl - bsl_ive = fill / ∫C, which uses the SLE's own ρ_i/ρ_w) holds
   !!       against bsl_ive formed here with that pair, and the ocean-mass residual of
   !!       every solve stays at round-off;
   !!   (3) flotation uses the pair: the reference ocean function of a 100 m ice column
   !!       on 92 m of water is grounded at (931, 1000) (93100 > 92000 kg m^-2) and afloat
   !!       at (910, 1028) (91000 < 94576).
   use vilma_precision, only: wp
   use vilma_constants, only: pi, rho_ice, rho_water
   use vilma_params,    only: vilma_param_class, vilma_par_load
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_surface_integral
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   implicit none

   integer, parameter :: LMAX = 16, NSTEP = 4
   real(wp), parameter :: TOL = 1.0e-8_wp               ! [m]
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
   call run_pair("benchmark ", rho_ice, rho_water, .true.)
   call run_pair("host pair ", 910.0_wp, 1028.0_wp, .false.)

   call sht_grid_destroy(sht)
   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: densities are parameters; defaults are the GIA benchmark values'
   else
      write(*,'(a)') ' FAIL: density parameters'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine check_load()
      !! (1) the nml route: defaults file, then a run config overriding both.
      type(vilma_param_class) :: p
      character(len=*), parameter :: fn = "test_density.nml"
      integer :: u
      open(newunit=u, file=fn, status="replace", action="write")
      write(u,'(a)') "&vilma"
      write(u,'(a)') "  lmax = 16"
      write(u,'(a)') "  rho_ice = 910.0"
      write(u,'(a)') "  rho_water = 1028.0"
      write(u,'(a)') "/"
      close(u)
      call vilma_par_load(p, "input/vilma_defaults.nml")
      write(*,'(a,2f9.1)') ' density: defaults file gives rho_ice, rho_water =', p%rho_ice, p%rho_water
      if (p%rho_ice /= rho_ice .or. p%rho_water /= rho_water) then
         write(*,'(a)') '   FAIL: defaults are not the benchmark values'; ok = .false.
      end if
      call vilma_par_load(p, fn, defaults_file="input/vilma_defaults.nml")
      write(*,'(a,2f9.1)') ' density: run config gives rho_ice, rho_water =', p%rho_ice, p%rho_water
      if (p%rho_ice /= 910.0_wp .or. p%rho_water /= 1028.0_wp) then
         write(*,'(a)') '   FAIL: run config did not override the densities'; ok = .false.
      end if
      call execute_command_line("rm -f "//fn)
   end subroutine check_load

   subroutine run_pair(label, rho_i, rho_w, is_default)
      character(len=*), intent(in) :: label
      real(wp),         intent(in) :: rho_i, rho_w
      logical,          intent(in) :: is_default      !! leave the densities unset
      type(vilma_param_class) :: p
      type(solid_earth)       :: se
      real(wp) :: c_int, ive_ref, fill, worst_tie, worst_ive, worst_mass, c0_shelf
      integer  :: step

      p%lmax = LMAX;  p%nlat = 2*LMAX+2;  p%nphi = 4*LMAX
      if (.not. is_default) then;  p%rho_ice = rho_i;  p%rho_water = rho_w;  end if
      se%par = p
      call solid_earth_init(se, z_eq, h_eq)
      if (se%sle%rho_ice /= rho_i .or. se%sle%rho_water /= rho_w) then
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
      write(*,'(2a,f7.1,a,f7.1,a,es9.2,a,es9.2,a,es9.2,a,f4.1)') ' density: ', label, rho_i, ' /', rho_w, &
           '  |bsl_ive err|', worst_ive, '  |mass tie|', worst_tie, '  mass resid', worst_mass, &
           '  C0(shelf)=', c0_shelf
      if (worst_ive > TOL .or. worst_tie > TOL) then
         write(*,'(a)') '   FAIL: sea-level budget does not use the parameter densities'; ok = .false.
      end if
      if (worst_mass > 1.0e-10_wp) then
         write(*,'(a)') '   FAIL: ocean mass not conserved'; ok = .false.
      end if
      ! (3) flotation: grounded (C = 0) at the benchmark pair, afloat (C = 1) at the host pair
      if (is_default .and. c0_shelf /= 0.0_wp) then
         write(*,'(a)') '   FAIL: shelf ice should ground at 931/1000'; ok = .false.
      end if
      if (.not. is_default .and. c0_shelf /= 1.0_wp) then
         write(*,'(a)') '   FAIL: shelf ice should float at 910/1028'; ok = .false.
      end if
      call solid_earth_finalize(se)
   end subroutine run_pair

   subroutine make_fields()
      !! Reference: polar land cap (colat < 50, +500 m), a shelf at -92 m (50 <= colat
      !! < 62), deep ocean elsewhere. Reference ice: 100 m on the shelf (grounded only at
      !! the benchmark pair). Load: that, plus 2 km on the cap (colat < 30).
      integer :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            if (thd < 50.0_wp) then
               z_eq(i,j) = 500.0_wp
            else if (thd < 62.0_wp) then
               z_eq(i,j) = -92.0_wp
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
