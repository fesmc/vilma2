program test_restart_3d
   !! Restart exactness with laterally varying viscosity (3-D, toroidal channels, rotation on),
   !! lmax 16: a run restored mid-way continues the uninterrupted one bit for bit, at every
   !! step after the restore, in the sea-level fields, the scalars, the fill integrals and the
   !! memory. Two cases:
   !!
   !!   (a) the plain 3-D run;
   !!   (b) the 3-D run entered through the 1-D pre-spin-up (&vilma pre_spinup_1d), which
   !!       enables the lateral mean of the viscosity field and then the field itself. The
   !!       restored run is a fresh model that never spun up, so the viscosity the response
   !!       carries must not depend on the enables before the last one. (It did: the rates
   !!       were rebuilt from the previous enable's, a round-off off the radial ones, and a
   !!       restart continued at round-off, not exactly. Checked here on the rates, too.)
   !!   (c) as (b) with the implicit (trapezoidal) scheme, whose adaptive stepper sub-steps by the
   !!       restored dt_try (which the restart used to round-trip through years, inexactly).
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi
   use vilma_params,          only: vilma_param_class
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_coupling,        only: solid_earth_finalize, solid_earth_update, solid_earth_init, solid_earth, &
                                    solid_earth_spinup
   use vilma_io,              only: vilma_restart_write, vilma_restart_read
   implicit none

   integer, parameter :: LMAX = 16, K1 = 2, K2 = 3
   character(len=*), parameter :: VISC = "libvilma/include/test_restart_3d_visc.nc"
   type(sht_grid), target :: sht
   real(wp), allocatable  :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
   logical  :: ok

   ok = .true.
   call write_synth_visc_fesmdata(VISC)
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   allocate(z_bed_eq(sht%nphi,sht%nlat), h_ice_eq(sht%nphi,sht%nlat), h_ice(sht%nphi,sht%nlat))
   call make_fields(z_bed_eq, h_ice_eq, h_ice)

   call restart_case("plain 3-D",                 .false., "fe",   "libvilma/include/test_restart_3d_a.nc", ok)
   call restart_case("1-D pre-spin-up, then 3-D", .true.,  "fe",   "libvilma/include/test_restart_3d_b.nc", ok)
   call restart_case("trap, pre-spin-up, 3-D",    .true.,  "trap", "libvilma/include/test_restart_3d_c.nc", ok)
   call sht_grid_destroy(sht)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: a restart with 3-D viscosity continues bit for bit, also after the 1-D pre-spin-up'
      call radial_fe_finalize()
   else
      write(*,'(a)') ' FAIL: restart with 3-D viscosity'
      call radial_fe_finalize()
      error stop 1
   end if

contains

   include "synth_visc_fesmdata.inc"

   function same(x, y) result(r)
      real(wp), intent(in) :: x(:,:), y(:,:)
      logical :: r
      r = all(x == y)
   end function same

   subroutine restart_case(label, spin, scheme, file, ok)
      character(len=*), intent(in)    :: label, scheme, file
      logical,          intent(in)    :: spin
      logical,          intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth) :: a, c
      real(wp), allocatable :: rsl(:,:,:), dzb(:,:,:), dzs(:,:,:), zb(:,:,:), cc(:,:,:)
      real(wp) :: sc(4,K1+K2), t1, f
      integer  :: k, nbad
      logical  :: rates_ok, mem_ok

      p%lmax = LMAX;  p%nlat = 2*LMAX;  p%nphi = 4*LMAX
      p%earth_response = "ve";  p%rotation = .true.
      p%l_visc_3d = .true.;  p%visc_3d_file = VISC;  p%visc_3d_structure = 6;  p%l_toroidal = .true.
      p%pre_spinup_1d = spin;  p%equil_time_max = 0.0_wp;  p%scheme = scheme
      allocate(rsl(sht%nphi,sht%nlat,K1+K2), dzb(sht%nphi,sht%nlat,K1+K2), dzs(sht%nphi,sht%nlat,K1+K2), &
               zb(sht%nphi,sht%nlat,K1+K2), cc(sht%nphi,sht%nlat,K1+K2))

      ! uninterrupted run, a snapshot after K1 steps
      a%par = p;  call solid_earth_init(a, z_bed_eq, h_ice_eq)
      if (spin) call solid_earth_spinup(a, h_ice)
      do k = 1, K1+K2
         f = 1.0_wp - 0.1_wp*real(k, wp)
         call solid_earth_update(a, f*h_ice, 1.0e3_wp)
         rsl(:,:,k) = a%rsl;  dzb(:,:,k) = a%dz_bed;  dzs(:,:,k) = a%dz_ss
         zb(:,:,k) = a%z_bed;  cc(:,:,k) = a%gg%C
         sc(:,k) = [a%bsl, a%dphi, a%fill_int, a%worst_mass_resid]
         if (k == K1) then
            t1 = a%time
            call vilma_restart_write(a, t1, filename=file, init=.true.)
         end if
      end do

      ! a fresh model (no spin-up) restored at K1, then the same K2 steps
      c%par = p;  call solid_earth_init(c, z_bed_eq, h_ice_eq)
      call vilma_restart_read(c, file, time=t1)
      rates_ok = all(c%resp%MkPerDt == a%resp%MkPerDt) .and. all(c%resp%MkPerDt3 == a%resp%MkPerDt3) &
                 .and. all(c%resp%active1d .eqv. a%resp%active1d)
      nbad = 0
      do k = K1+1, K1+K2
         f = 1.0_wp - 0.1_wp*real(k, wp)
         call solid_earth_update(c, f*h_ice, 1.0e3_wp)
         if (.not. (same(c%rsl, rsl(:,:,k)) .and. same(c%dz_bed, dzb(:,:,k)) .and. same(c%dz_ss, dzs(:,:,k)) &
                    .and. same(c%z_bed, zb(:,:,k)) .and. same(c%gg%C, cc(:,:,k)) &
                    .and. all([c%bsl, c%dphi, c%fill_int, c%worst_mass_resid] == sc(:,k)))) nbad = nbad + 1
      end do
      mem_ok = all(c%resp%Are == a%resp%Are) .and. all(c%resp%Bim == a%resp%Bim) .and. all(c%resp%Cre == a%resp%Cre) &
               .and. (c%rotation%m == a%rotation%m)

      write(*,'(a,a,a,i0,a,i0)') ' restart [', label, ']: lmax=', LMAX, '   channels=', a%resp%nlam
      write(*,'(a,es10.2,a,es10.2,a,es10.2)') '     last step: max|rsl|', maxval(abs(c%rsl - a%rsl)), &
           '   max|dz_bed|', maxval(abs(c%dz_bed - a%dz_bed)), '   |bsl|', abs(c%bsl - a%bsl)
      write(*,'(a,i0,a,l1,a,l1)') '     steps differing after the restore: ', nbad, '   viscosity rates identical: ', rates_ok, &
           '   memory identical: ', mem_ok
      if (a%resp%nlam /= 6) then
         write(*,'(a)') '     FAIL: the toroidal channels are not on, the case is not 3-D';  ok = .false.
      end if
      if (nbad /= 0 .or. .not. rates_ok .or. .not. mem_ok) then
         write(*,'(a)') '     FAIL: the restarted run does not continue the uninterrupted one bit for bit';  ok = .false.
      end if
      call solid_earth_finalize(a);  call solid_earth_finalize(c)
   end subroutine restart_case

   subroutine make_fields(z_bed_eq, h_ice_eq, h_ice)
      !! Polar land cap (colat<50°, +500 m) over deep ocean (−4000 m); no reference ice;
      !! a 2 km grounded ice load on colat<30°.
      real(wp), intent(out) :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
            h_ice(i,j)    = merge(2000.0_wp, 0.0_wp,     thd < 30.0_wp)
         end do
      end do
      h_ice_eq = 0.0_wp
   end subroutine make_fields

end program test_restart_3d
