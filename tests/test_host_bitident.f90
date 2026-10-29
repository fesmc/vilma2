program test_host_bitident
   !! Host terms that are zero change nothing: a run given w = 0 and zero I_g, C, F at
   !! init, update and spin-up is bit-for-bit (==) the run given no host arguments, because
   !! the coupling then dispatches to the host-free statements. Two runs are stepped in
   !! lockstep and every field compared after every stage: rsl, dz_bed, dz_ss, C, C0, the
   !! bed and ice, bsl, bsl_ive, dphi, the spectral load (the memory's image) and the
   !! memory norm.
   !!
   !! Scenario of test_bsl (polar land cap, a shelf ring that floods, ice on the cap and the
   !! shelf, then removed), lmax 16, rotation on; forward Euler and trapezoidal; init ->
   !! dt = 0 seed -> updates -> spin-up (1-D and full) -> update.
   use vilma_precision, only: wp
   use vilma_constants, only: pi, sec_per_year
   use vilma_params,    only: vilma_param_class
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_response,  only: response_memory_norm
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_spinup, &
                              solid_earth_finalize
   implicit none

   integer, parameter :: LMAX = 16
   real(wp), allocatable :: z_eq(:,:), h_eq(:,:), h_on(:,:), h_off(:,:), zero(:,:), w0(:,:)
   logical :: ok
   integer :: np, nl, i, j

   ok = .true.
   np = 4*LMAX;  nl = 2*LMAX + 2
   allocate(z_eq(np,nl), h_eq(np,nl), h_on(np,nl), h_off(np,nl), zero(np,nl), w0(np,nl))
   call fields()
   zero = 0.0_wp;  w0 = 0.0_wp

   call run("fe", ok)
   call run("trap", ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: zero host terms are bit-identical to no host terms'
      call radial_fe_finalize()
   else
      write(*,'(a)') ' FAIL: zero host terms differ from no host terms'
      call radial_fe_finalize()
      error stop 1
   end if

contains

   subroutine fields()
      real(wp) :: th
      do j = 1, nl
         do i = 1, np
            th = colat(j)
            if (th < 50.0_wp) then
               z_eq(i,j) = 500.0_wp
            else if (th < 62.0_wp) then
               z_eq(i,j) = -100.0_wp
            else
               z_eq(i,j) = -4000.0_wp
            end if
            h_on(i,j) = 0.0_wp
            if (th < 30.0_wp) h_on(i,j) = 2000.0_wp
            if (th >= 50.0_wp .and. th < 62.0_wp) h_on(i,j) = 300.0_wp
         end do
      end do
      h_eq = 0.0_wp;  h_off = 0.0_wp
   end subroutine fields

   real(wp) function colat(j)
      !! Gauss-Legendre colatitude [deg] of row j, from the grid the model builds.
      use vilma_sht, only: sht_grid, sht_grid_init, sht_grid_destroy
      integer, intent(in) :: j
      type(sht_grid), save :: g
      logical, save :: built = .false.
      if (.not. built) then
         call sht_grid_init(g, LMAX, nlat=nl, nphi=np);  built = .true.
      end if
      colat = g%colat(j)*180.0_wp/pi
   end function colat

   subroutine run(scheme, ok)
      character(len=*), intent(in)    :: scheme
      logical,          intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth)       :: a, b                 ! a: no host arguments; b: zero host terms
      integer :: nfail
      nfail = 0
      p%lmax = LMAX;  p%nlat = nl;  p%nphi = np
      p%scheme = scheme;  p%rotation = .true.
      a%par = p;  b%par = p
      write(*,'(a,a)') ' scheme = ', scheme

      call solid_earth_init(a, z_eq, h_eq)
      call solid_earth_init(b, z_eq, h_eq, w=w0, I_g=zero, C=zero)
      call compare(a, b, 'init', nfail)

      call solid_earth_update(a, h_on, 0.0_wp)
      call solid_earth_update(b, h_on, 0.0_wp, w0, zero, zero, zero)
      call compare(a, b, 'update dt=0', nfail)
      do i = 1, 4
         call solid_earth_update(a, merge(h_on, h_off, i <= 2), 500.0_wp)
         call solid_earth_update(b, merge(h_on, h_off, i <= 2), 500.0_wp, w0, zero, zero, zero)
         call compare(a, b, 'update', nfail)
      end do

      ! spin-up under the load: 1-D pre-equilibration, then the full model
      a%par%pre_spinup_1d = .true.;  a%par%equil_time_max = 2.0e4_wp*sec_per_year
      b%par%pre_spinup_1d = .true.;  b%par%equil_time_max = 2.0e4_wp*sec_per_year
      call solid_earth_spinup(a, h_on)
      call solid_earth_spinup(b, h_on, w0, zero, zero, zero)
      call compare(a, b, 'spinup', nfail)
      call solid_earth_update(a, h_off, 300.0_wp)
      call solid_earth_update(b, h_off, 300.0_wp, w0, zero, zero, zero)
      call compare(a, b, 'update after spinup', nfail)

      if (nfail > 0) ok = .false.
      call solid_earth_finalize(a);  call solid_earth_finalize(b)
   end subroutine run

   subroutine compare(a, b, stage, nfail)
      type(solid_earth),  intent(in)    :: a, b
      character(len=*),   intent(in)    :: stage
      integer,            intent(inout) :: nfail
      logical :: same
      same = all(a%gg%rsl == b%gg%rsl) .and. all(a%gg%dz_bed == b%gg%dz_bed) .and. all(a%gg%dz_ss == b%gg%dz_ss) &
             .and. all(a%gg%C == b%gg%C) .and. all(a%gg%C0 == b%gg%C0) .and. all(a%gg%z_bed == b%gg%z_bed) &
             .and. all(a%gg%h_ice == b%gg%h_ice) .and. all(a%rsl == b%rsl) .and. all(a%z_bed == b%z_bed) &
             .and. a%bsl == b%bsl .and. a%bsl_ive == b%bsl_ive .and. a%dphi == b%dphi &
             .and. a%worst_mass_resid == b%worst_mass_resid .and. a%time == b%time &
             .and. all(a%sigma_lm == b%sigma_lm) &
             .and. response_memory_norm(a%resp) == response_memory_norm(b%resp)
      write(*,'(a,a,a,l1,a,es10.2)') '   ', stage, ': identical = ', same, '   max|rsl| = ', maxval(abs(a%gg%rsl))
      if (.not. same) then
         write(*,'(a)') '     FAIL: fields differ';  nfail = nfail + 1
      end if
   end subroutine compare

end program test_host_bitident
