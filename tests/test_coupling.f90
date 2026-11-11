program test_coupling
   !! CLIMBER-X coupling contract (design.md item 1): init → update → finalize.
   !! Drives the solid_earth derived type with an ice load and checks that
   !!
   !!   (0) at the reference ice (h_ice = h_ice_eq) nothing moves: rsl ≈ 0 and
   !!       z_bed ≈ z_bed_eq (the supplied reference IS the relaxed state);
   !!   (1) ocean mass is conserved every internal Maxwell sub-step;
   !!   (2) a steady grounded-ice load relaxes sensibly: the bed subsides under
   !!       the ice immediately (elastic) and keeps subsiding with decreasing
   !!       increments toward an isostatic limit (viscoelastic relaxation), and
   !!       the ocean draws down (adding land ice removes ocean water).
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, sec_per_year
   use vilma_params,          only: vilma_param_class
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_coupling,        only: solid_earth_finalize, solid_earth_update, solid_earth_init, solid_earth
   use vilma_response,        only: response_radial_rate
   implicit none

   integer, parameter :: LMAX = 16, NSTEP = 15
   character(len=*), parameter :: VISC = "libvilma/include/test_coupling_visc.nc"
   type(sht_grid), target :: sht
   type(vilma_param_class)   :: p
   type(solid_earth)      :: se
   real(wp), allocatable  :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
   real(wp) :: dt_couple, bed_eq_ice, bed_prev, bed_now, d_first, d_last
   real(wp) :: ocean_rsl, worst_mass, net_subs
   integer  :: i, jice, jocean, step
   logical  :: ok, monotone

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)

   allocate(z_bed_eq(sht%nphi,sht%nlat), h_ice_eq(sht%nphi,sht%nlat), &
            h_ice(sht%nphi,sht%nlat))
   call make_reference(z_bed_eq, h_ice_eq)

   ! representative latitude rows: under the ice cap (colat ~ 15°) and open ocean
   jice   = nearest_row(15.0_wp)
   jocean = nearest_row(70.0_wp)
   bed_eq_ice = z_bed_eq(1,jice)

   dt_couple    = 2.0e3_wp                          ! interval per se%update call [years]
   p%lmax = LMAX;  p%nlat = 2*LMAX;  p%nphi = 4*LMAX ! model builds its own grid from par (matches local sht)
   se%par = p; call solid_earth_init(se, z_bed_eq, h_ice_eq)   ! no grid => Gauss-native (passthrough)

   write(*,'(a,i0,a)') ' coupling: lmax=', LMAX, &
                       '  (dt_couple=2 kyr, default fe scheme)'

   ! --- (0) reference consistency: no ice change -> no motion -----------------
   call solid_earth_update(se, h_ice_eq, dt_couple)
   write(*,'(a,es10.2,a,es10.2)') '   ref: max|rsl|=', maxval(abs(se%rsl)), &
                                  '   max|z_bed-z_bed_eq|=', maxval(abs(se%z_bed - z_bed_eq))
   if (maxval(abs(se%rsl)) > 1.0e-9_wp .or. maxval(abs(se%z_bed - z_bed_eq)) > 1.0e-9_wp) then
      write(*,'(a)') '   FAIL: reference state is not stationary'; ok = .false.
   end if

   ! --- (1)+(2) hold a 2 km grounded ice cap and watch it relax ---------------
   call make_load(h_ice)                           ! cap on colat<30 (on land)
   write(*,'(a)') '       step   z_bed under ice [m]   d(subs) [m]   worst mass resid'
   bed_prev = bed_eq_ice;  worst_mass = 0.0_wp;  monotone = .true.
   do step = 1, NSTEP
      call solid_earth_update(se, h_ice, dt_couple)
      bed_now = se%z_bed(1,jice)
      worst_mass = max(worst_mass, se%worst_mass_resid)
      if (step == 1)     d_first = bed_prev - bed_now
      if (step == NSTEP) d_last  = bed_prev - bed_now
      if (bed_now > bed_prev + 1.0e-6_wp) monotone = .false.   ! must not rebound
      if (mod(step,3) == 1 .or. step == NSTEP) &
         write(*,'(i9,f18.4,f15.4,es18.2)') step, bed_now, bed_prev-bed_now, &
                                            se%worst_mass_resid
      bed_prev = bed_now
   end do

   ocean_rsl = se%rsl(1,jocean)
   net_subs  = bed_eq_ice - se%z_bed(1,jice)

   write(*,'(a)') ''
   write(*,'(a,f12.4,a)')  '      net bed subsidence under ice =', net_subs, ' m'
   write(*,'(a,f12.4,a,f12.4,a)') '      subsidence increment: first =', d_first, &
                                  '   last =', d_last, ' m'
   write(*,'(a,f12.4,a)')  '      ocean rsl (drawdown, should be < 0) =', ocean_rsl, ' m'
   write(*,'(a,es11.2)')   '      worst ocean-mass residual over the run =', worst_mass

   ! --- assertions ------------------------------------------------------------
   if (worst_mass > 1.0e-10_wp) then
      write(*,'(a)') '      FAIL: ocean mass not conserved during stepping'; ok = .false.
   end if
   if (.not. monotone) then
      write(*,'(a)') '      FAIL: bed rebounds under a steady load'; ok = .false.
   end if
   if (net_subs < 10.0_wp) then
      write(*,'(a)') '      FAIL: no appreciable subsidence under the ice'; ok = .false.
   end if
   if (d_last >= d_first) then
      write(*,'(a)') '      FAIL: subsidence not slowing (no relaxation toward a limit)'
      ok = .false.
   end if
   if (ocean_rsl >= 0.0_wp) then
      write(*,'(a)') '      FAIL: building land ice should draw the ocean down'; ok = .false.
   end if

   call solid_earth_finalize(se)

   ! --- (3) rotational feedback through the coupling driver (off-axis load) ----
   ! An off-axis ice cap drives polar motion; rotation must run, conserve mass, and
   ! perceptibly change the bed vs the rotation-off run (the centrifugal sea-level
   ! feedback). On-axis loads (the cap above) carry no (2,1) inertia, so use a load
   ! centred away from the pole here.
   call rotation_check()
   call rotation_visc_check()
   call reinit_check()
   call dt_min_check()
   call elastic_rotation_check()
   call sle_cap_check()

   call sht_grid_destroy(sht)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: coupling init/update/finalize conserves ocean mass,'
      write(*,'(a)') '       holds the reference state, relaxes a load, + rotation'
   else
      write(*,'(a)') ' FAIL: coupling validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   include "synth_visc_fesmdata.inc"

   subroutine make_reference(z_bed_eq, h_ice_eq)
      !! Reference (equilibrium) topography: a polar land cap (colat<50°, +500 m)
      !! over deep ocean (−4000 m); no reference ice.
      real(wp), intent(out) :: z_bed_eq(:,:), h_ice_eq(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
         end do
      end do
      h_ice_eq = 0.0_wp
   end subroutine make_reference

   subroutine make_load(h_ice)
      !! A 2 km grounded ice cap on the land interior (colat<30°), held fixed.
      real(wp), intent(out) :: h_ice(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            h_ice(i,j) = merge(2000.0_wp, 0.0_wp, thd < 30.0_wp)
         end do
      end do
   end subroutine make_load

   subroutine rotation_check()
      !! Drive the coupling with an off-axis cap, rotation OFF then ON. Assert: the
      !! rotation-on run conserves ocean mass, develops a finite polar motion, and the
      !! bed differs from the rotation-off run (the feedback has an effect).
      type(vilma_param_class) :: pr
      type(solid_earth)    :: se_off, se_on
      real(wp), allocatable :: h_off(:,:), zoff(:,:)
      real(wp) :: wmass, mdeg, dzmax
      integer  :: s
      allocate(h_off(sht%nphi,sht%nlat), zoff(sht%nphi,sht%nlat))
      call make_load_offaxis(h_off)

      pr%lmax = LMAX;  pr%nlat = 2*LMAX;  pr%nphi = 4*LMAX
      pr%rotation  = .false.
      se_off%par = pr; call solid_earth_init(se_off, z_bed_eq, h_ice_eq)
      do s = 1, NSTEP;  call solid_earth_update(se_off, h_off, dt_couple);  end do
      zoff = se_off%z_bed
      call solid_earth_finalize(se_off)

      pr%rotation = .true.
      se_on%par = pr; call solid_earth_init(se_on, z_bed_eq, h_ice_eq)
      wmass = 0.0_wp
      do s = 1, NSTEP
         call solid_earth_update(se_on, h_off, dt_couple)
         wmass = max(wmass, se_on%worst_mass_resid)
      end do
      mdeg  = abs(se_on%rotation%m)*180.0_wp/(acos(-1.0_wp))
      dzmax = maxval(abs(se_on%z_bed - zoff))
      write(*,'(a)') ''
      write(*,'(a,f10.5,a,es10.2,a,f9.4,a)') '   rotation: |m|=', mdeg, &
           ' deg   worst mass=', wmass, '   max|Δz_bed on-off|=', dzmax, ' m'
      if (wmass > 1.0e-10_wp) then
         write(*,'(a)') '   FAIL: rotation broke ocean-mass conservation'; ok = .false.
      end if
      if (mdeg <= 0.0_wp .or. mdeg /= mdeg) then
         write(*,'(a)') '   FAIL: no (finite) polar motion under an off-axis load'; ok = .false.
      end if
      if (dzmax < 1.0e-3_wp) then
         write(*,'(a)') '   FAIL: rotational feedback had no effect on the bed'; ok = .false.
      end if
      call solid_earth_finalize(se_on)
   end subroutine rotation_check

   subroutine rotation_visc_check()
      !! With a viscosity file, the rotation channels must relax on the response's
      !! radial profile (response_radial_rate), not the layer table: equal rates in
      !! both channels, and at a genuinely 3-D element (where the response keeps the
      !! layer-table scalar) the lateral geometric mean, not that scalar.
      type(vilma_param_class) :: pr
      type(solid_earth)    :: se
      real(wp), allocatable :: rate(:)
      real(wp) :: dmax
      integer  :: e
      call write_synth_visc_fesmdata(VISC)
      pr%lmax = LMAX;  pr%nlat = 2*LMAX;  pr%nphi = 4*LMAX
      pr%rotation = .true.;  pr%l_visc_3d = .true.;  pr%visc_3d_file = VISC;  pr%visc_3d_structure = 6
      se%par = pr; call solid_earth_init(se, z_bed_eq, h_ice_eq)
      rate  = response_radial_rate(se%resp)
      dmax  = max(maxval(abs(se%rotation%load_ch%MkPerDt - rate)), &
                  maxval(abs(se%rotation%tidal_ch%MkPerDt - rate)))
      write(*,'(a)') ''
      write(*,'(a,es10.2,a,i0,a,i0)') '   rotation viscosity: max|rate_rot − rate_resp| = ', dmax, &
           '   3-D elements ', se%resp%ne3d, ' of ', se%resp%ne
      if (dmax /= 0.0_wp) then
         write(*,'(a)') '   FAIL: rotation channels do not relax on the response viscosity'; ok = .false.
      end if
      if (se%resp%ne3d == 0) then
         write(*,'(a)') '   FAIL: the Bagge field gave no 3-D element to test'; ok = .false.
      else
         e = se%resp%e3d(1)
         write(*,'(a,i0,a,es10.3,a,es10.3)') '   element ', e, ': rotation rate ', rate(e), &
              '   layer-table rate ', se%resp%MkPerDt(e)
         if (rate(e) == se%resp%MkPerDt(e)) then
            write(*,'(a)') '   FAIL: 3-D element carries the layer-table rate'; ok = .false.
         end if
      end if
      call solid_earth_finalize(se)
   end subroutine rotation_visc_check

   subroutine reinit_check()
      !! init -> update -> finalize -> init -> update on ONE object, with rotation and a
      !! viscosity file (the configuration in which a leftover rotation flag stopped the
      !! second init). Finalize must leave the object pristine, so the second run
      !! repeats the first bit for bit.
      type(vilma_param_class) :: pr
      type(solid_earth)    :: se
      real(wp), allocatable :: h_off(:,:), z1(:,:)
      complex(wp) :: m1
      integer :: run, s
      allocate(h_off(sht%nphi,sht%nlat))
      call make_load_offaxis(h_off)
      pr%lmax = LMAX;  pr%nlat = 2*LMAX;  pr%nphi = 4*LMAX
      pr%rotation = .true.;  pr%l_visc_3d = .true.;  pr%visc_3d_file = VISC;  pr%visc_3d_structure = 6
      write(*,'(a)') ''
      do run = 1, 2
         se%par = pr; call solid_earth_init(se, z_bed_eq, h_ice_eq)
         do s = 1, 2;  call solid_earth_update(se, h_off, 500.0_wp);  end do
         if (run == 1) then
            z1 = se%z_bed;  m1 = se%rotation%m
         else
            write(*,'(a,es10.2,a,es10.2)') '   re-init: max|z_bed run2 − run1| = ', maxval(abs(se%z_bed - z1)), &
                 '   |m run2 − run1| = ', abs(se%rotation%m - m1)
            if (any(se%z_bed /= z1) .or. se%rotation%m /= m1) then
               write(*,'(a)') '   FAIL: a re-initialised object does not repeat the first run'; ok = .false.
            end if
         end if
         call solid_earth_finalize(se)
         if (se%rotation%enabled .or. se%use_vilma_v1 .or. se%remap .or. associated(se%sht) &
             .or. allocated(se%rsl) .or. allocated(se%gg%h_ice) .or. allocated(se%sigma_lm) &
             .or. allocated(se%rotation%load_ch%Are) .or. se%time /= 0.0_wp &
             .or. se%stepper%n_accept /= 0 .or. se%sle%n_solve /= 0 .or. se%bsl /= 0.0_wp) then
            write(*,'(a)') '   FAIL: solid_earth_finalize left state behind'; ok = .false.
         end if
      end do
   end subroutine reinit_check

   subroutine dt_min_check()
      !! The explicit sub-step is the stability ceiling cfl/max(mu/eta) (278 yr for the
      !! default earth: 8 sub-steps in 2 kyr). A dt_min above it overrides it and must be
      !! flagged; one below it changes nothing.
      type(vilma_param_class) :: pr
      type(solid_earth)    :: se
      real(wp), allocatable :: h(:,:)
      integer :: pass, nsub(2)
      logical :: flag(2)
      allocate(h(sht%nphi,sht%nlat));  call make_load(h)
      pr%lmax = LMAX;  pr%nlat = 2*LMAX;  pr%nphi = 4*LMAX;  pr%rotation = .false.
      do pass = 1, 2
         pr%dt_min = merge(100.0_wp, 400.0_wp, pass == 1)*sec_per_year
         se%par = pr; call solid_earth_init(se, z_bed_eq, h_ice_eq)
         call solid_earth_update(se, h, dt_couple)
         nsub(pass) = se%stepper%n_accept;  flag(pass) = se%stepper%dt_min_overrides
         call solid_earth_finalize(se)
      end do
      write(*,'(a)') ''
      write(*,'(a,i0,a,l1,a,i0,a,l1)') '   dt_min 100 yr: ', nsub(1), ' sub-steps, overrides=', flag(1), &
           '   dt_min 400 yr: ', nsub(2), ' sub-steps, overrides=', flag(2)
      if (flag(1) .or. .not. flag(2) .or. nsub(1) /= 8 .or. nsub(2) /= 5) then
         write(*,'(a)') '   FAIL: dt_min over the stability sub-step is not flagged (or the sub-steps changed)'
         ok = .false.
      end if
   end subroutine dt_min_check

   subroutine elastic_rotation_check()
      !! earth_response = "elastic" with the rotation on (the default) must initialise and
      !! run: the rotation channels are then elastic as well (no relaxation), so the polar
      !! motion under a held load is the same after every interval.
      type(vilma_param_class) :: pr
      type(solid_earth)    :: se
      real(wp), allocatable :: h_off(:,:)
      complex(wp) :: m1, m2
      allocate(h_off(sht%nphi,sht%nlat));  call make_load_offaxis(h_off)
      pr%lmax = LMAX;  pr%nlat = 2*LMAX;  pr%nphi = 4*LMAX
      pr%earth_response = "elastic";  pr%rotation = .true.
      se%par = pr; call solid_earth_init(se, z_bed_eq, h_ice_eq)
      call solid_earth_update(se, h_off, dt_couple);  m1 = se%rotation%m
      call solid_earth_update(se, h_off, dt_couple);  m2 = se%rotation%m
      write(*,'(a)') ''
      write(*,'(a,es10.3,a,es10.2)') '   elastic + rotation: |m| = ', abs(m1), ' rad   |m(2) - m(1)| = ', abs(m2 - m1)
      if (.not. (abs(m1) > 0.0_wp) .or. abs(m2 - m1) > 1.0e-5_wp*abs(m1) &
          .or. maxval(se%rotation%load_ch%MkPerDt) /= 0.0_wp) then
         write(*,'(a)') '   FAIL: an elastic Earth with rotation does not give a steady polar motion'; ok = .false.
      end if
      call solid_earth_finalize(se)
   end subroutine elastic_rotation_check

   subroutine sle_cap_check()
      !! The convergence record of the sea-level solves: with the defaults none ends at a
      !! cap and the residual is within sle_tol; with one inner iteration every solve ends
      !! at that cap and the residual says how far it is from the fixed point.
      type(vilma_param_class) :: pr
      type(solid_earth)    :: se
      real(wp), allocatable :: h(:,:)
      integer :: pass
      allocate(h(sht%nphi,sht%nlat));  call make_load(h)
      pr%lmax = LMAX;  pr%nlat = 2*LMAX;  pr%nphi = 4*LMAX;  pr%rotation = .false.
      write(*,'(a)') ''
      do pass = 1, 2
         if (pass == 2) pr%sle_n_inner = 1
         se%par = pr; call solid_earth_init(se, z_bed_eq, h_ice_eq)
         call solid_earth_update(se, h, dt_couple)
         write(*,'(a,i0,a,es10.2,a,i0,a,i0,a,i0,a,i0)') '   sle_n_inner=', se%sle%n_inner, ': resid ', &
              se%sle_resid, '   solves ', se%sle_n_solve, '   capped inner/coast/mem ', se%sle_n_cap_inner, &
              '/', se%sle_n_cap_coast, '/', se%sle_n_cap_mem
         if (pass == 1 .and. (se%sle_n_cap_inner + se%sle_n_cap_coast + se%sle_n_cap_mem /= 0 &
                              .or. se%sle_resid > pr%sle_tol .or. se%sle_n_solve < 1)) then
            write(*,'(a)') '   FAIL: a default solve reports a cap or a residual above sle_tol'; ok = .false.
         end if
         if (pass == 2 .and. (se%sle_n_cap_inner < 1 .or. se%sle_resid <= pr%sle_tol)) then
            write(*,'(a)') '   FAIL: a solve cut at one inner iteration is not reported as capped'; ok = .false.
         end if
         call solid_earth_finalize(se)
      end do
   end subroutine sle_cap_check

   subroutine make_load_offaxis(h_ice)
      !! A 2 km grounded cap centred at colat 35°, lon 90° (off the rotation axis, on
      !! land), radius ~20°.
      real(wp), intent(out) :: h_ice(:,:)
      real(wp) :: tc, lc, cg, ca
      integer  :: i, j
      tc = 35.0_wp*pi/180.0_wp;  lc = 90.0_wp*pi/180.0_wp;  ca = cos(20.0_wp*pi/180.0_wp)
      h_ice = 0.0_wp
      do j = 1, sht%nlat
         do i = 1, sht%nphi
            cg = cos(tc)*cos(sht%colat(j)) + sin(tc)*sin(sht%colat(j))*cos(sht%lon(i)-lc)
            if (cg >= ca) h_ice(i,j) = 2000.0_wp
         end do
      end do
   end subroutine make_load_offaxis

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

end program test_coupling
