program test_host_terms
   !! Host-resolved terms through the coupling API (solid_earth_init / update with w, I_g, C,
   !! F), lmax 16, 1-D PREM-like default Earth, forward Euler, rotation on.
   !!
   !!   (1) solid_earth_host_check accepts admissible terms and rejects each violation;
   !!   (2) ramp: one update of 800 yr (40 sub-steps of 20 yr) == four updates of 200 yr with
   !!       the ice and the host terms interpolated to their ends;
   !!   (3) ocean mass balance with 0 < w < 1, own ice undefined (NaN) where w = 1:
   !!       ∫C·rsl = -(rho_i/rho_w) ∫ΔI_g + ∫fill, and bsl, bsl_ive tied by ∫fill / ∫C;
   !!   (4) passthrough == a global lon-lat host: the same terms on a 96 x 48 host grid, remapped
   !!       by the model, and handed over already remapped on the Gauss grid;
   !!   (5) two hosts (disjoint regions) summed == one host over both;
   !!   (6) round-off beyond [0,1] and [0,w] is clipped, for a Gauss host as for a remapped one..
   use vilma_precision, only: wp
   use vilma_constants, only: pi, sec_per_year, rho_ice, rho_water
   use vilma_params,    only: vilma_param_class
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_surface_integral
   use vilma_sle,       only: ocean_function
   use vilma_remap,     only: remap_ll_gauss, remap_init, remap_to_gauss
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize, &
                              solid_earth_host_check
   use coords,          only: grid_class, grid_init
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none

   integer,  parameter :: LMAX = 16
   integer,  parameter :: NLON = 96, NLAT = 48
   real(wp), parameter :: DEG = pi/180.0_wp
   type(sht_grid), target :: sht
   integer  :: np, nl
   real(wp), allocatable :: colat(:,:), lon(:,:)
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX+2, nphi=4*LMAX)
   np = sht%nphi;  nl = sht%nlat
   allocate(colat(np,nl), lon(np,nl))
   block
      integer :: i, j
      do j = 1, nl
         do i = 1, np
            colat(i,j) = sht%colat(j)/DEG;  lon(i,j) = sht%lon(i)/DEG
         end do
      end do
   end block

   call check_host_check(ok)
   call check_ramp(ok)
   call check_balance(ok)
   call check_passthrough(ok)
   call check_two_hosts(ok)
   call check_clip(ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: host terms through the coupling: checks, ramp, mass balance, remap, sums'
      call radial_fe_finalize()
   else
      write(*,'(a)') ' FAIL: host terms through the coupling'
      call radial_fe_finalize()
      error stop 1
   end if

contains

   ! --- scene ---------------------------------------------------------------

   subroutine scene(th, lo, z_eq, h1, w, I1, C1, F1, I0, C0, Cr, wsel)
      !! Bed (land cap, shelf, bare shallow ring, deep ocean), own ice at t1 (2 km cap, 300 m
      !! shelf), and host terms on the grid with colat th, lon lo [deg]: region w = 1 (lon < 90)
      !! and 1/2 (lon < 180) in colat 20..70 (wsel: only the lon-sectors given: 1 = [0,90),
      !! 2 = [90,180), 3 = [180,270), 0 = all of the first two). The host's coast at t1 is
      !! colat 54, at t0 58 (its ice a quarter, half the weight: I0, C0), and the reference's
      !! 58 (Cr, ice-free, full weight).
      real(wp), intent(in)  :: th(:,:), lo(:,:)
      real(wp), intent(out) :: z_eq(:,:), h1(:,:), w(:,:), I1(:,:), C1(:,:), F1(:,:), I0(:,:), C0(:,:), Cr(:,:)
      integer,  intent(in)  :: wsel
      integer  :: i, j
      real(wp) :: ww, hh, fg, fr
      do j = 1, size(th,2)
         do i = 1, size(th,1)
            if (th(i,j) < 50.0_wp) then
               z_eq(i,j) = 500.0_wp
            else if (th(i,j) < 62.0_wp) then
               z_eq(i,j) = -100.0_wp
            else if (th(i,j) < 70.0_wp) then
               z_eq(i,j) = -60.0_wp
            else
               z_eq(i,j) = -4000.0_wp
            end if
            h1(i,j) = 0.0_wp
            if (th(i,j) < 30.0_wp) h1(i,j) = 2000.0_wp
            if (th(i,j) >= 50.0_wp .and. th(i,j) < 62.0_wp) h1(i,j) = 300.0_wp
            ww = 0.0_wp
            if (th(i,j) >= 20.0_wp .and. th(i,j) < 70.0_wp) then
               select case (wsel)
               case (0);  if (lo(i,j) < 90.0_wp) then; ww = 1.0_wp; else if (lo(i,j) < 180.0_wp) then; ww = 0.5_wp; end if
               case (1);  if (lo(i,j) < 90.0_wp) ww = 1.0_wp
               case (2);  if (lo(i,j) >= 90.0_wp .and. lo(i,j) < 180.0_wp) ww = 0.5_wp
               case (3);  if (lo(i,j) >= 180.0_wp .and. lo(i,j) < 270.0_wp) ww = 0.75_wp
               case (4);  if (lo(i,j) < 90.0_wp) then; ww = 1.0_wp; else if (lo(i,j) < 180.0_wp) then; ww = 0.5_wp; &
                          else if (lo(i,j) < 270.0_wp) then; ww = 0.75_wp; end if
               end select
            end if
            w(i,j) = ww
            hh = merge(3000.0_wp, 200.0_wp, th(i,j) < 30.0_wp)
            fg = merge(1.0_wp, 0.0_wp, th(i,j) < 54.0_wp)
            fr = merge(1.0_wp, 0.0_wp, th(i,j) < 58.0_wp)
            I1(i,j) = ww*hh*fg
            C1(i,j) = ww*(1.0_wp - fg)
            Cr(i,j) = ww*(1.0_wp - fr)
            F1(i,j) = ww*merge(200.0_wp, -50.0_wp, fg > 0.5_wp)*((1.0_wp - fg) - (1.0_wp - fr))
            I0(i,j) = 0.5_wp*ww*0.25_wp*hh*fr
            C0(i,j) = 0.5_wp*ww*(1.0_wp - fr)
         end do
      end do
   end subroutine scene

   subroutine new_par(p, dtmax_yr)
      type(vilma_param_class), intent(out) :: p
      real(wp), intent(in), optional :: dtmax_yr
      p%lmax = LMAX;  p%nlat = 2*LMAX+2;  p%nphi = 4*LMAX
      if (present(dtmax_yr)) p%dt_max = dtmax_yr*sec_per_year
   end subroutine new_par

   real(wp) function diff(a, b)
      real(wp), intent(in) :: a(:,:), b(:,:)
      diff = maxval(abs(a - b))
   end function diff

   ! --- (1) host_check --------------------------------------------------------

   subroutine check_host_check(ok)
      logical, intent(inout) :: ok
      real(wp), allocatable :: w(:,:), I_g(:,:), C(:,:), F(:,:), W2(:,:), x(:,:)
      character(len=120) :: msg
      integer :: nfail
      nfail = 0
      allocate(w(4,3), I_g(4,3), C(4,3), F(4,3), W2(4,3), x(4,3))
      w = 0.5_wp;  I_g = 100.0_wp;  C = 0.25_wp;  F = -3.0_wp
      w(1,1) = 0.0_wp;  I_g(1,1) = 0.0_wp;  C(1,1) = 0.0_wp;  F(1,1) = 0.0_wp
      write(*,'(a)') ' (1) solid_earth_host_check'
      call expect(w, I_g, C, F, .true., 'admissible terms', nfail)
      C(2,2) = 0.5_wp + 0.5e-9_wp;   call expect(w, I_g, C, F, .true., 'C = w + 0.5e-9 (inside the tolerance)', nfail)
      C(2,2) = 0.5_wp + 2.0e-9_wp;   call expect(w, I_g, C, F, .false., 'C = w + 2e-9', nfail);       C(2,2) = 0.25_wp
      C(2,2) = -0.5e-9_wp;           call expect(w, I_g, C, F, .true.,  'C = -0.5e-9', nfail)
      C(2,2) = -2.0e-9_wp;           call expect(w, I_g, C, F, .false., 'C = -2e-9', nfail);          C(2,2) = 0.25_wp
      w(2,2) = 1.0_wp + 0.5e-9_wp;   call expect(w, I_g, C, F, .true.,  'w = 1 + 0.5e-9', nfail)
      w(2,2) = 1.0_wp + 2.0e-9_wp;   call expect(w, I_g, C, F, .false., 'w = 1 + 2e-9', nfail);       w(2,2) = 0.5_wp
      w(2,2) = -0.5e-9_wp;           C(2,2) = 0.0_wp
      call expect(w, I_g, C, F, .false., 'w = -0.5e-9 with I_g /= 0', nfail)
      I_g(2,2) = 0.0_wp;  F(2,2) = 0.0_wp
      call expect(w, I_g, C, F, .true.,  'w = -0.5e-9, terms 0 (clipped to w = 0)', nfail)
      w(2,2) = -2.0e-9_wp;           call expect(w, I_g, C, F, .false., 'w = -2e-9', nfail)
      w(2,2) = 0.5_wp;  C(2,2) = 0.25_wp;  I_g(2,2) = 100.0_wp;  F(2,2) = -3.0_wp
      I_g(3,3) = -1.0e-6_wp;         call expect(w, I_g, C, F, .false., 'I_g < 0', nfail);            I_g(3,3) = 100.0_wp
      C(3,3) = -1.0e-6_wp;           call expect(w, I_g, C, F, .false., 'C < 0', nfail);              C(3,3) = 0.25_wp
      I_g(1,1) = 1.0e-12_wp;         call expect(w, I_g, C, F, .false., 'I_g /= 0 where w = 0', nfail); I_g(1,1) = 0.0_wp
      C(1,1) = 1.0e-12_wp;           call expect(w, I_g, C, F, .true.,  'C = 1e-12 where w = 0 (clipped to 0)', nfail); C(1,1) = 0.0_wp
      C(1,1) = 1.0e-6_wp;            call expect(w, I_g, C, F, .false., 'C = 1e-6 where w = 0', nfail);  C(1,1) = 0.0_wp
      F(1,1) = 1.0e-12_wp;           call expect(w, I_g, C, F, .false., 'F /= 0 where w = 0', nfail);   F(1,1) = 0.0_wp
      F(4,2) = ieee_value(1.0_wp, ieee_quiet_nan);  call expect(w, I_g, C, F, .false., 'F = NaN', nfail);  F(4,2) = -3.0_wp
      I_g(4,2) = -ieee_value(1.0_wp, ieee_quiet_nan); call expect(w, I_g, C, F, .false., 'I_g = NaN', nfail); I_g(4,2) = 100.0_wp
      w(4,2) = ieee_value(1.0_wp, ieee_quiet_nan);  call expect(w, I_g, C, F, .false., 'w = NaN', nfail);  w(4,2) = 0.5_wp
      call expect(w, I_g, C, F, .true., 'admissible again', nfail)
      if (nfail > 0) ok = .false.

   end subroutine check_host_check

   subroutine expect(w, I_g, C, F, want, what, nfail)
      real(wp),         intent(in)    :: w(:,:), I_g(:,:), C(:,:), F(:,:)
      logical,          intent(in)    :: want
      character(len=*), intent(in)    :: what
      integer,          intent(inout) :: nfail
      character(len=120) :: msg
      logical :: got
      got = solid_earth_host_check(w, I_g, C, F, msg)
      if (got .neqv. want) then
         write(*,'(a,a,a,l1,a,a)') '     FAIL: ', what, ': returned ', got, '  ', trim(msg);  nfail = nfail + 1
      end if
   end subroutine expect

   ! --- (2) ramp ----------------------------------------------------------------

   subroutine check_ramp(ok)
      logical, intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth) :: a, b
      real(wp), allocatable :: z_eq(:,:), h_eq(:,:), h1(:,:), w1(:,:), I1(:,:), C1(:,:), F1(:,:), I0(:,:), C0(:,:), Cr(:,:)
      real(wp), allocatable :: w0(:,:), F0(:,:), hk(:,:), wk(:,:), Ik(:,:), Ck(:,:), Fk(:,:)
      real(wp) :: f, d, dd
      integer  :: k
      allocate(z_eq(np,nl), h_eq(np,nl), h1(np,nl), w1(np,nl), I1(np,nl), C1(np,nl), F1(np,nl), I0(np,nl), &
               C0(np,nl), Cr(np,nl), w0(np,nl), F0(np,nl), hk(np,nl), wk(np,nl), Ik(np,nl), Ck(np,nl), Fk(np,nl))
      call scene(colat, lon, z_eq, h1, w1, I1, C1, F1, I0, C0, Cr, 0)
      h_eq = 0.0_wp;  w0 = 0.5_wp*w1;  F0 = 0.0_wp
      call new_par(p, 20.0_wp)
      a%par = p;  b%par = p
      ! terms at t0 (half the weight): seed the ramp start with a dt = 0 update on both
      call solid_earth_init(a, z_eq, h_eq, w=w1, I_g=0.0_wp*I1, C=Cr)
      call solid_earth_init(b, z_eq, h_eq, w=w1, I_g=0.0_wp*I1, C=Cr)
      call solid_earth_update(a, h_eq, 0.0_wp, w0, I0, C0, F0)
      call solid_earth_update(b, h_eq, 0.0_wp, w0, I0, C0, F0)
      ! a: one update; b: four, with the ice and the terms at the quarter points
      call solid_earth_update(a, h1, 800.0_wp, w1, I1, C1, F1)
      do k = 1, 4
         f = real(k,wp)/4.0_wp
         hk = f*h1;  wk = w0 + f*(w1 - w0);  Ik = I0 + f*(I1 - I0);  Ck = C0 + f*(C1 - C0);  Fk = F0 + f*(F1 - F0)
         call solid_earth_update(b, hk, 200.0_wp, wk, Ik, Ck, Fk)
      end do
      d = diff(a%gg%rsl, b%gg%rsl);  dd = diff(a%gg%dIg, b%gg%dIg)
      write(*,'(a)') ' (2) ramp: 1 x 800 yr == 4 x 200 yr at interpolated ends'
      write(*,'(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') '     |drsl| =', d, ' m   |dIg| =', dd, &
           '   |dC| =', diff(a%gg%C, b%gg%C), '   |dbsl| =', abs(a%bsl - b%bsl)
      write(*,'(a,es10.2,a,i0,a,i0)') '     max|rsl| =', maxval(abs(a%gg%rsl)), '   sub-steps ', a%stepper%n_accept, &
           ' and ', b%stepper%n_accept
      if (d > 1.0e-8_wp .or. dd > 1.0e-8_wp .or. diff(a%gg%C, b%gg%C) > 1.0e-12_wp .or. abs(a%bsl - b%bsl) > 1.0e-8_wp) then
         write(*,'(a)') '     FAIL: ramp is not the interpolation of the terms';  ok = .false.
      end if
      if (maxval(abs(a%gg%rsl)) < 1.0_wp) then
         write(*,'(a)') '     FAIL: no signal';  ok = .false.
      end if
      call solid_earth_finalize(a);  call solid_earth_finalize(b)
   end subroutine check_ramp

   ! --- (3) mass balance --------------------------------------------------------

   subroutine check_balance(ok)
      logical, intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth) :: se
      real(wp), allocatable :: z_eq(:,:), h_eq(:,:), h1(:,:), w1(:,:), I1(:,:), C1(:,:), F1(:,:), I0(:,:), C0(:,:), Cr(:,:)
      real(wp), allocatable :: a(:,:), Cown(:,:), C0own(:,:), fill(:,:), hz(:,:)
      real(wp) :: lhs, melt, fi, cint, scale, a2, fo, ff, e_tot, e_own, e_hst, e_bal
      real(wp) :: nan
      integer  :: step
      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      allocate(z_eq(np,nl), h_eq(np,nl), h1(np,nl), w1(np,nl), I1(np,nl), C1(np,nl), F1(np,nl), I0(np,nl), &
               C0(np,nl), Cr(np,nl), a(np,nl), Cown(np,nl), C0own(np,nl), fill(np,nl), hz(np,nl))
      call scene(colat, lon, z_eq, h1, w1, I1, C1, F1, I0, C0, Cr, 0)
      h_eq = 0.0_wp
      where (w1 >= 1.0_wp)                 ! own ice is undefined where the host covers the whole cell
         h_eq = nan;  h1 = nan
      end where
      call new_par(p)
      se%par = p
      call solid_earth_init(se, z_eq, h_eq, w=w1, I_g=0.0_wp*I1, C=Cr)
      call solid_earth_update(se, h1, 500.0_wp, w1, I1, C1, F1)
      call solid_earth_update(se, h1, 500.0_wp, w1, I1, C1, F1)
      a = 1.0_wp - w1
      hz = merge(h_eq, 0.0_wp, a > 0.0_wp)
      call ocean_function(z_eq, hz, C0own, rho_ice, rho_water)
      where (a > 0.0_wp)
         Cown = (se%gg%C - C1)/a
         fill = F1 + a*z_eq*(Cown - C0own)
      elsewhere
         Cown = 0.0_wp
         fill = F1
      end where
      lhs  = sht_grid_surface_integral(sht, se%gg%C*se%gg%rsl)
      melt = -(rho_ice/rho_water)*sht_grid_surface_integral(sht, se%gg%dIg)
      fi   = sht_grid_surface_integral(sht, fill)
      cint = sht_grid_surface_integral(sht, se%gg%C)
      scale = max(abs(melt), abs(fi), 1.0e-30_wp)
      write(*,'(a)') ' (3) mass balance, 0 < w < 1, own ice NaN where w = 1'
      write(*,'(a,es11.3,a,es11.3,a,es11.3)') '     int C rsl =', lhs, '   melt =', melt, '   fill =', fi
      write(*,'(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') '     balance error / scale =', abs(lhs - (melt + fi))/scale, &
           '   |bsl-bsl_ive - fill/C| =', abs(se%bsl - se%bsl_ive - fi/cint), '   |bsl_ive - melt/C| =', &
           abs(se%bsl_ive - melt/cint), '   mass resid ', se%worst_mass_resid
      if (abs(lhs - (melt + fi))/scale > 1.0e-10_wp) then
         write(*,'(a)') '     FAIL: ocean mass balance';  ok = .false.
      end if
      if (abs(se%bsl - se%bsl_ive - fi/cint) > 1.0e-8_wp .or. abs(se%bsl_ive - melt/cint) > 1.0e-8_wp) then
         write(*,'(a)') '     FAIL: bsl and bsl_ive not tied to the melt and the fill';  ok = .false.
      end if
      if (any(se%gg%rsl /= se%gg%rsl) .or. any(se%gg%C /= se%gg%C)) then
         write(*,'(a)') '     FAIL: NaN from the undefined own ice leaked';  ok = .false.
      end if
      if (abs(fi) < 1.0e-3_wp*abs(melt)) then
         write(*,'(a)') '     FAIL: fill negligible, scenario does not test it';  ok = .false.
      end if
      ! the exported fill integral [m³], whole and own part's share, and the balance with it
      a2 = se%earth%r_earth**2
      ff = a2*sht_grid_surface_integral(sht, F1)
      fo = a2*sht_grid_surface_integral(sht, merge(a*z_eq*(Cown - C0own), 0.0_wp, a > 0.0_wp))
      e_tot = abs(se%fill_int - a2*fi)/abs(a2*fi)
      e_own = abs(se%fill_int_own - fo)/abs(fo)
      e_hst = abs((se%fill_int - se%fill_int_own) - ff)/abs(ff)
      e_bal = abs(a2*lhs - (a2*melt + se%fill_int))/(a2*scale)
      write(*,'(a,es11.3,a,es11.3,a)') '     fill_int =', se%fill_int, ' m3   own part =', se%fill_int_own, ' m3'
      write(*,'(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') '     rel error: total ', e_tot, '   own ', e_own, &
           '   host F ', e_hst, '   balance with fill_int ', e_bal
      if (max(e_tot, e_own, e_hst) > 1.0e-9_wp) then
         write(*,'(a)') '     FAIL: exported fill integrals differ from the recomputed ones';  ok = .false.
      end if
      if (e_bal > 1.0e-10_wp) then
         write(*,'(a)') '     FAIL: ocean mass balance does not close with fill_int';  ok = .false.
      end if
      if (abs(fo) < 1.0e-3_wp*abs(ff) .or. abs(ff) < 1.0e-3_wp*abs(fo)) then
         write(*,'(a)') '     FAIL: one part of the fill is negligible, scenario does not test the split';  ok = .false.
      end if
      call solid_earth_finalize(se)
      step = 0
   end subroutine check_balance

   ! --- (4) passthrough vs global lon-lat host -----------------------------------

   subroutine check_passthrough(ok)
      logical, intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth) :: a, b
      type(grid_class)  :: hgrid
      type(remap_ll_gauss) :: map
      real(wp), allocatable :: xl(:), yl(:), tl(:,:), ll(:,:)
      real(wp), allocatable :: z_eq(:,:), h_eq(:,:), h1(:,:), w1(:,:), I1(:,:), C1(:,:), F1(:,:), I0(:,:), C0(:,:), Cr(:,:), Z0(:,:)
      real(wp), allocatable :: zg(:,:), heg(:,:), h1g(:,:), wg(:,:), Ig(:,:), Cg(:,:), Fg(:,:), Crg(:,:), one(:,:)
      real(wp) :: d1, d2, d3
      integer  :: i, j
      allocate(xl(NLON), yl(NLAT), tl(NLON,NLAT), ll(NLON,NLAT))
      do i = 1, NLON;  xl(i) = -180.0_wp + (real(i,wp)-0.5_wp)*(360.0_wp/NLON);  end do
      do j = 1, NLAT;  yl(j) =  -90.0_wp + (real(j,wp)-0.5_wp)*(180.0_wp/NLAT);  end do
      do j = 1, NLAT
         do i = 1, NLON
            tl(i,j) = 90.0_wp - yl(j);  ll(i,j) = modulo(xl(i), 360.0_wp)
         end do
      end do
      call grid_init(hgrid, name="host", mtype="latlon", units="degrees", planet="Spherical Earth", &
                     lon180=.true., x=xl, y=yl)
      allocate(z_eq(NLON,NLAT), h_eq(NLON,NLAT), h1(NLON,NLAT), w1(NLON,NLAT), I1(NLON,NLAT), C1(NLON,NLAT), &
               F1(NLON,NLAT), I0(NLON,NLAT), C0(NLON,NLAT), Cr(NLON,NLAT), Z0(NLON,NLAT))
      call scene(tl, ll, z_eq, h1, w1, I1, C1, F1, I0, C0, Cr, 0)
      h_eq = 0.0_wp;  Z0 = 0.0_wp
      call new_par(p)
      a%par = p;  b%par = p
      ! a: the model remaps from the host grid
      call solid_earth_init(a, z_eq, h_eq, grid=hgrid, w=w1, I_g=Z0, C=Cr)
      call solid_earth_update(a, h1, 500.0_wp, w1, I1, C1, F1)
      ! b: the same fields remapped here, as the model does, handed over on the Gauss grid
      call remap_init(map, sht, xl, yl)
      allocate(zg(np,nl), heg(np,nl), h1g(np,nl), wg(np,nl), Ig(np,nl), Cg(np,nl), Fg(np,nl), Crg(np,nl), one(np,nl))
      call remap_to_gauss(map, sht, z_eq, zg, conserve_mass=.false.)
      call remap_to_gauss(map, sht, h_eq, heg, conserve_mass=.true.)
      call remap_to_gauss(map, sht, h1, h1g, conserve_mass=.true.)
      one = 1.0_wp
      call frac(map, w1, wg, one)                        ! w in [0,1]
      call frac(map, Cr, Crg, wg)
      call remap_to_gauss(map, sht, I1, Ig, conserve_mass=.true.)
      call frac(map, C1, Cg, wg)
      call remap_to_gauss(map, sht, F1, Fg, conserve_mass=.false.)
      call solid_earth_init(b, zg, heg, w=wg, I_g=0.0_wp*Ig, C=Crg)
      call solid_earth_update(b, h1g, 500.0_wp, wg, Ig, Cg, Fg)
      d1 = diff(a%gg%rsl, b%gg%rsl);  d2 = diff(a%gg%C, b%gg%C);  d3 = diff(a%gg%dIg, b%gg%dIg)
      write(*,'(a)') ' (4) passthrough == global lon-lat host (96 x 48)'
      write(*,'(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') '     |drsl| =', d1, ' m   |dC| =', d2, '   |dIg| =', d3, &
           '   max|rsl| =', maxval(abs(a%gg%rsl))
      if (d1 > 1.0e-12_wp .or. d2 > 1.0e-12_wp .or. d3 > 1.0e-12_wp .or. diff(a%gg%C0, b%gg%C0) > 1.0e-12_wp) then
         write(*,'(a)') '     FAIL: remap path differs from passthrough';  ok = .false.
      end if
      if (maxval(abs(a%gg%rsl)) < 1.0_wp) then
         write(*,'(a)') '     FAIL: no signal';  ok = .false.
      end if
      call solid_earth_finalize(a);  call solid_earth_finalize(b)

   end subroutine check_passthrough

   subroutine frac(map, f, g, up)
      !! The model's to_gauss_frac: remap without the mass factor, clip to [0, up].
      type(remap_ll_gauss), intent(in) :: map
      real(wp),             intent(in)  :: f(:,:)
      real(wp),             intent(out) :: g(:,:)
      real(wp),             intent(in)  :: up(:,:)
      call remap_to_gauss(map, sht, f, g, conserve_mass=.false.)
      g = min(max(g, 0.0_wp), up)
   end subroutine frac

   ! --- (5) two hosts summed == one -----------------------------------------------

   subroutine check_two_hosts(ok)
      logical, intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth) :: a, b
      real(wp), allocatable :: z_eq(:,:), h_eq(:,:), h1(:,:), wu(:,:), Iu(:,:), Cu(:,:), Fu(:,:), I0(:,:), C0(:,:), Cru(:,:)
      real(wp), allocatable :: w1(:,:), I1(:,:), C1(:,:), F1(:,:), Cr1(:,:), w2(:,:), I2(:,:), C2(:,:), F2(:,:), Cr2(:,:)
      real(wp), allocatable :: Z0(:,:)
      real(wp) :: d
      allocate(z_eq(np,nl), h_eq(np,nl), h1(np,nl), wu(np,nl), Iu(np,nl), Cu(np,nl), Fu(np,nl), I0(np,nl), C0(np,nl), &
               Cru(np,nl), w1(np,nl), I1(np,nl), C1(np,nl), F1(np,nl), Cr1(np,nl), w2(np,nl), I2(np,nl), C2(np,nl), &
               F2(np,nl), Cr2(np,nl), Z0(np,nl))
      ! host 1 in lon [0,90) (w = 1), host 2 in [180,270) (w = 3/4); the union in one go
      call scene(colat, lon, z_eq, h1, w1, I1, C1, F1, I0, C0, Cr1, 1)
      call scene(colat, lon, z_eq, h1, w2, I2, C2, F2, I0, C0, Cr2, 3)
      call scene(colat, lon, z_eq, h1, wu, Iu, Cu, Fu, I0, C0, Cru, 4)
      ! the union scene also has the [90,180) sector: take it out, so that the union is hosts 1 + 2
      where (lon >= 90.0_wp .and. lon < 180.0_wp)
         wu = 0.0_wp;  Iu = 0.0_wp;  Cu = 0.0_wp;  Fu = 0.0_wp;  Cru = 0.0_wp
      end where
      h_eq = 0.0_wp;  Z0 = 0.0_wp
      call new_par(p)
      a%par = p;  b%par = p
      call solid_earth_init(a, z_eq, h_eq, w=w1 + w2, I_g=Z0, C=Cr1 + Cr2)
      call solid_earth_init(b, z_eq, h_eq, w=wu, I_g=Z0, C=Cru)
      call solid_earth_update(a, h1, 500.0_wp, w1 + w2, I1 + I2, C1 + C2, F1 + F2)
      call solid_earth_update(b, h1, 500.0_wp, wu, Iu, Cu, Fu)
      d = diff(a%gg%rsl, b%gg%rsl)
      write(*,'(a)') ' (5) two hosts summed == one host over both'
      write(*,'(a,es10.2,a,es10.2,a,es10.2)') '     |drsl| =', d, ' m   |dC| =', diff(a%gg%C, b%gg%C), &
           '   max|rsl| =', maxval(abs(a%gg%rsl))
      if (any(a%gg%rsl /= b%gg%rsl) .or. any(a%gg%C /= b%gg%C) .or. a%bsl /= b%bsl) then
         write(*,'(a)') '     FAIL: summed hosts differ from the single host';  ok = .false.
      end if
      call solid_earth_finalize(a);  call solid_earth_finalize(b)
   end subroutine check_two_hosts

   ! --- (6) round-off in the fractions is clipped, also for a Gauss host -----------

   subroutine check_clip(ok)
      logical, intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth) :: a, b
      real(wp), allocatable :: z_eq(:,:), h_eq(:,:), h1(:,:), w1(:,:), I1(:,:), C1(:,:), F1(:,:), I0(:,:), C0(:,:), Cr(:,:)
      real(wp), allocatable :: wq(:,:), Cp(:,:), Z0(:,:)
      allocate(z_eq(np,nl), h_eq(np,nl), h1(np,nl), w1(np,nl), I1(np,nl), C1(np,nl), F1(np,nl), I0(np,nl), &
               C0(np,nl), Cr(np,nl), wq(np,nl), Cp(np,nl), Z0(np,nl))
      call scene(colat, lon, z_eq, h1, w1, I1, C1, F1, I0, C0, Cr, 0)
      h_eq = 0.0_wp;  Z0 = 0.0_wp
      ! b: w overshoots 1 and undershoots 0 by 1e-10, C overshoots w by 1e-10 where C = w (all inside the
      ! tolerance); a: the exact terms
      wq = w1
      where (w1 >= 1.0_wp)
         wq = 1.0_wp + 1.0e-10_wp
      elsewhere (w1 <= 0.0_wp)
         wq = -1.0e-10_wp
      end where
      Cp = C1;  where (w1 >= 1.0_wp .and. C1 >= w1) Cp = C1 + 1.0e-10_wp
      call new_par(p)
      a%par = p;  b%par = p
      call solid_earth_init(a, z_eq, h_eq, w=w1, I_g=Z0, C=Cr)
      call solid_earth_init(b, z_eq, h_eq, w=wq, I_g=Z0, C=Cr)
      call solid_earth_update(a, h1, 500.0_wp, w1, I1, C1, F1)
      call solid_earth_update(b, h1, 500.0_wp, wq, I1, Cp, F1)
      write(*,'(a)') ' (6) round-off beyond [0,1] and [0,w] is clipped (Gauss host)'
      write(*,'(a,es10.2)') '     |drsl| =', diff(a%gg%rsl, b%gg%rsl)
      if (any(a%gg%rsl /= b%gg%rsl) .or. any(a%gg%C /= b%gg%C) .or. any(b%gg%host%last%w > 1.0_wp) &
          .or. any(b%gg%host%last%w < 0.0_wp) .or. any(b%gg%host%last%C > b%gg%host%last%w)) then
         write(*,'(a)') '     FAIL: round-off not clipped';  ok = .false.
      end if
      call solid_earth_finalize(a);  call solid_earth_finalize(b)
   end subroutine check_clip

end program test_host_terms
