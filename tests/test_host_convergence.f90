program test_host_convergence
   !! Convergence of the host-resolved terms against the Gauss-grid own path. A global lon-lat
   !! host (0.5 deg: 720 x 360) resolves an analytic geometry -- a land cap, a marine basin with
   !! grounded ice over its deeper part, a shelf that crosses sea level, a steep continental
   !! slope, all modulated in longitude -- and hands its grounded columns I_g, ocean fraction C
   !! and fill F (w = 1) to the model; the own path judges flotation and the coastline on the
   !! Gauss cells. Gauss resolutions lmax = 8, 16, 32, 48 (nlat = 2 lmax + 2).
   !!
   !! Level A, prescribed sea level (analytic rsl): the water-budget integrals
   !!     Mi = ∫ΔI_g,   Mw = ∫C·rsl,   Mf = ∫fill
   !! and the L1 distance of C from the remapped host coastline, against the 0.5 deg truth, with
   !! the host's cell areas as the map weights have them (coords', up to 0.9% from the exact
   !! sphere at this resolution; the conservative remap preserves sum(f*area) in them). The
   !! own-path error falls with lmax, roughly as 1/lmax (the coast is a line, its cell-wise
   !! error is first order) but noisily, since a coarse cell flips between ocean and
   !! grounded; the host path stays at its floor, set by the host resolution and the product
   !! of means.
   !!
   !! Level B, the full elastic sea-level solve, one step (rho_ice, rho_water as given, in
   !! the host, the flotation and the solve alike): host-only (w = 1) against own, at each
   !! lmax. The host's coastline is that of the lmax-48 own solution. |bsl_host - bsl_own|
   !! and |dphi_host - dphi_own| fall with lmax, as the own path approaches the host-resolved
   !! geometry.
   !!
   !! Run for the two density pairs (910, 1028) and the GIA benchmark pair (931, 1000).
   !! The weights of the host -> Gauss maps are cached under maps/.
   use vilma_precision, only: wp
   use vilma_constants, only: pi, sec_per_year
   use vilma_params,    only: vilma_param_class
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_surface_integral
   use vilma_sle,       only: ocean_function
   use vilma_remap,     only: remap_ll_gauss, remap_init, remap_to_gauss
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   use coords,          only: grid_class, grid_init
   implicit none

   integer,  parameter :: NLON = 720, NLAT = 360, NL = 4
   integer,  parameter :: LS(NL) = [8, 16, 32, 48]
   real(wp), parameter :: DEG = pi/180.0_wp
   real(wp), parameter :: RHO(2,2) = reshape([910.0_wp, 1028.0_wp, 931.0_wp, 1000.0_wp], [2,2])
   type(grid_class) :: hgrid
   real(wp), allocatable :: xl(:), yl(:), th(:,:), la(:,:), area(:,:)
   real(wp), allocatable :: bed_f(:,:), ice_f(:,:), rsl_f(:,:), zero_f(:,:), one_f(:,:)
   real(wp), allocatable :: rsl_ref(:,:)
   logical :: ok

   ok = .true.
   allocate(xl(NLON), yl(NLAT), th(NLON,NLAT), la(NLON,NLAT), area(NLON,NLAT))
   allocate(bed_f(NLON,NLAT), ice_f(NLON,NLAT), rsl_f(NLON,NLAT), zero_f(NLON,NLAT), one_f(NLON,NLAT))
   call fine_grid()
   zero_f = 0.0_wp;  one_f = 1.0_wp
   call geometry(th, la, bed_f, ice_f, rsl_f)

   call level_a(RHO(:,1), ok)
   call level_a(RHO(:,2), ok)
   call level_b(RHO(:,1), ok)
   call level_b(RHO(:,2), ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: host-resolved terms converge to the host geometry; the own path to them'
      call radial_fe_finalize()
   else
      write(*,'(a)') ' FAIL: convergence of the host-resolved terms'
      call radial_fe_finalize()
      error stop 1
   end if

contains

   subroutine fine_grid()
      integer  :: i, j
      real(wp) :: dlon, dlat, s0, s1
      dlon = 360.0_wp/NLON;  dlat = 180.0_wp/NLAT
      do i = 1, NLON;  xl(i) = -180.0_wp + (real(i,wp)-0.5_wp)*dlon;  end do
      do j = 1, NLAT;  yl(j) =  -90.0_wp + (real(j,wp)-0.5_wp)*dlat;  end do
      call grid_init(hgrid, name="host05", mtype="latlon", units="degrees", lon180=.true., x=xl, y=yl)
      do j = 1, NLAT
         s0 = sin((yl(j) - 0.5_wp*dlat)*DEG);  s1 = sin((yl(j) + 0.5_wp*dlat)*DEG)
         do i = 1, NLON
            th(i,j) = 90.0_wp - yl(j);  la(i,j) = xl(i)*DEG
            area(i,j) = dlon*DEG*(s1 - s0)             ! steradians
         end do
      end do
      write(*,'(a,es10.2)') ' host cell areas, coords vs exact sphere: max relative difference =', &
           maxval(abs(hgrid%area/sum(hgrid%area)*4.0_wp*pi/area - 1.0_wp))
      ! the host's areas are those of the map weights (coords): the conservative remap preserves
      ! sum(f*area) in them, so the truth is formed with them (normalised to the sphere)
      area = hgrid%area/sum(hgrid%area)*4.0_wp*pi
   end subroutine fine_grid

   subroutine geometry(t, l, bed, ice, rsl)
      !! Analytic bed, grounded-ice thickness at t1 and prescribed sea-level change at colatitude t [deg]
      !! and longitude l [rad].
      real(wp), intent(in)  :: t(:,:), l(:,:)
      real(wp), intent(out) :: bed(:,:), ice(:,:), rsl(:,:)
      integer  :: i, j
      real(wp) :: x
      do j = 1, size(t,2)
         do i = 1, size(t,1)
            x = t(i,j)
            if (x < 25.0_wp) then
               bed(i,j) = 400.0_wp
            else if (x < 45.0_wp) then
               bed(i,j) = 400.0_wp - 900.0_wp*(x - 25.0_wp)/20.0_wp       ! marine basin, to -500 m
            else if (x < 65.0_wp) then
               bed(i,j) = -500.0_wp + 800.0_wp*(x - 45.0_wp)/20.0_wp      ! shelf: crosses sea level
            else if (x < 80.0_wp) then
               bed(i,j) = 300.0_wp
            else if (x < 85.0_wp) then
               bed(i,j) = 300.0_wp - 4300.0_wp*(x - 80.0_wp)/5.0_wp       ! continental slope
            else
               bed(i,j) = -4000.0_wp
            end if
            bed(i,j) = bed(i,j) + 150.0_wp*sin(3.0_wp*l(i,j))*exp(-((x - 55.0_wp)/15.0_wp)**2)
            ice(i,j) = 0.0_wp
            if (x < 42.0_wp) ice(i,j) = 2500.0_wp*(1.0_wp - (x/42.0_wp)**2)*(1.0_wp + 0.3_wp*cos(2.0_wp*l(i,j)))
            rsl(i,j) = -80.0_wp + 30.0_wp*cos(x*DEG) + 10.0_wp*sin(2.0_wp*l(i,j))*sin(x*DEG)
         end do
      end do
   end subroutine geometry

   subroutine host_terms(rho, rsl, bed, ice, I_g, C, F, Cref)
      !! The host's terms on the fine grid at sea-level change rsl: binary flotation at its resolution.
      real(wp), intent(in)  :: rho(2)
      real(wp), intent(in)  :: rsl(:,:), bed(:,:), ice(:,:)
      real(wp), intent(out) :: I_g(:,:), C(:,:), F(:,:), Cref(:,:)
      call ocean_function(bed - rsl, ice, C, rho(1), rho(2))
      call ocean_function(bed, 0.0_wp*ice, Cref, rho(1), rho(2))
      I_g = ice*(1.0_wp - C)
      F   = bed*(C - Cref)
   end subroutine host_terms

   real(wp) function fine_int(f)
      real(wp), intent(in) :: f(:,:)
      fine_int = sum(f*area)
   end function fine_int

   ! --- Level A ---------------------------------------------------------------------

   subroutine level_a(rho, ok)
      real(wp), intent(in)    :: rho(2)
      logical,  intent(inout) :: ok
      real(wp), allocatable :: I_f(:,:), C_f(:,:), F_f(:,:), Cr_f(:,:)
      real(wp), allocatable :: topo(:,:), H(:,:), rs(:,:), Cg_t(:,:), Ih(:,:), Ch(:,:), Fh(:,:)
      real(wp), allocatable :: Co(:,:), C0o(:,:), Io(:,:), Fo(:,:), tg(:,:), lg(:,:), dum(:,:)
      type(sht_grid), target :: sht
      type(remap_ll_gauss) :: map
      real(wp) :: Mi_t, Mw_t, Mf_t, sc, eo(3,NL), eh(3,NL), l1o(NL), l1h(NL), sa
      integer  :: k, L, np, nl_, i, j
      allocate(I_f(NLON,NLAT), C_f(NLON,NLAT), F_f(NLON,NLAT), Cr_f(NLON,NLAT))
      call host_terms(rho, rsl_f, bed_f, ice_f, I_f, C_f, F_f, Cr_f)
      Mi_t = fine_int(I_f);  Mw_t = fine_int(C_f*rsl_f);  Mf_t = fine_int(F_f)
      sc = abs(Mi_t)
      write(*,'(a,f6.0,a,f6.0,a)') ' Level A (prescribed rsl), rho_i, rho_w =', rho(1), ',', rho(2), ' kg m-3'
      write(*,'(a,3es12.4,a)') '   truth (0.5 deg): Mi, Mw, Mf =', Mi_t, Mw_t, Mf_t, '  [m sr]'
      write(*,'(a)') '   lmax   own: |dMi| |dMw| |dMf| / |Mi|              host: |dMi| |dMw| |dMf| / |Mi|        L1(C) own   host'
      do k = 1, NL
         L = LS(k)
         call sht_grid_init(sht, L, nlat=2*L+2, nphi=4*L)
         np = sht%nphi;  nl_ = sht%nlat
         allocate(topo(np,nl_), H(np,nl_), rs(np,nl_), Cg_t(np,nl_), Ih(np,nl_), Ch(np,nl_), Fh(np,nl_), &
                  Co(np,nl_), C0o(np,nl_), Io(np,nl_), Fo(np,nl_), tg(np,nl_), lg(np,nl_), dum(np,nl_))
         call remap_init(map, sht, xl, yl)
         ! own path: cell means of the bed and the ice, analytic rsl at the cell centres
         call remap_to_gauss(map, sht, bed_f, topo, conserve_mass=.false.)
         call remap_to_gauss(map, sht, ice_f, H,    conserve_mass=.true.)
         do j = 1, nl_
            do i = 1, np
               tg(i,j) = sht%colat(j)/DEG;  lg(i,j) = sht%lon(i)
               if (lg(i,j) > pi) lg(i,j) = lg(i,j) - 2.0_wp*pi         ! the fine grid's longitudes: -180..180
            end do
         end do
         call geometry(tg, lg, dum, dum, rs)
         call ocean_function(topo - rs, H, Co, rho(1), rho(2))
         call ocean_function(topo, 0.0_wp*H, C0o, rho(1), rho(2))
         Io = H*(1.0_wp - Co);  Fo = topo*(Co - C0o)
         ! host path: the fine terms remapped
         call remap_to_gauss(map, sht, C_f, Ch, conserve_mass=.false.)
         call remap_to_gauss(map, sht, I_f, Ih, conserve_mass=.true.)
         call remap_to_gauss(map, sht, F_f, Fh, conserve_mass=.false.)
         Cg_t = Ch
         eo(1,k) = abs(sht_grid_surface_integral(sht, Io) - Mi_t)/sc
         eo(2,k) = abs(sht_grid_surface_integral(sht, Co*rs) - Mw_t)/sc
         eo(3,k) = abs(sht_grid_surface_integral(sht, Fo) - Mf_t)/sc
         eh(1,k) = abs(sht_grid_surface_integral(sht, Ih) - Mi_t)/sc
         eh(2,k) = abs(sht_grid_surface_integral(sht, Ch*rs) - Mw_t)/sc
         eh(3,k) = abs(sht_grid_surface_integral(sht, Fh) - Mf_t)/sc
         sa = sht_grid_surface_integral(sht, 0.0_wp*Co + 1.0_wp)
         l1o(k) = sht_grid_surface_integral(sht, abs(Co - Cg_t))/sa
         l1h(k) = 0.0_wp
         write(*,'(i6,3es11.2,4x,3es11.2,2es11.2)') L, eo(:,k), eh(:,k), l1o(k), l1h(k)
         deallocate(topo, H, rs, Cg_t, Ih, Ch, Fh, Co, C0o, Io, Fo, tg, lg, dum)
         call sht_grid_destroy(sht)
      end do
      ! the own path converges, noisily (a coarse cell flips between ocean and grounded as lmax
      ! changes): the mean error of the two finest grids is below that of the two coarsest, the
      ! L1 distance of its coastline falls monotonically; the host path sits at a floor far below
      ! (Mi conserved by the remap to round-off, Mw and Mf to the product-of-means error)
      do i = 1, 3
         if (.not. (0.5_wp*(eo(i,3) + eo(i,4)) < 0.5_wp*(eo(i,1) + eo(i,2)))) then
            write(*,'(a,i0,a)') '   FAIL: own-path error of term ', i, ' does not fall with lmax';  ok = .false.
         end if
      end do
      if (.not. (l1o(2) < l1o(1) .and. l1o(3) < l1o(2) .and. l1o(4) < l1o(3) .and. l1o(4) < 0.5_wp*l1o(1))) then
         write(*,'(a)') '   FAIL: the own coastline does not converge to the host one';  ok = .false.
      end if
      if (.not. all(eh(1,:) < 1.0e-10_wp) .or. .not. all(eh(2,:) < 5.0e-4_wp) .or. .not. all(eh(3,:) < 5.0e-4_wp)) then
         write(*,'(a)') '   FAIL: the host-path floor is too high';  ok = .false.
      end if
      if (.not. (maxval(eh(2,:)) < minval(eo(2,:)) .and. maxval(eh(3,:)) < minval(eo(3,:)))) then
         write(*,'(a)') '   FAIL: the host path is not below the own path at every lmax';  ok = .false.
      end if
   end subroutine level_a

   ! --- Level B ---------------------------------------------------------------------

   subroutine level_b(rho, ok)
      real(wp), intent(in)    :: rho(2)
      logical,  intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth) :: own, host
      real(wp), allocatable :: I_f(:,:), C_f(:,:), F_f(:,:), Cr_f(:,:), Ir_f(:,:)
      real(wp) :: db(NL), dp(NL), t0
      integer  :: k, L
      allocate(I_f(NLON,NLAT), C_f(NLON,NLAT), F_f(NLON,NLAT), Cr_f(NLON,NLAT), Ir_f(NLON,NLAT))
      write(*,'(a,f6.0,a,f6.0,a)') ' Level B (elastic SLE, one step), rho_i, rho_w =', rho(1), ',', rho(2), ' kg m-3'
      ! the host's coastline: that of the own solution at the finest lmax
      call new_par(p, LS(NL), rho)
      own%par = p
      call solid_earth_init(own, bed_f, zero_f, grid=hgrid)
      call solid_earth_update(own, ice_f, 1.0e3_wp)
      rsl_ref = own%rsl
      call host_terms(rho, rsl_ref, bed_f, ice_f, I_f, C_f, F_f, Cr_f)
      Ir_f = 0.0_wp
      call solid_earth_finalize(own)
      write(*,'(a)') '   lmax   bsl own     bsl host    |dbsl|      dphi own    dphi host   |ddphi|'
      do k = 1, NL
         L = LS(k)
         call new_par(p, L, rho)
         own%par = p;  host%par = p
         call solid_earth_init(own, bed_f, zero_f, grid=hgrid)
         call solid_earth_update(own, ice_f, 1.0e3_wp)
         call solid_earth_init(host, bed_f, zero_f, grid=hgrid, w=one_f, I_g=Ir_f, C=Cr_f)
         call solid_earth_update(host, zero_f, 1.0e3_wp, one_f, I_f, C_f, F_f)
         db(k) = abs(host%bsl - own%bsl);  dp(k) = abs(host%dphi - own%dphi)
         write(*,'(i6,6f12.5)') L, own%bsl, host%bsl, db(k), own%dphi, host%dphi, dp(k)
         if (host%worst_mass_resid > 1.0e-10_wp) then
            write(*,'(a)') '     FAIL: host solve does not conserve ocean mass';  ok = .false.
         end if
         call solid_earth_finalize(own);  call solid_earth_finalize(host)
      end do
      t0 = 0.0_wp
      if (.not. (db(NL) < 0.5_wp*db(1) .and. dp(NL) < 0.5_wp*dp(1))) then
         write(*,'(a)') '   FAIL: |bsl_host - bsl_own| and |dphi_host - dphi_own| do not fall with lmax';  ok = .false.
      end if
   end subroutine level_b

   subroutine new_par(p, lmax, rho)
      type(vilma_param_class), intent(out) :: p
      integer,  intent(in) :: lmax
      real(wp), intent(in) :: rho(2)
      p%lmax = lmax;  p%nlat = 2*lmax + 2;  p%nphi = 4*lmax
      p%earth_response = "elastic"
      p%rotation = .false.
      p%rho_ice = rho(1);  p%rho_water = rho(2)
   end subroutine new_par

end program test_host_convergence
