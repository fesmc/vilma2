program test_remap
   !! Conservative lon-lat -> Gauss remap (vilma_remap): constant preservation, latitude
   !! orientation (the ascending<->SHTns-row flip), zonal-field accuracy, and the
   !! optional global mass-rescale (SHTns surface integral == source area integral).
   !! Then the source-grid cases: latitudes descending (the same result as ascending), and
   !! a regional source (the Northern Hemisphere): an ice field keeps its mass and its
   !! thickness and is zero outside the source, and the grid reports that it does not
   !! cover the sphere (a bed or a fraction cannot come from it).
   use vilma_precision, only: wp
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_surface_integral, sht_grid_destroy
   use vilma_remap,     only: remap_ll_gauss, remap_init, remap_to_gauss, remap_to_ll, remap_covers_sphere
   implicit none

   integer,  parameter :: LMAX = 32
   integer,  parameter :: NLON = 72, NLAT_S = 36       ! 5-degree source grid
   real(wp), parameter :: DEG = 3.141592653589793_wp/180.0_wp
   real(wp), parameter :: RAD2DEG = 57.295779513082323_wp

   type(sht_grid)      :: sht
   type(remap_ll_gauss)  :: rmap
   real(wp), allocatable :: lon_s(:), lat_s(:), fsrc(:,:), fdst(:,:)
   real(wp) :: latg, emax, e, src_mean, gauss_mean, isrc, idst
   integer  :: i, j
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)

   allocate(lon_s(NLON), lat_s(NLAT_S))
   do i = 1, NLON;   lon_s(i) = -180.0_wp + (real(i,wp)-0.5_wp)*(360.0_wp/NLON);  end do
   do j = 1, NLAT_S; lat_s(j) =  -90.0_wp + (real(j,wp)-0.5_wp)*(180.0_wp/NLAT_S); end do
   call remap_init(rmap, sht, lon_s, lat_s)

   allocate(fsrc(NLON,NLAT_S), fdst(sht%nphi,sht%nlat))

   ! (1) constant field preserved exactly
   fsrc = 7.0_wp
   call remap_to_gauss(rmap, sht, fsrc, fdst)
   emax = maxval(abs(fdst - 7.0_wp))
   write(*,'(a,es10.2)') '   constant-field max error      = ', emax
   if (emax > 1.0e-10_wp) then; write(*,*) 'FAIL: constant not preserved'; ok = .false.; end if

   ! (2) zonal field f = sin(lat): accuracy + orientation (north row positive)
   do j = 1, NLAT_S
      do i = 1, NLON;  fsrc(i,j) = sin(lat_s(j)*DEG);  end do
   end do
   call remap_to_gauss(rmap, sht, fsrc, fdst)
   emax = 0.0_wp
   do j = 1, sht%nlat
      latg = 90.0_wp - sht%colat(j)*RAD2DEG
      e = maxval(abs(fdst(:,j) - sin(latg*DEG)))
      emax = max(emax, e)
   end do
   write(*,'(a,es10.2)') '   zonal sin(lat) max error      = ', emax
   if (emax > 2.0e-2_wp) then; write(*,*) 'FAIL: zonal field beyond resolution tol'; ok = .false.; end if
   ! orientation: SHTns row 1 is the northern pole-most cell (lat>0 -> sin>0)
   if (fdst(1,1) <= 0.0_wp .or. fdst(1,sht%nlat) >= 0.0_wp) then
      write(*,*) 'FAIL: latitude orientation flipped'; ok = .false.
   end if
   write(*,'(a,f8.3,a,f8.3)') '   north row=', fdst(1,1), '   south row=', fdst(1,sht%nlat)

   ! (3) mass-rescale: a positive bump, conserve_mass=.true. -> SHTns integral matches
   !     the source area-integral (mean) to machine precision.
   do j = 1, NLAT_S
      do i = 1, NLON
         fsrc(i,j) = 1.0_wp + 0.5_wp*cos(lat_s(j)*DEG)*cos(lon_s(i)*DEG)
      end do
   end do
   call remap_to_gauss(rmap, sht, fsrc, fdst, conserve_mass=.true.)
   ! source area-weighted mean and the model's gauss-weighted mean (integral/4pi)
   isrc = sum(fsrc*rmap%ll%area)/sum(rmap%ll%area)
   idst = sht_grid_surface_integral(sht, fdst)/(16.0_wp*atan(1.0_wp))
   write(*,'(a,f12.8,a,f12.8,a,es10.2)') '   src mean=', isrc, '  gauss mean=', idst, &
        '  rel=', abs(isrc-idst)/abs(isrc)
   if (abs(isrc-idst) > 1.0e-10_wp*abs(isrc)) then
      write(*,*) 'FAIL: mass not conserved with conserve_mass'; ok = .false.
   end if

   ! (3b) with Gauss edges from the Gauss weights the cell areas are the SHTns quadrature
   !      weights, so the plain remap already preserves the integral: the conserve_mass
   !      factor is 1 to round-off.
   call remap_to_gauss(rmap, sht, fsrc, fdst)
   idst = sht_grid_surface_integral(sht, fdst)/(16.0_wp*atan(1.0_wp))
   write(*,'(a,es10.2)') '   unscaled remap integral vs source, rel = ', abs(isrc-idst)/abs(isrc)
   if (abs(isrc-idst) > 1.0e-10_wp*abs(isrc)) then
      write(*,*) 'FAIL: plain remap does not preserve the integral'; ok = .false.
   end if

   ! (4) reverse leg: bilinear Gauss -> host lon-lat reproduces a smooth zonal field
   block
      real(wp), allocatable :: gsrc(:,:), fll(:,:)
      real(wp) :: latg2
      allocate(gsrc(sht%nphi,sht%nlat), fll(NLON,NLAT_S))
      do j = 1, sht%nlat
         latg2 = 90.0_wp - sht%colat(j)*RAD2DEG
         gsrc(:,j) = sin(latg2*DEG)            ! smooth field on the Gauss grid (SHTns rows)
      end do
      call remap_to_ll(rmap, gsrc, fll)
      emax = 0.0_wp
      do j = 1, NLAT_S
         e = maxval(abs(fll(:,j) - sin(lat_s(j)*DEG)))
         emax = max(emax, e)
      end do
      write(*,'(a,es10.2)') '   reverse bilinear sin(lat) err = ', emax
      if (emax > 5.0e-2_wp) then; write(*,*) 'FAIL: reverse bilinear beyond tol'; ok = .false.; end if
   end block

   ! (5) source latitudes descending: the same Gauss field as the ascending source, and the
   !     reverse leg returns rows in the source's own (descending) order
   block
      type(remap_ll_gauss)  :: rdesc
      real(wp), allocatable :: fa(:,:), fd(:,:), la(:,:), ld(:,:)
      allocate(fa(sht%nphi,sht%nlat), fd(sht%nphi,sht%nlat), la(NLON,NLAT_S), ld(NLON,NLAT_S))
      call remap_init(rdesc, sht, lon_s, lat_s(NLAT_S:1:-1))
      call remap_to_gauss(rmap,  sht, fsrc, fa, conserve_mass=.true.)
      call remap_to_gauss(rdesc, sht, fsrc(:, NLAT_S:1:-1), fd, conserve_mass=.true.)
      call remap_to_ll(rmap,  fa, la)
      call remap_to_ll(rdesc, fa, ld)
      write(*,'(a,es10.2,a,es10.2)') '   descending source: max|to_gauss diff| = ', maxval(abs(fd - fa)), &
           '   max|to_ll diff| = ', maxval(abs(ld(:, NLAT_S:1:-1) - la))
      if (maxval(abs(fd - fa)) > 1.0e-12_wp .or. maxval(abs(ld(:, NLAT_S:1:-1) - la)) > 1.0e-12_wp) then
         write(*,*) 'FAIL: a source with descending latitudes is not remapped as the ascending one'; ok = .false.
      end if
      if (.not. (remap_covers_sphere(rmap) .and. remap_covers_sphere(rdesc))) then
         write(*,*) 'FAIL: a global source reported as regional'; ok = .false.
      end if
   end block

   ! (6) regional source, the Northern Hemisphere on the same 5-degree cells: an ice cap of
   !     1000 m north of 60N. Its mass on the Gauss grid is the source's (not 4 pi R^2 / A_src
   !     = 2 times it), the thickness is not inflated, and south of the source there is none.
   block
      integer, parameter :: NLAT_N = NLAT_S/2
      type(remap_ll_gauss)  :: rnh
      real(wp), allocatable :: latn(:), ice(:,:), g(:,:)
      real(wp) :: m_src, m_gauss, r2, south
      allocate(latn(NLAT_N), ice(NLON,NLAT_N), g(sht%nphi,sht%nlat))
      do j = 1, NLAT_N;  latn(j) = (real(j,wp)-0.5_wp)*(180.0_wp/NLAT_S);  end do
      call remap_init(rnh, sht, lon_s, latn)
      do j = 1, NLAT_N
         ice(:,j) = merge(1000.0_wp, 0.0_wp, latn(j) > 60.0_wp)
      end do
      call remap_to_gauss(rnh, sht, ice, g, conserve_mass=.true.)
      r2      = sum(rnh%gauss%area)/(16.0_wp*atan(1.0_wp))          ! R^2 of the coords sphere
      m_src   = sum(ice*rnh%ll%area)                                ! source volume [m^3]
      m_gauss = sht_grid_surface_integral(sht, g)*r2
      south = 0.0_wp
      do j = 1, sht%nlat
         if (90.0_wp - sht%colat(j)*RAD2DEG < -3.0_wp) south = max(south, maxval(abs(g(:,j))))
      end do
      write(*,'(a,es10.2,a,f9.3,a,es9.2,a,l1)') '   regional (NH) ice: volume rel. error = ', &
           abs(m_gauss - m_src)/m_src, '   max thickness = ', maxval(g), ' m   max south of it = ', south, &
           '   covers sphere = ', remap_covers_sphere(rnh)
      if (abs(m_gauss - m_src) > 1.0e-10_wp*m_src) then
         write(*,*) 'FAIL: a regional ice field does not keep its mass'; ok = .false.
      end if
      if (maxval(g) > 1000.0_wp*(1.0_wp + 1.0e-3_wp) .or. maxval(g) < 999.0_wp) then
         write(*,*) 'FAIL: a regional ice field changes its thickness'; ok = .false.
      end if
      if (south /= 0.0_wp) then
         write(*,*) 'FAIL: ice outside the regional source'; ok = .false.
      end if
      if (remap_covers_sphere(rnh)) then
         write(*,*) 'FAIL: a hemisphere reported as covering the sphere'; ok = .false.
      end if
      ! ice up to the edge of the source: the Gauss row the edge cuts carries its share
      ice = 1000.0_wp
      call remap_to_gauss(rnh, sht, ice, g, conserve_mass=.true.)
      m_src   = sum(ice*rnh%ll%area)
      m_gauss = sht_grid_surface_integral(sht, g)*r2
      write(*,'(a,es10.2,a,f9.3,a)') '   regional (NH) ice to the edge: volume rel. error = ', &
           abs(m_gauss - m_src)/m_src, '   max thickness = ', maxval(g), ' m'
      if (abs(m_gauss - m_src) > 1.0e-10_wp*m_src .or. maxval(g) > 1000.0_wp*(1.0_wp + 1.0e-3_wp) &
          .or. minval(g(:, 1:sht%nlat/2 - 2)) < 1000.0_wp*(1.0_wp - 1.0e-3_wp)) then
         write(*,*) 'FAIL: regional ice reaching the edge of the source is not conserved'; ok = .false.
      end if
   end block

   call sht_grid_destroy(sht)
   if (ok) then
      write(*,'(a)') ' PASS: vilma_remap conservative lon-lat -> Gauss + bilinear Gauss -> lon-lat'
   else
      write(*,'(a)') ' FAIL: vilma_remap'
      error stop 1
   end if
end program test_remap
