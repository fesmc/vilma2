program test_host_orient
   !! The row order of a host on the model's Gauss grid. VILMA2 works in the SHTns layout
   !! (rows north first); a host on the same grid built by coords runs south first. The
   !! passthrough (no remap) must be an exact row-reversing copy both ways, and every
   !! host-grid output in the host's order.
   !!
   !! (1) Round trip, a hemispherically asymmetric field: with a south-first Gauss host the
   !!     reference fields and the host terms reach se%gg rows reversed, and rsl, dz_bed,
   !!     dz_ss, C come back rows reversed (z_bed = z_bed_eq - rsl); with a north-first
   !!     Gauss host (and with no grid) both are plain copies.
   !! (2) A northern ice cap on a hemispherically asymmetric 3-D viscosity (stiff north,
   !!     weak south): the same load through the south-first host, through the north-first
   !!     host, through no grid and through the remap from a lon-lat host gives the same
   !!     response, while the mirrored load (a southern cap) differs, so a mirrored layout
   !!     would show. The remap run fixes the geography: the others take its Gauss-grid
   !!     inputs, so all four solve the same problem.
   use vilma_precision, only: wp
   use vilma_constants, only: pi, rad2deg
   use vilma_params,    only: vilma_param_class
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_remap,     only: ROWS_NORTH_FIRST, ROWS_SOUTH_FIRST
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   use coords,          only: grid_class, grid_init, grid_init_gaussian
   implicit none

   integer,  parameter :: LMAX = 12, NL = 32, NP = 4*LMAX          ! NL: a size SHTns keeps
   integer,  parameter :: NLON = 90, NLATH = 45                 ! the lon-lat host of (2), 4 degrees
   integer,  parameter :: NSTEP = 2
   real(wp), parameter :: DT_YR = 1.0e3_wp
   character(len=*), parameter :: VISC = "libvilma/include/test_host_orient_visc.nc"

   type(grid_class) :: g_south, g_north, g_ll
   real(wp), allocatable :: lat_s(:), lon_g(:), lon_ll(:), lat_ll(:)
   logical :: ok

   ok = .true.
   call grid_init_gaussian(g_south, "orient_gauss_s", nlon=NP, nlat=NL, lon180=.false.)
   allocate(lon_g(NP), lat_s(NL))
   lon_g = g_south%G%x;  lat_s = g_south%G%y
   call grid_init(g_north, name="orient_gauss_n", mtype="latlon", units="degrees", planet="Spherical Earth", &
                  lon180=.false., x=lon_g, y=lat_s(NL:1:-1))
   allocate(lon_ll(NLON), lat_ll(NLATH))
   call axes_ll()
   call grid_init(g_ll, name="orient_ll", mtype="latlon", units="degrees", planet="Spherical Earth", &
                  lon180=.true., x=lon_ll, y=lat_ll)

   call round_trip()
   call write_visc_asym(VISC)
   call response_orient()

   write(*,'(a)') ''
   call radial_fe_finalize()
   if (ok) then
      write(*,'(a)') ' PASS: host row order: exact reversing passthrough, same response as the remap path'
   else
      write(*,'(a)') ' FAIL: host row order'
      error stop 1
   end if

contains

   subroutine axes_ll()
      integer :: i
      do i = 1, NLON;   lon_ll(i) = -180.0_wp + (real(i,wp) - 0.5_wp)*(360.0_wp/NLON);  end do
      do i = 1, NLATH;  lat_ll(i) =  -90.0_wp + (real(i,wp) - 0.5_wp)*(180.0_wp/NLATH); end do
   end subroutine axes_ll

   real(wp) function cap(lon, lat, lat0) result(h)
      !! A 2 km cap of 20 degree radius (cosine taper) centred at (60 E, lat0).
      real(wp), intent(in) :: lon, lat, lat0
      real(wp) :: a, b, dl, d
      a = lat*pi/180.0_wp;  b = lat0*pi/180.0_wp;  dl = (lon - 60.0_wp)*pi/180.0_wp
      d = acos(max(-1.0_wp, min(1.0_wp, sin(a)*sin(b) + cos(a)*cos(b)*cos(dl))))*rad2deg
      h = 0.0_wp
      if (d < 20.0_wp) h = 2000.0_wp*0.5_wp*(1.0_wp + cos(pi*d/20.0_wp))
   end function cap

   real(wp) function bed(lon, lat) result(z)
      !! Land north of 35 N and in a southern box, deep ocean elsewhere: not symmetric.
      real(wp), intent(in) :: lon, lat
      z = -4000.0_wp
      if (lat > 35.0_wp) z = 400.0_wp
      if (lat < -20.0_wp .and. lat > -50.0_wp .and. lon > 100.0_wp .and. lon < 160.0_wp) z = 200.0_wp
   end function bed

   function rev(a) result(b)
      real(wp), intent(in) :: a(:,:)
      real(wp) :: b(size(a,1), size(a,2))
      b = a(:, size(a,2):1:-1)
   end function rev

   subroutine expect(cond, what)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: what
      if (.not. cond) then
         write(*,'(a)') '   FAIL: '//what;  ok = .false.
      end if
   end subroutine expect

   subroutine round_trip()
      !! (1) an elastic Earth (quick), a northern cap, host terms in an asymmetric region.
      type(vilma_param_class) :: p
      type(solid_earth) :: s, n
      real(wp), allocatable :: z(:,:), h0(:,:), h(:,:), w(:,:), Ig(:,:), C(:,:), F(:,:)
      integer :: i, j
      allocate(z(NP,NL), h0(NP,NL), h(NP,NL), w(NP,NL), Ig(NP,NL), C(NP,NL), F(NP,NL))
      do j = 1, NL
         do i = 1, NP
            z(i,j) = bed(lon_g(i), lat_s(j));  h(i,j) = cap(lon_g(i), lat_s(j), 65.0_wp)
            w(i,j) = merge(1.0_wp, 0.0_wp, lat_s(j) > 50.0_wp .and. lon_g(i) < 120.0_wp)
         end do
      end do
      h0 = 0.0_wp;  Ig = w*h;  C = 0.0_wp;  F = 0.0_wp
      p%lmax = LMAX;  p%nlat = NL;  p%nphi = NP
      p%earth_response = "elastic";  p%rotation = .false.

      ! the south-first host (south-first host fields: the arrays as laid out by lat_s)
      s%par = p
      call solid_earth_init(s, z, h0, grid=g_south, w=w, I_g=h0, C=C)
      call expect(.not. s%remap .and. s%host_rows == ROWS_SOUTH_FIRST, 'south-first Gauss host not recognised')
      call expect(all(s%gg%z_bed_eq == rev(z)) .and. all(s%gg%host%ref%w == rev(w)), &
                  'south-first host: z_bed_eq or the reference w not rows reversed on the Gauss grid')
      call expect(all(s%C == rev(s%gg%C)) .and. all(s%z_bed == z), 'south-first host: reference C or z_bed')
      call solid_earth_update(s, h0, 100.0_wp, w=w, I_g=Ig, C=C, F=F)
      call expect(all(s%gg%host%last%w == rev(w)) .and. all(s%gg%host%last%I_g == rev(Ig)), &
                  'south-first host: host terms not rows reversed on the Gauss grid')
      call expect(all(s%rsl == rev(s%gg%rsl)) .and. all(s%dz_bed == rev(s%gg%dz_bed)) .and. &
                  all(s%dz_ss == rev(s%gg%dz_ss)) .and. all(s%C == rev(s%gg%C)), &
                  'south-first host: rsl, dz_bed, dz_ss or C not rows reversed on the host grid')
      call expect(all(s%z_bed == z - s%rsl), 'south-first host: z_bed /= z_bed_eq - rsl')
      call expect(minval(s%dz_bed(:, NL/2+1:)) < 10.0_wp*minval(s%dz_bed(:, :NL/2)) .and. &
                  minval(s%dz_bed(:, NL/2+1:)) < -1.0_wp, 'south-first host: the subsidence is not under the northern cap')
      write(*,'(a,2es11.3)') ' (1) south-first host: min dz_bed north, south of the equator [m]', &
           minval(s%dz_bed(:, NL/2+1:)), minval(s%dz_bed(:, :NL/2))

      ! the north-first host: the same problem, the arrays reversed, a plain copy
      n%par = p
      call solid_earth_init(n, rev(z), h0, grid=g_north, w=rev(w), I_g=h0, C=C)
      call expect(.not. n%remap .and. n%host_rows == ROWS_NORTH_FIRST, 'north-first Gauss host not recognised')
      call expect(all(n%gg%z_bed_eq == rev(z)), 'north-first host: z_bed_eq not copied')
      call solid_earth_update(n, h0, 100.0_wp, w=rev(w), I_g=rev(Ig), C=C, F=F)
      call expect(all(n%rsl == n%gg%rsl) .and. all(n%dz_bed == n%gg%dz_bed) .and. all(n%C == n%gg%C), &
                  'north-first host: outputs not copied')
      call expect(all(n%gg%dz_bed == s%gg%dz_bed) .and. all(n%gg%rsl == s%gg%rsl), &
                  'north-first and south-first hosts solve different problems')
      call solid_earth_finalize(s);  call solid_earth_finalize(n)
   end subroutine round_trip

   subroutine response_orient()
      !! (2) the viscoelastic response on the asymmetric 3-D viscosity, four ways and mirrored.
      type(vilma_param_class) :: p
      type(solid_earth) :: r, s, n, g, m
      real(wp), allocatable :: zl(:,:), hl(:,:), zg(:,:), hg(:,:), zero_l(:,:), zero_g(:,:)
      real(wp) :: big, d_s, d_n, d_g, d_m
      integer  :: i, j, k

      p%lmax = LMAX;  p%nlat = NL;  p%nphi = NP
      p%earth_response = "ve";  p%rotation = .false.
      p%l_visc_3d = .true.;  p%visc_3d_file = VISC;  p%visc_3d_structure = 1;  p%l_toroidal = .true.

      ! the remap run: the geography on a lon-lat host (latitude ascending)
      allocate(zl(NLON,NLATH), hl(NLON,NLATH), zero_l(NLON,NLATH), zero_g(NP,NL))
      do j = 1, NLATH
         do i = 1, NLON
            zl(i,j) = bed(modulo(lon_ll(i), 360.0_wp), lat_ll(j));  hl(i,j) = cap(lon_ll(i), lat_ll(j), 65.0_wp)
         end do
      end do
      zero_l = 0.0_wp;  zero_g = 0.0_wp
      r%par = p
      call solid_earth_init(r, zl, zero_l, grid=g_ll)
      call expect(r%remap, 'lon-lat host: no remap built')
      zg = r%gg%z_bed_eq                                   ! the remapped inputs, SHTns layout
      allocate(hg(NP,NL))

      s%par = p;  call solid_earth_init(s, rev(zg), zero_g, grid=g_south)
      n%par = p;  call solid_earth_init(n, zg, zero_g, grid=g_north)
      g%par = p;  call solid_earth_init(g, zg, zero_g)
      m%par = p;  call solid_earth_init(m, rev(zg), zero_g)  ! mirrored: no grid, the south-first arrays
      do k = 1, NSTEP
         call solid_earth_update(r, hl, DT_YR)
         if (k == 1) hg = r%gg%h_ice
         call solid_earth_update(s, rev(hg), DT_YR)
         call solid_earth_update(n, hg, DT_YR)
         call solid_earth_update(g, hg, DT_YR)
         call solid_earth_update(m, rev(hg), DT_YR)
      end do

      big = maxval(abs(r%gg%dz_bed))
      d_s = maxval(abs(rev(s%dz_bed) - r%gg%dz_bed))
      d_n = maxval(abs(n%dz_bed - r%gg%dz_bed))
      d_g = maxval(abs(g%dz_bed - r%gg%dz_bed))
      d_m = maxval(abs(rev(m%gg%dz_bed) - r%gg%dz_bed))
      write(*,'(a,es11.3,a)') ' (2) max|dz_bed| (remap run)  =', big, ' m;  max |difference to it| [m]:'
      write(*,'(a,3es11.3,a,es11.3)') '     south-first, north-first, no grid =', d_s, d_n, d_g, &
           ';  mirrored load =', d_m
      call expect(d_s <= 1.0e-10_wp*big .and. d_n <= 1.0e-10_wp*big .and. d_g <= 1.0e-10_wp*big, &
                  'the Gauss hosts do not reproduce the remap run')
      call expect(maxval(abs(rev(s%rsl) - r%gg%rsl)) <= 1.0e-10_wp*maxval(abs(r%gg%rsl)), &
                  'south-first host: rsl differs from the remap run')
      call expect(d_m > 0.02_wp*big, 'the mirrored load gives the same response: the viscosity is not asymmetric enough')
      call solid_earth_finalize(r);  call solid_earth_finalize(s);  call solid_earth_finalize(n)
      call solid_earth_finalize(g);  call solid_earth_finalize(m)
   end subroutine response_orient

   subroutine write_visc_asym(file)
      !! A 3-D viscosity in the FesmData layout (log10_visc on lon, lat, depth, structure; one
      !! structure): a mantle profile plus, in the upper 700 km, +-0.8 dex with the sign of the
      !! latitude (stiff north, weak south).
      use ncio, only: nc_create, nc_write_dim, nc_write
      character(len=*), intent(in) :: file
      integer,  parameter :: NLO = 72, NLA = 37, NDEP = 30
      real(wp), allocatable :: lon(:), lat(:), dep(:), eta(:,:,:,:)
      integer  :: i, j, k
      allocate(lon(NLO), lat(NLA), dep(NDEP), eta(NLO,NLA,NDEP,1))
      do i = 1, NLO;   lon(i) = -180.0_wp + 360.0_wp*real(i-1,wp)/real(NLO,wp);  end do
      do j = 1, NLA;   lat(j) = -90.0_wp + 180.0_wp*real(j-1,wp)/real(NLA-1,wp);  end do
      do k = 1, NDEP;  dep(k) = 2891.0e3_wp*real(k-1,wp)/real(NDEP-1,wp);          end do
      do k = 1, NDEP
         do j = 1, NLA
            do i = 1, NLO
               eta(i,j,k,1) = 20.7_wp + 1.0_wp*dep(k)/2891.0e3_wp &
                            + merge(0.8_wp*tanh(lat(j)/15.0_wp), 0.0_wp, dep(k) < 700.0e3_wp)
            end do
         end do
      end do
      call nc_create(file, overwrite=.true.)
      call nc_write_dim(file, "lon",       x=lon,  units="degrees_east")
      call nc_write_dim(file, "lat",       x=lat,  units="degrees_north")
      call nc_write_dim(file, "depth",     x=dep,  units="m")
      call nc_write_dim(file, "structure", x=[1],  units="1")
      call nc_write(file, "log10_visc", eta, dim1="lon", dim2="lat", dim3="depth", dim4="structure")
   end subroutine write_visc_asym

end program test_host_orient
