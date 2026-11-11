module vilma_remap
   !! Bidirectional conservative/bilinear remapping between a host lon-lat grid and
   !! the model's SHTns Gauss-Legendre grid, built on the fesm-utils `coords` library
   !! (in-package SCRIP-style weights, gen="coords", no CDO).
   !!
   !! Two legs, mirroring CLIMBER-X's VILMA1 coupling:
   !!   host lon-lat --(conservative)--> Gauss   [remap_to_gauss]   mass-bearing fields (ice)
   !!   Gauss --(bilinear)--> host lon-lat        [remap_to_ll]      smooth fields (rsl)
   !!
   !! Build the map pair once with remap_init, then apply per field / time slice. The
   !! coupling layer (vilma_coupling) holds one of these when the host grid differs from
   !! the model Gauss grid, and drives h_ice in / rsl out through it; the standalone
   !! driver (vilma_drive) uses it to remap lon-lat forcing onto the Gauss grid.
   !!
   !! Caching: map_init writes each weight set to a SCRIP(-superset) NetCDF cache under
   !! `fldr` (default "maps") and reloads it on the next run, keyed by the grid names —
   !! which we tag with the grid dimensions so a resolution change invalidates the
   !! cache. The conservative weight build (polygon clipping) is the costly step; the
   !! cache lets a restart skip it. A regional grid's name carries its extent as well; global
   !! grids of the same dimensions but different axes would collide on the name — clear
   !! `fldr` if a global host grid changes at fixed resolution. The key
   !! does not include the planet either: both grids are spherical (planet="Spherical Earth",
   !! as the model), so a cache built before that was set (WGS84 areas) must be cleared.
   !!
   !! Source grid. Longitudes ascending; latitudes ascending or descending (a descending
   !! axis is reversed on the way in and out). The grid may be regional. What that means
   !! depends on the field, and remap_to_gauss takes it from conserve_mass:
   !!   mass (ice thickness, conserve_mass=.true.): the source's mass goes where the source
   !!     has it. A Gauss cell under the source gets the area mean of the source over it; a
   !!     cell a REGIONAL source covers only in part gets that mean times the covered
   !!     fraction of its area (the part's mass over the whole cell), and a cell outside the
   !!     source is zero: ice outside the file is no ice. The field is then scaled so that
   !!     its Gauss-grid integral is the source integral, for a regional source as for a
   !!     global one. (The covered fraction is coords' and good to ~1e-3, so a cell covered
   !!     to 99 % counts as covered. A source that reaches every Gauss cell is global, and
   !!     its cells are not scaled: a polar gap narrower than a cell is not an ice-free cap.)
   !!   anything else (bed, fractions; conserve_mass=.false.): the area mean of the source
   !!     over the cell. There is no value to give a cell outside the source, so the source
   !!     must cover the sphere (remap_covers_sphere) and the call stops if it does not.
   !!
   !! Mass note: the Gauss target is built with coords' grid_init_gaussian, whose cell edges
   !! come from the Gauss weights, so its cell areas are the SHTns quadrature weights
   !! (R^2 dlon w_j) and the conservative map preserves the integral SHTns takes of the
   !! remapped field to round-off. conserve_mass=.true. then rescales by a factor of 1
   !! (to ~1e-16, tests/test_remap) for a global source: the source total in steradians is
   !! sum(f*area)/R^2, with R^2 from the Gauss grid's total area 4*pi*R^2.
   use vilma_precision, only: wp, is_finite
   use vilma_sht,       only: sht_grid, sht_grid_surface_integral
   use coords,       only: grid_class, grid_init, grid_init_gaussian, map_class, map_init, map_field
   implicit none
   private

   public :: remap_ll_gauss
   public :: remap_init, remap_to_gauss, remap_to_ll, remap_covers_sphere
   public :: gauss_axes_rows, rows_copy
   public :: ROWS_NONE, ROWS_NORTH_FIRST, ROWS_SOUTH_FIRST

   !! Row order of a field on the model's Gauss grid (gauss_axes_rows). The model's own layout
   !! (SHTns, se%gg) is north first; a lon-lat host on that grid usually runs south first.
   integer, parameter :: ROWS_NONE        =  0   !! the axes are not those of the Gauss grid
   integer, parameter :: ROWS_NORTH_FIRST =  1   !! latitude descending: the SHTns layout
   integer, parameter :: ROWS_SOUTH_FIRST = -1   !! latitude ascending: rows reversed

   real(wp), parameter :: RAD2DEG = 57.295779513082323_wp
   real(wp), parameter :: FOURPI  = 12.566370614359172_wp
   integer,  parameter :: BIL_NEIGHBORS = 8   !! neighbour pool for the bilinear leg
   !! A Gauss cell a regional source covers beyond this fraction counts as covered: coords'
   !! overlap areas are good to ~1e-3 of a cell, so a covered cell's fraction is not exactly 1.
   real(wp), parameter :: COVER_FULL = 0.99_wp

   type :: remap_ll_gauss
      !! A precomputed map pair between a host lon-lat grid and the model's SHTns
      !! Gauss grid. Build once with remap_init; remap_to_gauss / remap_to_ll apply it.
      type(grid_class) :: ll                  !! host lon-lat grid
      type(grid_class) :: gauss               !! Gauss target (matches sht), lat ascending
      type(map_class)  :: to_gauss            !! ll -> Gauss, conservative
      type(map_class)  :: to_ll               !! Gauss -> ll, bilinear
      integer :: nlon = 0, nlat_ll = 0        !! host dimensions
      integer :: nphi = 0, nlat = 0           !! Gauss (SHTns) dimensions
      logical :: flip_lat = .false.           !! host latitudes descend: rows reversed in and out
      logical :: global   = .false.           !! every Gauss cell has host cells under it
      real(wp), allocatable :: cover(:,:)     !! (nphi,nlat), latitude ascending: fraction of each Gauss
                                              !! cell under a regional host grid (1 for a global one)
   end type remap_ll_gauss

contains

   subroutine remap_init(self, sht, lon_src, lat_src, fldr)
      !! Build (or load from cache) the conservative (ll->Gauss) and bilinear
      !! (Gauss->ll) maps for source cell-centre axes lon_src(:), lat_src(:) [degrees]
      !! against the SHTns Gauss grid. Longitudes must ascend; latitudes may ascend or
      !! descend (strictly), and the grid may be regional (see the module header). The
      !! Gauss grid is built with ASCENDING latitude (south first), which coords expects;
      !! the apply routines flip into/out of the SHTns (north-first) row order. fldr is
      !! the weight-cache directory.
      type(remap_ll_gauss), intent(out) :: self
      type(sht_grid),       intent(in)  :: sht
      real(wp),             intent(in)  :: lon_src(:)   !! source longitudes [deg]
      real(wp),             intent(in)  :: lat_src(:)   !! source latitudes  [deg]
      character(len=*), optional, intent(in) :: fldr    !! weight-cache dir (default "maps")
      character(len=:), allocatable :: cache
      character(len=64) :: ll_name, gauss_name
      real(wp), allocatable :: lat_asc(:), ones(:,:), mean1(:,:)
      logical,  allocatable :: m2(:,:)
      real(wp) :: dx, dy, span
      integer :: j

      self%nlon = size(lon_src);  self%nlat_ll = size(lat_src)
      if (self%nlon < 2 .or. self%nlat_ll < 2) error stop 'vilma_remap: the source grid needs >= 2 points per axis'
      if (.not. (all(is_finite(lon_src)) .and. all(is_finite(lat_src)))) &
         error stop 'vilma_remap: a source longitude or latitude is not finite'
      if (any(lon_src(2:) <= lon_src(:self%nlon-1))) &
         error stop 'vilma_remap: source longitudes must be strictly ascending'
      self%flip_lat = lat_src(2) < lat_src(1)
      lat_asc = lat_src
      if (self%flip_lat) lat_asc = lat_src(self%nlat_ll:1:-1)
      if (any(lat_asc(2:) <= lat_asc(:self%nlat_ll-1))) &
         error stop 'vilma_remap: source latitudes must be strictly ascending or strictly descending'
      if (lat_asc(1) < -90.0_wp .or. lat_asc(self%nlat_ll) > 90.0_wp) &
         error stop 'vilma_remap: source latitudes outside [-90, 90] degrees'
      self%nphi = sht%nphi;       self%nlat    = sht%nlat
      cache = "maps";  if (present(fldr)) cache = trim(fldr)
      call execute_command_line("mkdir -p '"//cache//"'")

      if (sht%mres /= 1) error stop 'vilma_remap: the Gauss target is global (needs mres = 1)'

      ! tag grid names with their dimensions so the weight cache is resolution-specific
      write(ll_name,    '(a,i0,a,i0)') "ll_",    self%nlon, "x", self%nlat_ll
      ! a regional grid is named with its extent too [1/100 degree], so two regions of the
      ! same dimensions do not share a cache; global grids keep the plain name
      dx = lon_src(2) - lon_src(1);  dy = lat_asc(2) - lat_asc(1)
      span = lon_src(self%nlon) - lon_src(1) + dx
      if (abs(span - 360.0_wp) > 0.5_wp*dx .or. lat_asc(1) - dy > -90.0_wp &
          .or. lat_asc(self%nlat_ll) + dy < 90.0_wp) &
         write(ll_name, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') "ll_", self%nlon, "x", self%nlat_ll, "_r", &
              nint(100.0_wp*lon_src(1)), "_", nint(100.0_wp*lon_src(self%nlon)), "_", &
              nint(100.0_wp*lat_asc(1)), "_", nint(100.0_wp*lat_asc(self%nlat_ll))
      write(gauss_name, '(a,i0,a,i0)') "gauss_", self%nphi, "x", self%nlat
      ! a reversed axis gets its own name: its weights are those of the ascending grid
      if (self%flip_lat) ll_name = trim(ll_name)//"_f"
      call grid_init(self%ll, name=trim(ll_name), mtype="latlon", units="degrees", &
                     planet="Spherical Earth", lon180=.true., x=lon_src, y=lat_asc)
      ! Gauss target with edges from the Gauss weights, so its cell areas are the SHTns
      ! quadrature weights (R^2 dlon w_j); lat ascending (south first), as coords expects.
      call grid_init_gaussian(self%gauss, trim(gauss_name), nlon=self%nphi, nlat=self%nlat, &
                              planet="Spherical Earth", lon180=.false.)
      ! the axes must be those of the SHTns grid (the apply routines flip rows only)
      do j = 1, self%nphi
         if (abs(self%gauss%G%x(j) - sht%lon(j)*RAD2DEG) > 1.0e-8_wp) &
            error stop 'vilma_remap: coords Gauss longitudes differ from the SHTns grid'
      end do
      do j = 1, self%nlat
         ! SHTns colat runs north -> south, so lat descends; coords' axis is ascending.
         if (abs(self%gauss%G%y(j) - (90.0_wp - sht%colat(self%nlat - j + 1)*RAD2DEG)) > 1.0e-8_wp) &
            error stop 'vilma_remap: coords Gauss latitudes differ from the SHTns grid'
      end do

      ! in-package weights (gen="coords"), auto-cached under fldr and reloaded if present
      call map_init(self%to_gauss, self%ll, self%gauss, method="con", gen="coords", fldr=cache)
      call map_init(self%to_ll, self%gauss, self%ll, method="bilinear", &
                    max_neighbors=BIL_NEIGHBORS, gen="coords", fldr=cache)

      ! Coverage. Does the source reach every Gauss cell (a field that is not a mass needs
      ! it)? If not, the fraction of each cell it covers, for the mass fields.
      allocate(ones(self%nlon, self%nlat_ll), source=1.0_wp)
      allocate(mean1(self%nphi, self%nlat), m2(self%nphi, self%nlat), self%cover(self%nphi, self%nlat))
      call map_field(self%to_gauss, "cover", ones, mean1, stat="mean", mask2=m2, frac2=self%cover)
      self%global = all(m2)
      if (self%global) then
         self%cover = 1.0_wp
      else
         where (.not. m2)
            self%cover = 0.0_wp
         elsewhere (self%cover > COVER_FULL)
            self%cover = 1.0_wp
         end where
      end if
   end subroutine remap_init

   pure logical function remap_covers_sphere(self) result(g)
      !! Does the source grid reach every Gauss cell? .false. for a regional source, which
      !! can supply a mass field (ice) but not a bed or a fraction (remap_to_gauss).
      type(remap_ll_gauss), intent(in) :: self
      g = self%global
   end function remap_covers_sphere

   subroutine remap_to_gauss(self, sht, f_ll, f_gauss, conserve_mass)
      !! Conservatively remap f_ll(nlon, nlat_ll) onto f_gauss(nphi, nlat) in the SHTns
      !! spatial layout. conserve_mass says what the field is (see the module header):
      !! .true., a mass per area (ice thickness): each Gauss cell gets the source mass that
      !! falls into it, cells outside a regional source are zero, and the SHTns surface
      !! integral equals the source area-integral. .false. (default), anything else (bed,
      !! fractions): the area mean over the cell; stops unless the source covers the sphere.
      type(remap_ll_gauss), intent(in)  :: self
      type(sht_grid),       intent(in)  :: sht
      real(wp),             intent(in)  :: f_ll(:,:)    !! (nlon, nlat_ll)
      real(wp),             intent(out) :: f_gauss(:,:) !! (nphi, nlat)
      logical, optional,    intent(in)  :: conserve_mass
      real(wp), allocatable :: vt(:,:), src(:,:)
      logical,  allocatable :: m2(:,:)
      real(wp) :: src_int_sr, gauss_int
      integer  :: j
      logical  :: do_mass

      if (size(f_ll,1) /= self%nlon .or. size(f_ll,2) /= self%nlat_ll) &
         error stop 'vilma_remap: source field shape /= map source grid'
      if (size(f_gauss,1) /= self%nphi .or. size(f_gauss,2) /= self%nlat) &
         error stop 'vilma_remap: target field shape /= SHTns grid'
      do_mass = .false.;  if (present(conserve_mass)) do_mass = conserve_mass
      if (.not. do_mass .and. .not. self%global) &
         error stop 'vilma_remap: a field that is not a mass (bed, fraction) needs a source grid '// &
                    'that covers the sphere; this one is regional'

      allocate(src(self%nlon, self%nlat_ll))              ! the source with latitude ascending
      if (self%flip_lat) then
         do j = 1, self%nlat_ll
            src(:, j) = f_ll(:, self%nlat_ll - j + 1)
         end do
      else
         src = f_ll
      end if

      allocate(vt(self%nphi, self%nlat), m2(self%nphi, self%nlat))
      call map_field(self%to_gauss, "f", src, vt, stat="mean", mask2=m2)
      if (do_mass) then
         ! a mass: the covered part's mean over the whole cell, none outside the source
         ! (cover is 1 throughout for a global source, so nothing changes there)
         where (m2)
            vt = vt*self%cover
         elsewhere
            vt = 0.0_wp
         end where
      end if

      do j = 1, self%nlat                                 ! ascending-lat -> SHTns north-first
         f_gauss(:, j) = vt(:, self%nlat - j + 1)
      end do

      if (do_mass) then
         ! the source total as a steradian integral, sum(f*area)/R^2 with 4*pi*R^2 the
         ! Gauss grid's total area (the source's own area may be regional). The factor is
         ! 1 to round-off for a global source and to the accuracy of the covered
         ! fractions (~1e-3 of the edge cells' mass) for a regional one.
         src_int_sr = sum(src * self%ll%area) / sum(self%gauss%area) * FOURPI
         gauss_int  = sht_grid_surface_integral(sht, f_gauss)
         if (abs(gauss_int) > tiny(1.0_wp)) f_gauss = f_gauss * (src_int_sr/gauss_int)
      end if
   end subroutine remap_to_gauss

   subroutine remap_to_ll(self, f_gauss, f_ll)
      !! Bilinearly remap f_gauss(nphi, nlat) in the SHTns spatial layout onto the host
      !! lon-lat grid f_ll(nlon, nlat_ll). Used for smooth fields (relative sea level)
      !! returned to the host; the host reconstructs z_bed on its own high-resolution
      !! bed as z_bed_eq - f_ll, so no mass rescale is applied here. Rows in the host's
      !! own latitude order.
      type(remap_ll_gauss), intent(in)  :: self
      real(wp),             intent(in)  :: f_gauss(:,:) !! (nphi, nlat), SHTns north-first
      real(wp),             intent(out) :: f_ll(:,:)    !! (nlon, nlat_ll)
      real(wp), allocatable :: g_asc(:,:), row(:)
      logical,  allocatable :: m2(:,:)
      integer :: j

      if (size(f_gauss,1) /= self%nphi .or. size(f_gauss,2) /= self%nlat) &
         error stop 'vilma_remap: source field shape /= SHTns grid'
      if (size(f_ll,1) /= self%nlon .or. size(f_ll,2) /= self%nlat_ll) &
         error stop 'vilma_remap: target field shape /= map target grid'

      allocate(g_asc(self%nphi, self%nlat), m2(self%nlon, self%nlat_ll))
      do j = 1, self%nlat                                 ! SHTns north-first -> ascending-lat
         g_asc(:, j) = f_gauss(:, self%nlat - j + 1)
      end do

      call map_field(self%to_ll, "f", g_asc, f_ll, method="bilinear", mask2=m2)
      where (.not. m2) f_ll = 0.0_wp                      ! uncovered host cells (full-sphere Gauss src: none)
      if (self%flip_lat) then                             ! back to the host's descending latitudes
         allocate(row(self%nlon))
         do j = 1, self%nlat_ll/2
            row = f_ll(:, j);  f_ll(:, j) = f_ll(:, self%nlat_ll - j + 1);  f_ll(:, self%nlat_ll - j + 1) = row
         end do
      end if
   end subroutine remap_to_ll

   integer function gauss_axes_rows(sht, lon, lat) result(rows)
      !! Are lon(:), lat(:) [degrees] the axes of the Gauss grid sht (same sizes, every value
      !! within 1e-6 degree; longitudes as SHTns has them, from 0), and in which row order:
      !! ROWS_NORTH_FIRST (latitude descending, the SHTns layout), ROWS_SOUTH_FIRST (ascending)
      !! or ROWS_NONE (another grid).
      type(sht_grid), intent(in) :: sht
      real(wp),       intent(in) :: lon(:), lat(:)
      real(wp), parameter :: TOL = 1.0e-6_wp           ! degrees
      real(wp) :: lat_n(sht%nlat)                      ! the SHTns latitudes, north first
      rows = ROWS_NONE
      if (size(lon) /= sht%nphi .or. size(lat) /= sht%nlat) return
      if (any(abs(lon - sht%lon*RAD2DEG) > TOL)) return
      lat_n = 90.0_wp - sht%colat*RAD2DEG
      if (all(abs(lat - lat_n) <= TOL)) then
         rows = ROWS_NORTH_FIRST
      else if (all(abs(lat - lat_n(sht%nlat:1:-1)) <= TOL)) then
         rows = ROWS_SOUTH_FIRST
      end if
   end function gauss_axes_rows

   subroutine rows_copy(a, b, rows)
      !! b = a between the SHTns layout and a field on the Gauss grid in row order rows: a copy
      !! (ROWS_NORTH_FIRST) or the rows reversed (ROWS_SOUTH_FIRST); the same either way round.
      real(wp), intent(in)  :: a(:,:)
      real(wp), intent(out) :: b(:,:)
      integer,  intent(in)  :: rows
      integer :: j, n
      if (any(shape(a) /= shape(b))) error stop 'vilma_remap: rows_copy between fields of different shape'
      n = size(a,2)
      select case (rows)
      case (ROWS_NORTH_FIRST)
         b = a
      case (ROWS_SOUTH_FIRST)
         do j = 1, n
            b(:, j) = a(:, n - j + 1)
         end do
      case default
         error stop 'vilma_remap: rows_copy needs ROWS_NORTH_FIRST or ROWS_SOUTH_FIRST'
      end select
   end subroutine rows_copy

end module vilma_remap
