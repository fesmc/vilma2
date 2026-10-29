   subroutine write_synth_visc_fesmdata(file)
      !! A small synthetic 3-D viscosity file in the FesmData layout (log10_visc on
      !! lon -180:180, lat, depth [m] from the surface down, structure), standing in
      !! for the FesmData Bagge2021 product so that `make check` needs no data file.
      !! Structure 6 (the default visc_3d_structure) carries a smooth mantle profile
      !! plus a lateral anomaly of ~0.7 dex in the upper mantle, which makes elements
      !! genuinely 3-D (6 memory channels with l_toroidal); the other structures hold
      !! the mean profile. Included by the tests that need a lateral viscosity.
      use ncio, only: nc_create, nc_write_dim, nc_write
      character(len=*), intent(in) :: file
      integer,  parameter :: NLON = 90, NLAT = 61, NDEP = 30, NSTR = 6
      real(wp), parameter :: DEG = 3.14159265358979323846_wp/180.0_wp
      real(wp), allocatable :: lon(:), lat(:), dep(:), eta(:,:,:,:)
      integer,  allocatable :: istr(:)
      real(wp) :: prof, lat_anom
      integer  :: i, j, k, s
      allocate(lon(NLON), lat(NLAT), dep(NDEP), istr(NSTR), eta(NLON,NLAT,NDEP,NSTR))
      do i = 1, NLON;  lon(i) = -180.0_wp + 360.0_wp*real(i-1,wp)/real(NLON,wp);  end do
      do j = 1, NLAT;  lat(j) = -90.0_wp + 180.0_wp*real(j-1,wp)/real(NLAT-1,wp);  end do
      do k = 1, NDEP;  dep(k) = 2891.0e3_wp*real(k-1,wp)/real(NDEP-1,wp);          end do
      istr = [(s, s = 1, NSTR)]
      do s = 1, NSTR
         do k = 1, NDEP
            prof = 20.5_wp + 1.0_wp*dep(k)/2891.0e3_wp                  ! 20.5 ... 21.5
            do j = 1, NLAT
               do i = 1, NLON
                  lat_anom = 0.0_wp
                  if (s == NSTR) lat_anom = 0.7_wp*cos(2.0_wp*lon(i)*DEG)*cos(lat(j)*DEG)**2 &
                                                  *exp(-dep(k)/600.0e3_wp)
                  eta(i,j,k,s) = prof + lat_anom
               end do
            end do
         end do
      end do
      call nc_create(file, overwrite=.true.)
      call nc_write_dim(file, "lon",       x=lon,  units="degrees_east")
      call nc_write_dim(file, "lat",       x=lat,  units="degrees_north")
      call nc_write_dim(file, "depth",     x=dep,  units="m")
      call nc_write_dim(file, "structure", x=istr, units="1")
      call nc_write(file, "log10_visc", eta, dim1="lon", dim2="lat", dim3="depth", dim4="structure")
      deallocate(lon, lat, dep, istr, eta)
   end subroutine write_synth_visc_fesmdata
