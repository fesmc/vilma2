program test_visc_fesmdata
   !! The FesmData Bagge2021_visc.nc product against the isostasy_data file that
   !! CLIMBER-X reads today (input/vilma2/visc3d_Bagge2021.nc = isostasy_data
   !! earth_structure/viscosity/bagge2021.nc): structure 6 (v_1.0_s16) of the
   !! product must give the same log10(eta) on the model grid as that file.
   !!
   !! Both files are read by vilma_read_visc_3d, onto the lmax = 170 Gauss grid at the
   !! radii of the reference file plus a few off-layer radii. Needs the two data
   !! files, which are not in the repository, so it is NOT in `make check`:
   !!
   !!   make test_visc_fesmdata
   !!   bin/test_visc_fesmdata.x <Bagge2021_visc.nc> <bagge2021.nc> [structure=6]
   use vilma_precision,       only: wp
   use vilma_earth_structure, only: vilma_read_visc_3d
   use vilma_sht,             only: sht_grid, sht_grid_init
   use ncio,                  only: nc_read, nc_size
   implicit none

   integer,  parameter :: LMAX = 170
   real(wp), parameter :: TOL = 1.0e-6_wp          !! [dex]
   character(len=512) :: f_fesm, f_ref, arg
   integer            :: structure, nr_ref, nr, k, i
   real(wp), allocatable :: r_ref(:), r_node(:)
   real(wp), allocatable :: v_fesm(:,:), v_ref(:,:)
   type(sht_grid) :: sht
   real(wp) :: maxdiff, vmin, vmax

   if (command_argument_count() < 2) then
      write(*,'(a)') ' usage: test_visc_fesmdata.x <Bagge2021_visc.nc> <bagge2021.nc> [structure]'
      error stop 'test_visc_fesmdata: files not given'
   end if
   call get_command_argument(1, f_fesm)
   call get_command_argument(2, f_ref)
   structure = 6
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg);  read(arg,*) structure
   end if

   call sht_grid_init(sht, LMAX, nlat=2*LMAX+2, nphi=4*LMAX, mmax=LMAX)

   ! radii: the reference file's layers, and the mid-point of every 8th pair
   nr_ref = nc_size(f_ref, 'r')
   allocate(r_ref(nr_ref))
   call nc_read(f_ref, 'r', r_ref)
   nr = nr_ref + (nr_ref-1)/8
   allocate(r_node(nr))
   r_node(1:nr_ref) = r_ref
   i = nr_ref
   do k = 1, nr_ref-1, 8
      i = i + 1
      r_node(i) = 0.5_wp*(r_ref(k) + r_ref(k+1))
   end do

   call vilma_read_visc_3d(f_fesm, sht, r_node, v_fesm, structure=structure)
   call vilma_read_visc_3d(f_ref,  sht, r_node, v_ref)

   maxdiff = maxval(abs(v_fesm - v_ref))
   vmin = minval(v_fesm);  vmax = maxval(v_fesm)
   write(*,'(a,i0,a,f6.2,a,f6.2)') ' FesmData structure ', structure, ': log10(eta) min=', vmin, ' max=', vmax
   write(*,'(a,es10.2,a)') ' max |FesmData - reference| = ', maxdiff, ' dex'
   if (vmin /= vmin .or. vmax /= vmax) error stop 'test_visc_fesmdata FAILED: NaN in the FesmData field'
   if (maxdiff > TOL) error stop 'test_visc_fesmdata FAILED: fields differ'
   write(*,'(a)') ' test_visc_fesmdata PASSED'
end program test_visc_fesmdata
