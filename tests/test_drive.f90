program test_drive
   !! Standalone driver (vilma_drive): generate a synthetic reference state and an
   !! ice-thickness forcing on the model Gauss grid, write a sparse &vilma config
   !! overlaid on the shipped defaults, run vilma_run, and check the output:
   !!   (1) one output slice per forcing slice,
   !!   (2) the bed subsides under the growing ice cap (z_bed < z_bed_eq),
   !!   (3) the ocean draws down (rsl < 0 in the far field).
   !! Then the legacy Gauss-grid reader (remap_input=.false.) with the same files stored
   !! with latitude ascending (rows reversed): read with i_eq = 0 (file_ref) and i_eq = 2
   !! (z_bed_eq_file, h_ice_eq_file), the output must equal that of the north-first files.
   !! Then the restart clock: a run stopped at the middle slice and resumed from its final
   !! restart continues there (one slice left, the same end state as the uninterrupted
   !! run), and forcing_window refuses every window that would replay or skip forcing.
   use vilma_precision,       only: wp
   use vilma_constants,       only: rad2deg, pi
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_drive,           only: vilma_run, forcing_window, forcing_time_ok
   use vilma_radial_fe,       only: radial_fe_finalize
   use ncio,               only: nc_create, nc_write_dim, nc_write, nc_read, nc_size
   implicit none

   integer, parameter :: LMAX = 8, NLAT = 16, NPHI = 32, NT = 3
   character(len=*), parameter :: REF   = "libvilma/include/test_drive_ref.nc"
   character(len=*), parameter :: FORCE = "libvilma/include/test_drive_force.nc"
   character(len=*), parameter :: OUT   = "libvilma/include/test_drive_out.nc"
   character(len=*), parameter :: HOR   = "libvilma/include/test_drive_hor.nc"
   character(len=*), parameter :: CFG   = "libvilma/include/test_drive.nml"
   character(len=*), parameter :: REFA  = "libvilma/include/test_drive_ref_asc.nc"
   character(len=*), parameter :: FORCEA = "libvilma/include/test_drive_force_asc.nc"
   character(len=*), parameter :: OUTA0 = "libvilma/include/test_drive_out_asc0.nc"
   character(len=*), parameter :: OUTA2 = "libvilma/include/test_drive_out_asc2.nc"
   character(len=*), parameter :: CFGA  = "libvilma/include/test_drive_asc.nml"
   character(len=*), parameter :: DEFS  = "input/vilma_defaults.nml"
   character(len=*), parameter :: OUTL1 = "libvilma/include/test_drive_out_leg1.nc"
   character(len=*), parameter :: OUTL2 = "libvilma/include/test_drive_out_leg2.nc"
   character(len=*), parameter :: RST   = "libvilma/include/test_drive_rst.nc"

   type(sht_grid), target :: sht
   real(wp), allocatable  :: lon_deg(:), lat_deg(:), z_bed_eq(:,:), h_ice_eq(:,:)
   real(wp), allocatable  :: h_ice(:,:,:), tyr(:), zb(:,:), rsl(:,:)
   real(wp), allocatable  :: ue(:,:), un(:,:), uet(:,:), unt(:,:)
   real(wp), allocatable  :: a3(:,:,:), b3(:,:,:)
   real(wp) :: thd, sub
   integer  :: i, j, k, jice, jocean, nout
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=NLAT, nphi=NPHI)

   allocate(lon_deg(sht%nphi), lat_deg(sht%nlat))
   lon_deg = sht%lon*rad2deg
   lat_deg = 90.0_wp - sht%colat*rad2deg

   ! --- synthetic reference state (polar land cap over deep ocean) -------------
   allocate(z_bed_eq(sht%nphi,sht%nlat), h_ice_eq(sht%nphi,sht%nlat))
   do j = 1, sht%nlat
      thd = sht%colat(j)*rad2deg
      do i = 1, sht%nphi
         z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
      end do
   end do
   h_ice_eq = 0.0_wp

   call nc_create(REF, overwrite=.true.)
   call nc_write_dim(REF, "lon", x=lon_deg, units="degrees_east")
   call nc_write_dim(REF, "lat", x=lat_deg, units="degrees_north")
   call nc_write(REF, "z_bed_eq",  z_bed_eq,  dim1="lon", dim2="lat")
   call nc_write(REF, "h_ice_eq", h_ice_eq, dim1="lon", dim2="lat")

   ! --- forcing: a grounded cap (colat<30) growing 0 -> 1 -> 2 km --------------
   allocate(h_ice(sht%nphi,sht%nlat,NT), tyr(NT))
   tyr = [0.0_wp, 1000.0_wp, 2000.0_wp]            ! years
   do k = 1, NT
      do j = 1, sht%nlat
         thd = sht%colat(j)*rad2deg
         do i = 1, sht%nphi
            h_ice(i,j,k) = merge(1000.0_wp*real(k-1,wp), 0.0_wp, thd < 30.0_wp)
         end do
      end do
   end do

   call nc_create(FORCE, overwrite=.true.)
   call nc_write_dim(FORCE, "lon",  x=lon_deg, units="degrees_east")
   call nc_write_dim(FORCE, "lat",  x=lat_deg, units="degrees_north")
   call nc_write_dim(FORCE, "time", x=tyr,     units="years", unlimited=.true.)
   call nc_write(FORCE, "h_ice", h_ice, dim1="lon", dim2="lat", dim3="time")

   ! --- run config: sparse &vilma (overlaid on the physics defaults) + a COMPLETE
   ! --- &ctl group (the run config owns &ctl in full; it has no defaults file).
   call write_cfg(CFG, FORCE, REF, OUT, HOR, 0)
   call vilma_run(CFG, defaults_file=DEFS)

   ! --- checks -----------------------------------------------------------------
   jice   = nearest_row(15.0_wp)
   jocean = nearest_row(70.0_wp)

   nout = nc_size(OUT, "time")
   write(*,'(a,i0)') '   output slices: ', nout
   if (nout /= NT) then
      write(*,'(a)') '   FAIL: expected one output slice per forcing slice'; ok = .false.
   end if

   allocate(zb(sht%nphi,sht%nlat), rsl(sht%nphi,sht%nlat))
   call nc_read(OUT, "z_bed", zb,  start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   call nc_read(OUT, "rsl",   rsl, start=[1,1,nout], count=[sht%nphi,sht%nlat,1])

   sub = z_bed_eq(1,jice) - zb(1,jice)
   write(*,'(a,f10.3,a)') '   bed subsidence under ice =', sub, ' m'
   write(*,'(a,f10.3,a)') '   far-field ocean rsl      =', rsl(1,jocean), ' m'
   if (sub < 10.0_wp) then
      write(*,'(a)') '   FAIL: no appreciable subsidence under the growing ice cap'; ok = .false.
   end if
   if (rsl(1,jocean) >= 0.0_wp) then
      write(*,'(a)') '   FAIL: building land ice should draw the ocean down'; ok = .false.
   end if

   ! --- the horizontal-displacement file (&ctl file_hor) --------------------------
   ! An axisymmetric cap on a radially symmetric Earth: the displacement is purely
   ! meridional (u_east ≡ 0) and there is no toroidal field at all, so the
   ! toroidal components must be exactly zero, not merely small.
   write(*,'(a,i0)') '   horizontal slices: ', nc_size(HOR, "time")
   if (nc_size(HOR, "time") /= NT) then
      write(*,'(a)') '   FAIL: expected one horizontal slice per forcing slice'; ok = .false.
   end if
   allocate(ue(sht%nphi,sht%nlat), un(sht%nphi,sht%nlat), uet(sht%nphi,sht%nlat), unt(sht%nphi,sht%nlat))
   call nc_read(HOR, "u_east",      ue,  start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   call nc_read(HOR, "u_north",     un,  start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   call nc_read(HOR, "u_east_tor",  uet, start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   call nc_read(HOR, "u_north_tor", unt, start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   write(*,'(a,f10.4,a,es10.2,a)') '   max|u_north| =', maxval(abs(un)), ' m,  max|u_east| =', &
        maxval(abs(ue)), ' m'
   write(*,'(a,2es10.2)') '   max|toroidal| (east, north) =', maxval(abs(uet)), maxval(abs(unt))
   if (maxval(abs(un)) < 0.01_wp) then
      write(*,'(a)') '   FAIL: no horizontal displacement under the growing cap'; ok = .false.
   end if
   if (maxval(abs(ue)) > 1.0e-12_wp*maxval(abs(un)) .or. maxval(abs(uet)) /= 0.0_wp &
       .or. maxval(abs(unt)) /= 0.0_wp) then
      write(*,'(a)') '   FAIL: an axisymmetric load on a radial Earth moved east or toroidally'; ok = .false.
   end if

   ! --- the legacy reader with latitude ascending: the same fields, rows reversed -------
   call nc_create(REFA, overwrite=.true.)
   call nc_write_dim(REFA, "lon", x=lon_deg, units="degrees_east")
   call nc_write_dim(REFA, "lat", x=lat_deg(sht%nlat:1:-1), units="degrees_north")
   call nc_write(REFA, "z_bed_eq", z_bed_eq(:, sht%nlat:1:-1), dim1="lon", dim2="lat")
   call nc_write(REFA, "h_ice_eq", h_ice_eq(:, sht%nlat:1:-1), dim1="lon", dim2="lat")
   call nc_create(FORCEA, overwrite=.true.)
   call nc_write_dim(FORCEA, "lon",  x=lon_deg, units="degrees_east")
   call nc_write_dim(FORCEA, "lat",  x=lat_deg(sht%nlat:1:-1), units="degrees_north")
   call nc_write_dim(FORCEA, "time", x=tyr,     units="years", unlimited=.true.)
   call nc_write(FORCEA, "h_ice", h_ice(:, sht%nlat:1:-1, :), dim1="lon", dim2="lat", dim3="time")
   allocate(a3(sht%nphi,sht%nlat,NT), b3(sht%nphi,sht%nlat,NT))
   call write_cfg(CFGA, FORCEA, REFA, OUTA0, "", 0)
   call vilma_run(CFGA, defaults_file=DEFS)
   call same_output(OUTA0, 'ascending files, i_eq = 0')
   call write_cfg(CFGA, FORCEA, REFA, OUTA2, "", 2)
   call vilma_run(CFGA, defaults_file=DEFS)
   call same_output(OUTA2, 'ascending files, i_eq = 2')

   call restart_clock()
   call time_axis()

   call sht_grid_destroy(sht)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: vilma_run reads ref+forcing, marches the model,'
      write(*,'(a)') '       and writes a sensible subsiding/draw-down output'
   else
      write(*,'(a)') ' FAIL: standalone driver did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine restart_clock()
      !! A restart continues where it left off. Leg 1 stops at the middle slice (t = 1000 yr);
      !! leg 2 resumes its final restart with time_init unset and must march the one
      !! remaining interval, 1000 -> 2000 yr, to the end state of the uninterrupted run.
      !! (Before the clock was checked, leg 2 replayed the forcing from t = 0 on the
      !! restart's memory.) Then the window rule itself, case by case.
      real(wp), allocatable :: t2(:), z2(:,:), zf(:,:)
      character(len=400) :: msg
      integer :: k0, k1
      logical :: w
      call write_cfg(CFGA, FORCE, REF, OUTL1, "", 2, time_end=1000.0_wp)
      call vilma_run(CFGA, defaults_file=DEFS)
      call execute_command_line("cp libvilma/include/final/vilma_restart.nc "//RST)
      call write_cfg(CFGA, FORCE, REF, OUTL2, "", 2, restart=RST)
      call vilma_run(CFGA, defaults_file=DEFS)
      allocate(t2(nc_size(OUTL2, "time")), z2(sht%nphi,sht%nlat), zf(sht%nphi,sht%nlat))
      call nc_read(OUTL2, "time", t2)
      write(*,'(a,i0,a,2f10.1)') '   resumed run: ', size(t2), ' slices, t first/last =', t2(1), t2(size(t2))
      if (size(t2) /= 2 .or. t2(1) /= 1000.0_wp .or. t2(size(t2)) /= 2000.0_wp) then
         write(*,'(a)') '   FAIL: the resumed run does not continue at the restart time'; ok = .false.
      else
         call nc_read(OUTL2, "z_bed", z2, start=[1,1,2],  count=[sht%nphi,sht%nlat,1])
         call nc_read(OUT,   "z_bed", zf, start=[1,1,NT], count=[sht%nphi,sht%nlat,1])
         write(*,'(a,es10.2,a)') '   resumed run: max|z_bed - uninterrupted run| at the end =', maxval(abs(z2 - zf)), ' m'
         if (maxval(abs(z2 - zf)) > 1.0e-5_wp) then
            write(*,'(a)') '   FAIL: the resumed run does not reach the uninterrupted end state'; ok = .false.
         end if
      end if

      ! the window rule (tyr = 0, 1000, 2000 yr; -1e30 = time_init unset)
      w = forcing_window(tyr, -1.0e30_wp, 1.0e30_wp, k0, k1, msg)
      call expect(w .and. k0 == 1 .and. k1 == 3, 'no restart: the whole axis')
      w = forcing_window(tyr, -1.0e30_wp, 1.0e30_wp, k0, k1, msg, t_restart=1000.0_wp)
      call expect(w .and. k0 == 2 .and. k1 == 3, 'restart at 1000, time_init unset: starts at 1000')
      w = forcing_window(tyr, 1000.0_wp, 1.0e30_wp, k0, k1, msg, t_restart=1000.0_wp)
      call expect(w .and. k0 == 2 .and. k1 == 3, 'restart at 1000, time_init = 1000')
      w = forcing_window(tyr, -1.0e30_wp, 1.0e30_wp, k0, k1, msg, t_restart=2000.0_wp)
      call expect(.not. w, 'final restart (t = 2000), time_init unset: refused, nothing left to run')
      w = forcing_window(tyr, 0.0_wp, 1.0e30_wp, k0, k1, msg, t_restart=1000.0_wp)
      call expect(.not. w, 'restart at 1000, time_init = 0: refused')
      w = forcing_window(tyr, 0.0_wp, 1.0e30_wp, k0, k1, msg, t_restart=2000.0_wp)
      call expect(.not. w, 'final restart, time_init = 0 (the replay): refused')
      w = forcing_window(tyr, -1.0e30_wp, 1.0e30_wp, k0, k1, msg, t_restart=500.0_wp)
      call expect(.not. w, 'restart at 500 (not a slice): refused')
      w = forcing_window(tyr, 0.0_wp, 1.0e30_wp, k0, k1, msg, t_restart=2000.0_wp, restart_reset=.true.)
      call expect(w .and. k0 == 1 .and. k1 == 3, 'restart_time_reset: the window time_init gives')
   end subroutine restart_clock

   subroutine time_axis()
      !! The forcing time axis must be in years (or a unit that converts to them without a
      !! guess) and strictly ascending.
      character(len=400) :: msg
      real(wp) :: sc
      logical  :: w
      w = forcing_time_ok(tyr, 'years', sc, msg);   call expect(w .and. sc == 1.0_wp, 'units "years"')
      w = forcing_time_ok(tyr, '', sc, msg);        call expect(w .and. sc == 1.0_wp, 'no units: years')
      w = forcing_time_ok(tyr, 'Years since 0000-01-01', sc, msg)
      call expect(w .and. sc == 1.0_wp, 'units "years since ..."')
      w = forcing_time_ok(tyr, 'kyr', sc, msg);     call expect(w .and. sc == 1000.0_wp, 'units "kyr": x 1000')
      w = forcing_time_ok(tyr, 'ka BP', sc, msg);   call expect(.not. w, 'units "ka BP": refused')
      w = forcing_time_ok(tyr, 'ka', sc, msg);      call expect(.not. w, 'units "ka": refused')
      w = forcing_time_ok(tyr, 'days since 1950-01-01', sc, msg);  call expect(.not. w, 'units "days since": refused')
      w = forcing_time_ok(tyr, 'seconds', sc, msg); call expect(.not. w, 'units "seconds": refused')
      w = forcing_time_ok(tyr(NT:1:-1), 'years', sc, msg);  call expect(.not. w, 'descending axis: refused')
      w = forcing_time_ok([0.0_wp, 1000.0_wp, 1000.0_wp], 'years', sc, msg)
      call expect(.not. w, 'repeated time: refused')
   end subroutine time_axis

   subroutine expect(cond, what)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: what
      if (cond) then
         write(*,'(a)') '   ok:   '//what
      else
         write(*,'(a)') '   FAIL: '//what;  ok = .false.
      end if
   end subroutine expect

   subroutine write_cfg(cfg, force, ref, out, hor, i_eq, time_end, restart)
      !! A sparse &vilma and a complete &ctl: the legacy Gauss-grid input (remap_input=.false.),
      !! the reference from file_ref (i_eq = 0, the ice-free start slice) or from
      !! z_bed_eq_file / h_ice_eq_file (i_eq = 2, the same fields). time_end [years] and the
      !! restart to resume are optional (default: the whole axis, no restart).
      character(len=*), intent(in) :: cfg, force, ref, out, hor
      integer,          intent(in) :: i_eq
      real(wp),         intent(in), optional :: time_end
      character(len=*), intent(in), optional :: restart
      integer :: u
      open(newunit=u, file=cfg, status="replace", action="write")
      write(u,'(a)')    "&vilma"
      write(u,'(a,i0)') "    lmax = ", LMAX
      write(u,'(a,i0)') "    nlat = ", NLAT
      write(u,'(a,i0)') "    nphi = ", NPHI
      write(u,'(a)')    "/"
      write(u,'(a)')    "&ctl"
      write(u,'(a)')    '    file_forcing = "'//force//'"'
      write(u,'(a)')    '    name_ice     = "h_ice"'
      write(u,'(a)')    '    name_time    = "time"'
      write(u,'(a)')    '    file_ref     = "'//ref//'"'
      write(u,'(a)')    '    name_zbed_eq = "z_bed_eq"'
      write(u,'(a)')    '    name_hice_ref = "h_ice_eq"'
      write(u,'(a)')    '    file_out     = "'//out//'"'
      write(u,'(a)')    '    file_hor     = "'//hor//'"'
      write(u,'(a)')    "    write_lossless = .true."
      write(u,'(a)')    "    time_init    = -1.0e30"
      if (present(time_end)) then
         write(u,'(a,es12.4)') "    time_end     = ", time_end
      else
         write(u,'(a)') "    time_end     =  1.0e30"
      end if
      write(u,'(a)')    "    remap_input  = .false."   ! input already on the Gauss grid
      write(u,'(a)')    '    name_lon     = "lon"'
      write(u,'(a)')    '    name_lat     = "lat"'
      write(u,'(a,i0)') "    i_eq         = ", i_eq
      write(u,'(a)')    '    z_bed_ref_file   = ""'
      write(u,'(a)')    '    h_ice_ref_file   = ""'
      if (i_eq == 2) then
         write(u,'(a)') '    z_bed_eq_file    = "'//ref//'"'
         write(u,'(a)') '    h_ice_eq_file    = "'//ref//'"'
         write(u,'(a)') '    name_z_bed_ref = "z_bed_eq"'
         write(u,'(a)') '    name_h_ice_ref = "h_ice_eq"'
      else
         write(u,'(a)') '    z_bed_eq_file    = ""'
         write(u,'(a)') '    h_ice_eq_file    = ""'
         write(u,'(a)') '    name_z_bed_ref = "bedrock_topography"'
         write(u,'(a)') '    name_h_ice_ref = "ice_thickness"'
      end if
      write(u,'(a)')    '    rsl_restart_file = ""'
      write(u,'(a)')    '    name_rsl       = "rsl"'
      if (present(restart)) then
         write(u,'(a)') '    restart_in_file = "'//restart//'"'
      else
         write(u,'(a)') '    restart_in_file = ""'
      end if
      write(u,'(a)')    "/"
      close(u)
   end subroutine write_cfg

   subroutine same_output(file, what)
      !! rsl, z_bed and dz_bed of file (all slices) equal those of OUT, bit for bit.
      character(len=*), intent(in) :: file, what
      character(len=8), parameter :: VARS(3) = [character(len=8) :: "rsl", "z_bed", "dz_bed"]
      integer :: q
      if (nc_size(file, "time") /= NT) then
         write(*,'(a)') '   FAIL: '//what//': expected one output slice per forcing slice'; ok = .false.;  return
      end if
      do q = 1, size(VARS)
         call nc_read(OUT,  trim(VARS(q)), a3)
         call nc_read(file, trim(VARS(q)), b3)
         write(*,'(a,es10.2)') '   '//what//': max|'//trim(VARS(q))//' - north-first run| =', maxval(abs(b3 - a3))
         if (any(b3 /= a3)) then
            write(*,'(a)') '   FAIL: '//what//' differs from the north-first files'; ok = .false.
         end if
      end do
   end subroutine same_output

   integer function nearest_row(colat_deg) result(jbest)
      real(wp), intent(in) :: colat_deg
      integer  :: j
      real(wp) :: d, dbest
      jbest = 1;  dbest = huge(1.0_wp)
      do j = 1, sht%nlat
         d = abs(sht%colat(j)*rad2deg - colat_deg)
         if (d < dbest) then;  dbest = d;  jbest = j;  end if
      end do
   end function nearest_row

end program test_drive
