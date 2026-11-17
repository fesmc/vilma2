program test_drive
   !! Standalone driver (vilma_drive): generate a synthetic reference state and an
   !! ice-thickness forcing on the model Gauss grid, write a sparse &vilma config
   !! overlaid on the shipped defaults, run vilma_run, and check the output:
   !!   (1) one output slice per forcing slice,
   !!   (2) the bed subsides under the growing ice cap (z_bed < z_bed_eq),
   !!   (3) the ocean draws down (rsl < 0 in the far field),
   !!   (4) the uplift rate dz_bed_rate of a slice is the change of dz_bed since the slice
   !!       before over the interval (the file is lossless: to round-off), < 0 under the cap.
   !! Then the legacy Gauss-grid reader (remap_input=.false.) with the same files stored
   !! with latitude ascending (rows reversed): read with i_eq = 0 (file_ref) and i_eq = 2
   !! (z_bed_eq_file, h_ice_eq_file), the output must equal that of the north-first files.
   !! Then the restart clock: a run stopped at the middle slice and resumed from its final
   !! restart continues there (one slice left), and forcing_window refuses every window that
   !! would replay or skip forcing. The two legs are the uninterrupted run bit for bit, every
   !! output variable and slice (the resumed leg's first slice is the first leg's last). So
   !! are the legs of a run with a spin-up (under a deglaciating cap): the resumed leg runs
   !! neither the spin-up nor the seed again, and a run resumed from the spin-up restart is
   !! the whole run. Then the wall-time stop: with SLURM_JOB_END_TIME already past, the run
   !! stops after its first slice with status VILMA_WALL_STOP and a final restart, and the
   !! run resumed from it is the uninterrupted one bit for bit; with the variable not a
   !! number or far ahead, the run goes to the end.
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use vilma_precision,       only: wp
   use vilma_constants,       only: rad2deg, pi
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
   use, intrinsic :: iso_fortran_env, only: int64
   use vilma_drive,           only: vilma_run, forcing_window, forcing_time_ok, input_field_ok, unix_now, &
                                    VILMA_DONE, VILMA_WALL_STOP
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
   character(len=*), parameter :: OUTL0 = "libvilma/include/test_drive_out_leg0.nc"
   character(len=*), parameter :: FORCED = "libvilma/include/test_drive_force_deglac.nc"
   character(len=*), parameter :: OUTS0 = "libvilma/include/test_drive_out_spin0.nc"
   character(len=*), parameter :: OUTS1 = "libvilma/include/test_drive_out_spin1.nc"
   character(len=*), parameter :: OUTS2 = "libvilma/include/test_drive_out_spin2.nc"
   character(len=*), parameter :: OUTS3 = "libvilma/include/test_drive_out_spin3.nc"
   character(len=*), parameter :: OUTSN = "libvilma/include/test_drive_out_nospin.nc"
   character(len=*), parameter :: OUTW0 = "libvilma/include/test_drive_out_wall0.nc"
   character(len=*), parameter :: OUTW1 = "libvilma/include/test_drive_out_wall1.nc"
   character(len=*), parameter :: OUTW2 = "libvilma/include/test_drive_out_wall2.nc"
   character(len=*), parameter :: RSTW  = "libvilma/include/test_drive_rst_wall.nc"
   character(len=*), parameter :: RSTS0 = "libvilma/include/test_drive_rst_spin0.nc"
   character(len=*), parameter :: RSTS1 = "libvilma/include/test_drive_rst_spin1.nc"
   !! every output variable of the native solver with rotation: (lon, lat, time) fields and scalars
   character(len=11), parameter :: OUT2D(7) = [character(len=11) :: "h_ice", "rsl", "z_bed", "C_ocean", &
        "dz_bed", "dz_bed_rate", "dz_ss"]
   character(len=11), parameter :: OUT1D(8) = [character(len=11) :: "time", "bsl", "bsl_ive", "dphi", &
        "sle_resid", "sle_n_cap", "rot_m_re", "rot_m_im"]

   type(sht_grid), target :: sht
   real(wp), allocatable  :: lon_deg(:), lat_deg(:), z_bed_eq(:,:), h_ice_eq(:,:)
   real(wp), allocatable  :: h_ice(:,:,:), tyr(:), zb(:,:), rsl(:,:)
   real(wp), allocatable  :: ue(:,:), un(:,:), uet(:,:), unt(:,:)
   real(wp), allocatable  :: a3(:,:,:), b3(:,:,:)
   real(wp), allocatable  :: ub0(:,:), ub1(:,:), ubr(:,:)
   real(wp) :: thd, sub, err
   integer  :: i, j, k, jice, jocean, nout
   logical  :: ok

   interface   ! POSIX: the test sets SLURM_JOB_END_TIME itself
      integer(c_int) function c_setenv(name, value, overwrite) bind(C, name="setenv")
         import :: c_int, c_char
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
      end function c_setenv
      integer(c_int) function c_unsetenv(name) bind(C, name="unsetenv")
         import :: c_int, c_char
         character(kind=c_char), intent(in) :: name(*)
      end function c_unsetenv
   end interface

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

   ! --- (4) the uplift rate of the last slice: d(dz_bed)/dt over its interval ------------
   allocate(ub0(sht%nphi,sht%nlat), ub1(sht%nphi,sht%nlat), ubr(sht%nphi,sht%nlat))
   call nc_read(OUT, "dz_bed",      ub0, start=[1,1,nout-1], count=[sht%nphi,sht%nlat,1])
   call nc_read(OUT, "dz_bed",      ub1, start=[1,1,nout],   count=[sht%nphi,sht%nlat,1])
   call nc_read(OUT, "dz_bed_rate", ubr, start=[1,1,nout],   count=[sht%nphi,sht%nlat,1])
   err = maxval(abs(ubr - (ub1 - ub0)/(tyr(NT) - tyr(NT-1))))/maxval(abs(ubr))
   write(*,'(a,es10.2,a,es10.2)') '   dz_bed_rate under ice =', ubr(1,jice), ' m/yr; |rate - d(dz_bed)/dt|/max|rate| =', err
   if (err > 1.0e-13_wp .or. ubr(1,jice) >= 0.0_wp) then
      write(*,'(a)') '   FAIL: dz_bed_rate is not the subsidence rate of the last interval'; ok = .false.
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
   call wall_stop()
   call time_axis()
   call field_values()

   call sht_grid_destroy(sht)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: vilma_run reads ref+forcing, marches the model,'
      write(*,'(a)') '       and writes a sensible subsiding/draw-down output'
   else
      write(*,'(a)') ' FAIL: standalone driver did not all pass'
      error stop 1
   end if

contains

   subroutine restart_clock()
      !! A restart continues where it left off. Leg 1 stops at the middle slice (t = 1000 yr);
      !! leg 2 resumes its final restart with time_init unset and must march the one
      !! remaining interval, 1000 -> 2000 yr, to the end state of the uninterrupted run.
      !! (Before the clock was checked, leg 2 replayed the forcing from t = 0 on the
      !! restart's memory.) The legs are the uninterrupted run bit for bit. Then the same
      !! with a spin-up, and the window rule itself, case by case.
      real(wp), allocatable :: t2(:)
      character(len=400) :: msg
      integer :: k0, k1
      logical :: w
      call write_cfg(CFGA, FORCE, REF, OUTL0, "", 2)
      call vilma_run(CFGA, defaults_file=DEFS)
      call write_cfg(CFGA, FORCE, REF, OUTL1, "", 2, time_end=1000.0_wp)
      call vilma_run(CFGA, defaults_file=DEFS)
      call execute_command_line("cp libvilma/include/final/vilma_restart.nc "//RST)
      call write_cfg(CFGA, FORCE, REF, OUTL2, "", 2, restart=RST)
      call vilma_run(CFGA, defaults_file=DEFS)
      allocate(t2(nc_size(OUTL2, "time")))
      call nc_read(OUTL2, "time", t2)
      write(*,'(a,i0,a,2f10.1)') '   resumed run: ', size(t2), ' slices, t first/last =', t2(1), t2(size(t2))
      if (size(t2) /= 2 .or. t2(1) /= 1000.0_wp .or. t2(size(t2)) /= 2000.0_wp) then
         write(*,'(a)') '   FAIL: the resumed run does not continue at the restart time'; ok = .false.
      else
         call same_slices(OUTL0, 1, OUTL1, 1, 2, 'leg 1')
         call same_slices(OUTL0, 2, OUTL2, 1, 2, 'leg 2 (resumed)')
      end if

      ! the same with a spin-up under the start-slice ice: a cap deglaciating 2 -> 1 -> 0 km
      call nc_create(FORCED, overwrite=.true.)
      call nc_write_dim(FORCED, "lon",  x=lon_deg, units="degrees_east")
      call nc_write_dim(FORCED, "lat",  x=lat_deg, units="degrees_north")
      call nc_write_dim(FORCED, "time", x=tyr,     units="years", unlimited=.true.)
      call nc_write(FORCED, "h_ice", h_ice(:, :, NT:1:-1), dim1="lon", dim2="lat", dim3="time")
      call write_cfg(CFGA, FORCED, REF, OUTS0, "", 2, spinup=.true.)
      call vilma_run(CFGA, defaults_file=DEFS)
      call write_cfg(CFGA, FORCED, REF, OUTS1, "", 2, time_end=1000.0_wp, spinup=.true.)
      call vilma_run(CFGA, defaults_file=DEFS)
      call execute_command_line("cp libvilma/include/spinup/vilma_restart.nc "//RSTS0)
      call execute_command_line("cp libvilma/include/final/vilma_restart.nc "//RSTS1)
      call write_cfg(CFGA, FORCED, REF, OUTS2, "", 2, restart=RSTS1, spinup=.true.)
      call vilma_run(CFGA, defaults_file=DEFS)
      call write_cfg(CFGA, FORCED, REF, OUTS3, "", 2, restart=RSTS0, spinup=.true.)
      call vilma_run(CFGA, defaults_file=DEFS)
      call same_slices(OUTS0, 1, OUTS1, 1, 2, 'spin-up, leg 1')
      call same_slices(OUTS0, 2, OUTS2, 1, 2, 'spin-up, leg 2 (resumed)')
      call same_slices(OUTS0, 1, OUTS3, 1, 3, 'spin-up, resumed from the spin-up restart')
      ! the spin-up did something: under the cap the bed starts deeper than without it (elastic)
      call write_cfg(CFGA, FORCED, REF, OUTSN, "", 2)
      call vilma_run(CFGA, defaults_file=DEFS)
      block
         real(wp), allocatable :: z0(:,:), ze(:,:)
         integer :: jc
         allocate(z0(sht%nphi,sht%nlat), ze(sht%nphi,sht%nlat))
         jc = nearest_row(15.0_wp)
         call nc_read(OUTS0, "dz_bed", z0, start=[1,1,1], count=[sht%nphi,sht%nlat,1])
         call nc_read(OUTSN, "dz_bed", ze, start=[1,1,1], count=[sht%nphi,sht%nlat,1])
         write(*,'(a,2f10.3,a)') '   dz_bed under the cap at the start, spun up / elastic =', z0(1,jc), ze(1,jc), ' m'
         call expect(z0(1,jc) < ze(1,jc) - 1.0_wp, 'spin-up: the bed starts deeper under the cap than the elastic one')
      end block

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

   subroutine wall_stop()
      !! The wall-time stop (vilma_drive): SLURM_JOB_END_TIME set by the test. The variable of
      !! the job the test may run in is put back at the end.
      character(len=64) :: job, v
      integer :: st, n, u, ios
      integer(int64) :: t_date
      logical :: had
      call get_environment_variable("SLURM_JOB_END_TIME", job, length=n, status=st)
      had = (st == 0)
      ! the clock: unix_now is date +%s
      call execute_command_line("date +%s > libvilma/include/test_drive_date.txt")
      open(newunit=u, file="libvilma/include/test_drive_date.txt", action="read")
      read(u, *, iostat=ios) t_date;  close(u)
      write(*,'(a,i0,a,i0)') '   unix_now = ', unix_now(), ', date +%s = ', t_date
      call expect(ios == 0 .and. abs(unix_now() - t_date) <= 2, 'unix_now is the Unix time (date +%s)')
      ! the job ends now: the run stops after its first slice, with the final restart
      write(v,'(i0)') unix_now()
      call set_end(trim(v))
      call write_cfg(CFGA, FORCE, REF, OUTW1, "", 2)
      call vilma_run(CFGA, defaults_file=DEFS, status=st)
      call expect(st == VILMA_WALL_STOP .and. nc_size(OUTW1, "time") == 2, &
                  'job end past: stops after the first slice (status VILMA_WALL_STOP)')
      call same_slices(OUTL0, 1, OUTW1, 1, 2, 'wall stop')
      call execute_command_line("cp libvilma/include/final/vilma_restart.nc "//RSTW)
      ! not a number: no stop
      call set_end("soon")
      call write_cfg(CFGA, FORCE, REF, OUTW0, "", 2)
      call vilma_run(CFGA, defaults_file=DEFS, status=st)
      call expect(st == VILMA_DONE .and. nc_size(OUTW0, "time") == NT, 'SLURM_JOB_END_TIME not a number: runs to the end')
      ! the job ends in a day: the resumed run goes to the end, and is the uninterrupted one
      write(v,'(i0)') unix_now() + 86400_int64
      call set_end(trim(v))
      call write_cfg(CFGA, FORCE, REF, OUTW2, "", 2, restart=RSTW)
      call vilma_run(CFGA, defaults_file=DEFS, status=st)
      call expect(st == VILMA_DONE .and. nc_size(OUTW2, "time") == 2, 'job end ahead: the resumed run goes to the end')
      call same_slices(OUTL0, 2, OUTW2, 1, 2, 'resumed after the wall stop')
      if (had) then
         call set_end(job(1:n))
      else
         st = c_unsetenv("SLURM_JOB_END_TIME"//c_null_char)
      end if
   end subroutine wall_stop

   subroutine set_end(v)
      character(len=*), intent(in) :: v
      if (c_setenv("SLURM_JOB_END_TIME"//c_null_char, v//c_null_char, 1_c_int) /= 0) &
         error stop 'test_drive: setenv failed'
   end subroutine set_end

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

   subroutine field_values()
      !! An input field with a fill value, a NaN or a value outside its range is refused.
      real(wp), allocatable :: f(:,:)
      character(len=200) :: msg
      logical :: w
      allocate(f(4,3))
      f = 1500.0_wp;                             w = input_field_ok(f, 0.0_wp, 1.0e4_wp, msg)
      call expect(w, 'ice 1500 m: accepted')
      f = 1500.0_wp;  f(2,2) = -9999.0_wp;       w = input_field_ok(f, 0.0_wp, 1.0e4_wp, msg)
      call expect(.not. w, 'ice with -9999 (ncio missing value): refused')
      f = 1500.0_wp;  f(3,1) = 9.96921e36_wp;    w = input_field_ok(f, 0.0_wp, 1.0e4_wp, msg)
      call expect(.not. w, 'ice with a netCDF fill value: refused')
      f = 1500.0_wp;  f(1,3) = ieee_value(1.0_wp, ieee_quiet_nan);  w = input_field_ok(f, 0.0_wp, 1.0e4_wp, msg)
      call expect(.not. w, 'ice with NaN: refused')
      f = -4000.0_wp;                            w = input_field_ok(f, -2.0e4_wp, 2.0e4_wp, msg)
      call expect(w, 'bed -4000 m: accepted')
   end subroutine field_values

   subroutine expect(cond, what)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: what
      if (cond) then
         write(*,'(a)') '   ok:   '//what
      else
         write(*,'(a)') '   FAIL: '//what;  ok = .false.
      end if
   end subroutine expect

   subroutine write_cfg(cfg, force, ref, out, hor, i_eq, time_end, restart, spinup)
      !! A sparse &vilma and a complete &ctl: the legacy Gauss-grid input (remap_input=.false.),
      !! the reference from file_ref (i_eq = 0, the ice-free start slice) or from
      !! z_bed_eq_file / h_ice_eq_file (i_eq = 2, the same fields). time_end [years], the
      !! restart to resume and a spin-up under the start-slice ice (equil_time_max) are
      !! optional (default: the whole axis, no restart, no spin-up).
      character(len=*), intent(in) :: cfg, force, ref, out, hor
      integer,          intent(in) :: i_eq
      real(wp),         intent(in), optional :: time_end
      character(len=*), intent(in), optional :: restart
      logical,          intent(in), optional :: spinup
      integer :: u
      open(newunit=u, file=cfg, status="replace", action="write")
      write(u,'(a)')    "&vilma"
      write(u,'(a,i0)') "    lmax = ", LMAX
      write(u,'(a,i0)') "    nlat = ", NLAT
      write(u,'(a,i0)') "    nphi = ", NPHI
      if (present(spinup)) then
         if (spinup) write(u,'(a)') "    equil_time_max = 2.0e4"
      end if
      write(u,'(a)')    "/"
      write(u,'(a)')    "&ctl"
      write(u,'(a)')    '    file_forcing = "'//force//'"'
      write(u,'(a)')    '    name_h_ice   = "h_ice"'
      write(u,'(a)')    '    name_time    = "time"'
      write(u,'(a)')    '    file_ref     = "'//ref//'"'
      write(u,'(a)')    '    name_z_bed_start = "z_bed_eq"'
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
      !! rsl, z_bed, dz_bed and dz_bed_rate of file (all slices) equal those of OUT, bit for bit.
      character(len=*), intent(in) :: file, what
      character(len=11), parameter :: VARS(4) = [character(len=11) :: "rsl", "z_bed", "dz_bed", "dz_bed_rate"]
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

   subroutine same_slices(fa, ia, fb, ib, n, what)
      !! Slices ib .. ib+n-1 of fb equal slices ia .. ia+n-1 of fa, bit for bit, in every
      !! output variable (OUT2D, OUT1D).
      character(len=*), intent(in) :: fa, fb, what
      integer,          intent(in) :: ia, ib, n
      real(wp), allocatable :: a(:,:,:), b(:,:,:), a1(:), b1(:)
      integer :: q, nbad
      nbad = 0
      allocate(a(sht%nphi,sht%nlat,n), b(sht%nphi,sht%nlat,n), a1(n), b1(n))
      do q = 1, size(OUT2D)
         call nc_read(fa, trim(OUT2D(q)), a, start=[1,1,ia], count=[sht%nphi,sht%nlat,n])
         call nc_read(fb, trim(OUT2D(q)), b, start=[1,1,ib], count=[sht%nphi,sht%nlat,n])
         if (any(a /= b)) then
            write(*,'(a,es10.2)') '   '//what//': '//trim(OUT2D(q))//' differs, max|d| =', maxval(abs(a - b))
            nbad = nbad + 1
         end if
      end do
      do q = 1, size(OUT1D)
         call nc_read(fa, trim(OUT1D(q)), a1, start=[ia], count=[n])
         call nc_read(fb, trim(OUT1D(q)), b1, start=[ib], count=[n])
         if (any(a1 /= b1)) then
            write(*,'(a,es10.2)') '   '//what//': '//trim(OUT1D(q))//' differs, max|d| =', maxval(abs(a1 - b1))
            nbad = nbad + 1
         end if
      end do
      call expect(nbad == 0, what//': the uninterrupted run bit for bit (every variable and slice)')
   end subroutine same_slices

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
