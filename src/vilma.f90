module vilma
   !! VILMA2 host API: the one module a host model `use`s. It re-exports the model's
   !! contract and nothing of its internals:
   !!
   !!   the parameter record   vilma_param_class, vilma_par_load, vilma_par_print, DEFAULTS_FILE
   !!   the model              solid_earth and its life cycle: solid_earth_init, _update,
   !!                          _finalize, _spinup, _set_time, _sync_host, _enable_visc_3d
   !!   checks                 solid_earth_check_solver, solid_earth_host_check
   !!   restarts               vilma_restart_write, vilma_restart_read, vilma_restart_time,
   !!                          vilma_restart_compatible, and the step writers
   !!   constants              wp (the real kind), VILMA_UNSET (a value the backend did not set)
   !!   the standalone driver  vilma_run
   !!
   !!   use vilma
   !!   type(solid_earth) :: se
   !!   call vilma_par_load(se%par, "my.nml", defaults_file=DEFAULTS_FILE)
   !!   call solid_earth_init(se, z_bed_eq, h_ice_eq, grid=host_grid, time=t0)
   !!   call solid_earth_update(se, h_ice, dt_yr)      ! se%rsl, se%z_bed, se%C, se%bsl, ...
   !!   call solid_earth_finalize(se)
   use vilma_precision, only: wp
   use vilma_params,    only: vilma_param_class, vilma_par_load, vilma_par_print
   use vilma_control,   only: DEFAULTS_FILE
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize, &
                              solid_earth_spinup, solid_earth_set_time, solid_earth_sync_host, &
                              solid_earth_enable_visc_3d, solid_earth_check_solver, solid_earth_host_check, &
                              VILMA_UNSET
   use vilma_io,        only: vilma_restart_write, vilma_restart_read, vilma_restart_time, &
                              vilma_restart_compatible, vilma_write_step, vilma_write_horizontal, &
                              vilma_io_set_table
   use vilma_drive,     only: vilma_run
   implicit none
   private

   public :: wp, VILMA_UNSET, DEFAULTS_FILE, vilma_version
   public :: vilma_param_class, vilma_par_load, vilma_par_print
   public :: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   public :: solid_earth_spinup, solid_earth_set_time, solid_earth_sync_host, solid_earth_enable_visc_3d
   public :: solid_earth_check_solver, solid_earth_host_check
   public :: vilma_restart_write, vilma_restart_read, vilma_restart_time, vilma_restart_compatible
   public :: vilma_write_step, vilma_write_horizontal, vilma_io_set_table
   public :: vilma_run

#ifndef VERSION
#define VERSION "unknown"
#endif
   character(len=*), parameter :: vilma_version = VERSION

end module vilma
