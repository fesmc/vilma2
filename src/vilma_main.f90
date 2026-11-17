program vilma_main
   !! Standalone VILMA2 driver.
   !!
   !!   ./bin/vilma.x [run-config.nml [defaults.nml]]   (default: vilma.nml, input/vilma_defaults.nml)
   !!
   !! The run config supplies the &ctl run-control group and any &vilma physics
   !! overrides; its &vilma is overlaid on the complete physics defaults in
   !! input/vilma_defaults.nml (loaded automatically — see DEFAULTS_FILE).
   !!
   !! Exit code 75 (EX_TEMPFAIL): the run stopped before the end of its batch job
   !! (SLURM_JOB_END_TIME, vilma_drive), with final/vilma_restart.nc to resume it from.
   use vilma, only: vilma_run, DEFAULTS_FILE, VILMA_WALL_STOP
   implicit none
   character(len=512) :: cfg, defs
   integer :: status

   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "vilma.nml"
   end if

   defs = DEFAULTS_FILE
   if (command_argument_count() >= 2) call get_command_argument(2, defs)
   call vilma_run(trim(cfg), trim(defs), status)
   if (status == VILMA_WALL_STOP) stop 75
end program vilma_main
