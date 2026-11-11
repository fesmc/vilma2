program vilma_main
   !! Standalone VILMA2 driver.
   !!
   !!   ./bin/vilma.x [run-config.nml [defaults.nml]]   (default: vilma.nml, input/vilma_defaults.nml)
   !!
   !! The run config supplies the &ctl run-control group and any &vilma physics
   !! overrides; its &vilma is overlaid on the complete physics defaults in
   !! input/vilma_defaults.nml (loaded automatically — see DEFAULTS_FILE).
   use vilma, only: vilma_run, DEFAULTS_FILE
   implicit none
   character(len=512) :: cfg, defs

   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "vilma.nml"
   end if

   defs = DEFAULTS_FILE
   if (command_argument_count() >= 2) call get_command_argument(2, defs)
   call vilma_run(trim(cfg), trim(defs))
end program vilma_main
