program test_params
   !! vilma_params: load a &vilma namelist into vilma_param_class and build the earth
   !! model from it. A sparse user file is overlaid on the complete, shipped
   !! defaults file (input/vilma_defaults.nml) — the yelmo defaults_file
   !! convention. Checks
   !! (1) scalar / string / logical overrides, (2) the YEARS->seconds conversion
   !! of the time fields, (3) custom per-layer earth assembly, (4) the named
   !! built-in path (defaults), (5) an un-overridden value falling through to
   !! the defaults file, (6) path expansion, and (7) an unknown key stopping the load.
   use vilma_precision,       only: wp
   use vilma_constants,       only: sec_per_year
   use vilma_params,          only: vilma_param_class, vilma_par_load, expand_path, scheme_name_valid, cfl_valid
   use vilma_earth_structure, only: earth_n_layers, earth_model, build_earth, RHEOL_MAXWELL, RHEOL_FLUID
   implicit none

   character(len=*), parameter :: NML  = "libvilma/include/test_params.nml"
   character(len=*), parameter :: NML_BAD = "libvilma/include/test_params_bad.nml"
   character(len=*), parameter :: DEFS = "input/vilma_defaults.nml"   ! shipped complete &vilma defaults
   type(vilma_param_class) :: p, pdef, pfile
   type(earth_model)    :: em, emdef
   integer :: u, ex
   character(len=256) :: arg, exe
   logical :: ok
   character(len=512) :: home
   real(wp) :: tol

   ok = .true.;  tol = 1.0e-6_wp

   ! a child run: load a user file with a key the defaults do not declare; reaching
   ! the end is the failure (exit status 0)
   call get_command_argument(1, arg)
   if (trim(arg) == "unknown-key") then
      open(newunit=u, file=NML_BAD, status="replace", action="write")
      write(u,'(a)') "&vilma"
      write(u,'(a)') "    lmax = 32"
      write(u,'(a)') "    no_such_parameter = 1"
      write(u,'(a)') "/"
      close(u)
      call vilma_par_load(p, NML_BAD, defaults_file=DEFS)
      stop 0
   end if

   ! --- write a SPARSE user file (custom 2-layer earth) overriding the defaults -
   open(newunit=u, file=NML, status="replace", action="write")
   write(u,'(a)') "&vilma"
   write(u,'(a)') "    lmax    = 32"
   write(u,'(a)') '    earth   = "custom"'
   write(u,'(a)') "    n_layer = 2"
   write(u,'(a)') "    r_bot    = 5701.0e3, 0.0"
   write(u,'(a)') "    r_top    = 6371.0e3, 5701.0e3"
   write(u,'(a)') "    rho      = 4000.0, 10000.0"
   write(u,'(a)') "    mu       = 1.0e11, 0.0"
   write(u,'(a)') "    eta      = 1.0e21, 0.0"
   write(u,'(a)') "    rheology = 1, 2"
   write(u,'(a)') '    scheme       = "trap"'
   write(u,'(a)') "    sle_n_outer  = 5"
   write(u,'(a)') "    sle_subgrid  = .false."
   write(u,'(a)') "    equil_time_max = 2000.0"
   write(u,'(a)') "    dt_init      = 500.0"
   write(u,'(a)') "    l_rotation   = .true."
   write(u,'(a)') "    visc_3d_structure = 6"
   write(u,'(a)') "/"
   close(u)

   call vilma_par_load(p, NML, defaults_file=DEFS)

   ! --- (1) scalar / string / logical overrides -------------------------------
   call check_int("lmax",        p%lmax,        32)
   call check_str("earth",       trim(p%earth), "custom")
   call check_int("n_layer",     p%n_layer,     2)
   call check_str("scheme",      trim(p%scheme), "trap")
   ! the schemes the model runs; the kernel-only rules are not a valid setting
   call check_log("scheme fe valid",    scheme_name_valid("fe"),   .true.)
   call check_log("scheme trap valid",  scheme_name_valid("trap"), .true.)
   call check_log("scheme etd1 valid",  scheme_name_valid("etd1"), .false.)
   call check_log("scheme be valid",    scheme_name_valid("be"),   .false.)
   ! the fe Maxwell-number ceiling: stable only for 0 < M < 2
   call check_log("cfl 1 valid",    cfl_valid(1.0_wp),  .true.)
   call check_log("cfl 1.9 valid",  cfl_valid(1.9_wp),  .true.)
   call check_log("cfl 2 valid",    cfl_valid(2.0_wp),  .false.)
   call check_log("cfl 2.5 valid",  cfl_valid(2.5_wp),  .false.)
   call check_log("cfl 0 valid",    cfl_valid(0.0_wp),  .false.)
   call check_int("sle_n_outer", p%sle_n_outer, 5)
   call check_log("sle_subgrid", p%sle_subgrid, .false.)
   call check_log("l_rotation",  p%l_rotation,    .true.)
   call check_int("visc_3d_structure", p%visc_3d_structure, 6)

   ! --- (2) YEARS -> seconds conversion of the time fields ---------------------
   call check_real("equil_time_max [s]", p%equil_time_max, 2000.0_wp*sec_per_year)
   call check_real("dt_init [s]",        p%dt_init,         500.0_wp*sec_per_year)

   ! --- (3) custom per-layer earth assembly ------------------------------------
   em = build_earth(p)
   call check_int("custom n_layers", earth_n_layers(em), 2)
   call check_str("custom name",     em%name,       "custom")
   call check_int("layer1 rheology", em%layers(1)%rheology, RHEOL_MAXWELL)
   call check_int("layer2 rheology", em%layers(2)%rheology, RHEOL_FLUID)
   call check_real("layer1 mu [Pa]", em%layers(1)%mu, 1.0e11_wp)

   ! --- (4) named built-in path (in-code defaults) -----------------------------
   emdef = build_earth(pdef)                     ! pdef left at defaults => "M3-L70-V01"
   call check_str("default earth",    trim(pdef%earth), "M3-L70-V01")
   call check_int("M3-L70-V01 layers", earth_n_layers(emdef), 5)

   ! --- (5) un-overridden values fall through to the defaults file -------------
   call check_int("fallthrough sle_max_mem_iter", p%sle_max_mem_iter, 20)
   call check_int("fallthrough sle_n_inner",     p%sle_n_inner,     20)
   call check_int("default visc_3d_structure",   pdef%visc_3d_structure, 6)

   ! --- (5b) one source of truth: the shipped defaults file IS the in-code defaults ----
   ! (a host fills the record in memory and never reads the file). lmax (required) and
   ! the custom-layer arrays (inert unless earth = "custom") are the file's own.
   call vilma_par_load(pfile, DEFS)
   call same_defaults(pfile, pdef)

   ! --- (6) ~/ and $HOME/ expansion in file paths ------------------------------
   call get_environment_variable("HOME", home)
   call check_str("expand ~/",     trim(expand_path("~/x/y.nc")),     trim(home)//"/x/y.nc")
   call check_str("expand $HOME/", trim(expand_path("$HOME/x/y.nc")), trim(home)//"/x/y.nc")
   call check_str("expand rel",    trim(expand_path("input/y.nc")),   "input/y.nc")
   call check_str("expand abs",    trim(expand_path("/data/y.nc")),   "/data/y.nc")
   call check_str("expand empty",  trim(expand_path("")),             "")

   ! --- (7) a key the defaults file does not declare stops the load (child process) ---
   call get_command_argument(0, exe)
   call execute_command_line(trim(exe)//' unknown-key > /dev/null 2>&1', exitstat=ex)
   call check_log("unknown key stops", ex /= 0, .true.)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: vilma_par_load fills the record (years->s) and build_earth'
      write(*,'(a)') '       assembles both custom and named earth models'
   else
      write(*,'(a)') ' FAIL: vilma_params did not all pass'
      error stop 1
   end if

contains

   subroutine same_defaults(a, b)
      !! Every parameter that acts without earth = "custom" agrees between a and b.
      type(vilma_param_class), intent(in) :: a, b
      call check_str("defaults: solver", trim(a%solver), trim(b%solver))
      call check_int("defaults: nlat", a%nlat, b%nlat);  call check_int("defaults: nphi", a%nphi, b%nphi)
      call check_int("defaults: mmax", a%mmax, b%mmax);  call check_int("defaults: mres", a%mres, b%mres)
      call check_real("defaults: eps_polar", a%eps_polar, b%eps_polar)
      call check_str("defaults: earth", trim(a%earth), trim(b%earth))
      call check_str("defaults: scheme", trim(a%scheme), trim(b%scheme))
      call check_str("defaults: earth_response", trim(a%earth_response), trim(b%earth_response))
      call check_int("defaults: sle_n_outer", a%sle_n_outer, b%sle_n_outer)
      call check_int("defaults: sle_n_inner", a%sle_n_inner, b%sle_n_inner)
      call check_real("defaults: sle_tol", a%sle_tol, b%sle_tol)
      call check_int("defaults: sle_max_mem_iter", a%sle_max_mem_iter, b%sle_max_mem_iter)
      call check_log("defaults: sle_fixed_ocean", a%sle_fixed_ocean, b%sle_fixed_ocean)
      call check_log("defaults: sle_subgrid", a%sle_subgrid, b%sle_subgrid)
      call check_real("defaults: rho_ice", a%rho_ice, b%rho_ice)
      call check_real("defaults: rho_water", a%rho_water, b%rho_water)
      call check_real("defaults: dt_init", a%dt_init, b%dt_init)
      call check_real("defaults: dt_min", a%dt_min, b%dt_min)
      call check_real("defaults: rtol", a%rtol, b%rtol);  call check_real("defaults: atol", a%atol, b%atol)
      call check_real("defaults: safety", a%safety, b%safety)
      call check_real("defaults: grow_max", a%grow_max, b%grow_max)
      call check_real("defaults: shrink_min", a%shrink_min, b%shrink_min)
      call check_real("defaults: cfl", a%cfl, b%cfl)
      call check_log("defaults: l_rotation", a%l_rotation, b%l_rotation)
      call check_real("defaults: rotation_k_s", a%rotation_k_s, b%rotation_k_s)
      call check_real("defaults: rotation_c_minus_a", a%rotation_c_minus_a, b%rotation_c_minus_a)
      call check_real("defaults: equil_time_max", a%equil_time_max, b%equil_time_max)
      call check_real("defaults: equil_rate_tol", a%equil_rate_tol, b%equil_rate_tol)
      call check_log("defaults: pre_spinup_1d", a%pre_spinup_1d, b%pre_spinup_1d)
      call check_str("defaults: deg1_frame", trim(a%deg1_frame), trim(b%deg1_frame))
      call check_log("defaults: l_visc_3d", a%l_visc_3d, b%l_visc_3d)
      call check_str("defaults: visc_3d_file", trim(a%visc_3d_file), trim(b%visc_3d_file))
      call check_int("defaults: visc_3d_structure", a%visc_3d_structure, b%visc_3d_structure)
      call check_str("defaults: name_visc", trim(a%name_visc), trim(b%name_visc))
      call check_str("defaults: name_visc_lon", trim(a%name_visc_lon), trim(b%name_visc_lon))
      call check_str("defaults: name_visc_lat", trim(a%name_visc_lat), trim(b%name_visc_lat))
      call check_str("defaults: name_visc_r", trim(a%name_visc_r), trim(b%name_visc_r))
      call check_str("defaults: name_visc_sd", trim(a%name_visc_sd), trim(b%name_visc_sd))
      call check_real("defaults: f_visc_sd", a%f_visc_sd, b%f_visc_sd)
      call check_real("defaults: f_visc_rel", a%f_visc_rel, b%f_visc_rel)
      call check_str("defaults: visc_shift_file", trim(a%visc_shift_file), trim(b%visc_shift_file))
      call check_real("defaults: visc_log10_min", a%visc_log10_min, b%visc_log10_min)
      call check_real("defaults: visc_log10_max", a%visc_log10_max, b%visc_log10_max)
      call check_real("defaults: visc_3d_lid_depth", a%visc_3d_lid_depth, b%visc_3d_lid_depth)
      call check_real("defaults: visc_3d_lid_log10max", a%visc_3d_lid_log10max, b%visc_3d_lid_log10max)
      call check_real("defaults: visc_3d_tol", a%visc_3d_tol, b%visc_3d_tol)
      call check_log("defaults: l_toroidal", a%l_toroidal, b%l_toroidal)
      call check_int("defaults: vilma_v1_jmax", a%vilma_v1_jmax, b%vilma_v1_jmax)
      call check_int("defaults: vilma_v1_l_prem", a%vilma_v1_l_prem, b%vilma_v1_l_prem)
      call check_int("defaults: vilma_v1_nsub", a%vilma_v1_nsub, b%vilma_v1_nsub)
   end subroutine same_defaults

   subroutine check_int(name, got, want)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: got, want
      if (got /= want) then
         write(*,'(a,a,a,i0,a,i0)') '   FAIL ', name, ': got ', got, ' want ', want
         ok = .false.
      end if
   end subroutine check_int

   subroutine check_real(name, got, want)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, want
      if (abs(got - want) > tol*max(abs(want), 1.0_wp)) then
         write(*,'(a,a,a,es12.4,a,es12.4)') '   FAIL ', name, ': got ', got, ' want ', want
         ok = .false.
      end if
   end subroutine check_real

   subroutine check_str(name, got, want)
      character(len=*), intent(in) :: name, got, want
      if (got /= want) then
         write(*,'(a,a,a,a,a,a)') '   FAIL ', name, ': got "', got, '" want "', want//'"'
         ok = .false.
      end if
   end subroutine check_str

   subroutine check_log(name, got, want)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: got, want
      if (got .neqv. want) then
         write(*,'(a,a,a,l1,a,l1)') '   FAIL ', name, ': got ', got, ' want ', want
         ok = .false.
      end if
   end subroutine check_log

end program test_params
