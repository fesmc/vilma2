program test_io_prec
   !! Output precision of vilma_io (the `prec` column of input/vilma-variables.md, ncio
   !! bit rounding). A small ve run with rotation, written three ways:
   !!
   !!   (1) diagnostics (vilma_write_step, vilma_write_horizontal): every rounded
   !!       variable within its table precision, max|err| <= max(prec_abs, prec_rel*|x|),
   !!       rounding actually happened, and the masks / scalars with prec "-" are exact;
   !!   (2) the lossless switches (the `lossless` argument; ncio's nc_set_lossless, as
   !!       &ctl write_lossless) write the model state bit for bit, no precision attributes;
   !!   (3) restarts are lossless whatever the table says, with the global switch off:
   !!       every variable of the file equals the model state exactly, the restore is
   !!       exact and the restarted run continues bit for bit.
   !!   (4) &ctl write_lossless is read; file sizes of the rounded and the lossless
   !!       diagnostics of the run are reported (rounded must be smaller).
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi
   use vilma_params,          only: vilma_param_class
   use vilma_control,         only: vilma_ctl_class, vilma_ctl_load
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_coupling,        only: solid_earth_finalize, solid_earth_update, solid_earth_init, solid_earth
   use vilma_io,              only: vilma_restart_write, vilma_restart_read, vilma_write_step, &
                                    vilma_write_horizontal
   use ncio,                  only: nc_read, nc_size, nc_set_lossless, nc_exists_attr
   use variable_io,           only: var_io_type, load_var_io_table, find_var_io_in_table
   implicit none

   integer, parameter :: LMAX = 16, NSTEP = 6, K2 = 3
   character(len=*), parameter :: F_RND = "libvilma/include/test_io_prec_rnd.nc"
   character(len=*), parameter :: F_LSL = "libvilma/include/test_io_prec_lsl.nc"
   character(len=*), parameter :: F_ARG = "libvilma/include/test_io_prec_arg.nc"
   character(len=*), parameter :: F_HRND = "libvilma/include/test_io_prec_hor_rnd.nc"
   character(len=*), parameter :: F_HLSL = "libvilma/include/test_io_prec_hor_lsl.nc"
   character(len=*), parameter :: F_RST = "libvilma/include/test_io_prec_rst.nc"
   character(len=*), parameter :: F_NML = "libvilma/include/test_io_prec.nml"
   ! diagnostics with a rounded table entry, and those that stay exact
   character(len=8), parameter :: ROUNDED(5) = [character(len=8) :: &
        "h_ice", "rsl", "z_bed", "dz_bed", "dz_ss"]
   character(len=11), parameter :: HOR(4) = [character(len=11) :: &
        "u_east", "u_north", "u_east_tor", "u_north_tor"]
   character(len=8), parameter :: DIAG(10) = [character(len=8) :: &
        "h_ice", "rsl", "z_bed", "dz_bed", "dz_ss", "C_ocean", "bsl", "bsl_ive", "dphi", "rot_m_re"]

   type(sht_grid), target :: sht
   type(var_io_type), allocatable :: tab(:)
   real(wp), allocatable  :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
   type(vilma_param_class) :: p
   type(solid_earth) :: a, b, c
   integer  :: np, nl, step, nt, q
   real(wp) :: dt_yr
   logical  :: ok
   real(wp), allocatable :: rnd(:,:), lsl(:,:), arg(:,:), ref(:,:)

   ok = .true.
   dt_yr = 1.0e3_wp
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   np = sht%nphi;  nl = sht%nlat
   allocate(z_bed_eq(np,nl), h_ice_eq(np,nl), h_ice(np,nl), rnd(np,nl), lsl(np,nl), arg(np,nl), ref(np,nl))
   call make_fields(z_bed_eq, h_ice_eq, h_ice)
   call load_var_io_table(tab, "input/vilma-variables.md")

   p%lmax = LMAX;  p%nlat = 2*LMAX;  p%nphi = 4*LMAX
   p%earth_response = "ve"
   p%rotation = .true.
   a%par = p; call solid_earth_init(a, z_bed_eq, h_ice_eq)

   ! ---- run, writing the diagnostics every step three ways ---------------------------
   do step = 0, NSTEP
      if (step > 0) call solid_earth_update(a, h_ice*(1.0_wp - 0.1_wp*real(step,wp)), dt_yr)
      call nc_set_lossless(.false.)
      call vilma_write_step(a, F_RND, a%time, nms=DIAG, init=(step == 0))
      call vilma_write_step(a, F_ARG, a%time, nms=DIAG, init=(step == 0), lossless=.true.)
      call vilma_write_horizontal(a, F_HRND, a%time, init=(step == 0))
      call nc_set_lossless(.true.)
      call vilma_write_step(a, F_LSL, a%time, nms=DIAG, init=(step == 0))
      call vilma_write_horizontal(a, F_HLSL, a%time, init=(step == 0))
      call nc_set_lossless(.false.)
   end do
   nt = nc_size(F_RND, "time")
   write(*,'(a,i0,a,i0,a)') ' io_prec: lmax=', LMAX, '  ', nt, ' output steps'
   if (nt /= NSTEP + 1) then;  write(*,'(a)') '   FAIL: wrong number of time steps';  ok = .false.;  end if

   ! ---- (1)+(2) 2-D diagnostics: rounded within prec; lossless exact -----------------
   do q = 1, size(ROUNDED)
      call model_field(a, trim(ROUNDED(q)), ref)
      call nc_read(F_RND, trim(ROUNDED(q)), rnd, start=[1,1,nt], count=[np,nl,1])
      call nc_read(F_LSL, trim(ROUNDED(q)), lsl, start=[1,1,nt], count=[np,nl,1])
      call nc_read(F_ARG, trim(ROUNDED(q)), arg, start=[1,1,nt], count=[np,nl,1])
      call check_rounded(trim(ROUNDED(q)), rnd, ref, ok)
      call check_exact(trim(ROUNDED(q))//" (nc_set_lossless)", lsl, ref, ok)
      call check_exact(trim(ROUNDED(q))//" (lossless arg)",    arg, ref, ok)
   end do
   call model_field(a, "C_ocean", ref)
   call nc_read(F_RND, "C_ocean", rnd, start=[1,1,nt], count=[np,nl,1])
   call check_exact("C_ocean (mask, rounded file)", rnd, ref, ok)

   ! ---- scalars: lossless in the table, so exact in the rounded file too -------------
   call check_scalar(F_RND, "bsl",      a%bsl,  nt, ok)
   call check_scalar(F_RND, "bsl_ive",  a%bsl_ive, nt, ok)
   call check_scalar(F_RND, "dphi",     a%dphi, nt, ok)
   call check_scalar(F_RND, "rot_m_re", real(a%rotation%m, wp), nt, ok)

   ! ---- horizontal displacement: rounded file vs the lossless file -------------------
   do q = 1, size(HOR)
      call nc_read(F_HRND, trim(HOR(q)), rnd, start=[1,1,nt], count=[np,nl,1])
      call nc_read(F_HLSL, trim(HOR(q)), lsl, start=[1,1,nt], count=[np,nl,1])
      ! the toroidal part is exactly zero without lateral viscosity: nothing to round
      call check_rounded(trim(HOR(q)), rnd, lsl, ok, must_change=(index(HOR(q), "tor") == 0))
   end do

   ! precision attributes: on the rounded variables of the rounded file only
   if (.not. nc_exists_attr(F_RND, "rsl", "precision_abs") .or. nc_exists_attr(F_LSL, "rsl", "precision_abs") &
       .or. nc_exists_attr(F_ARG, "rsl", "precision_abs")) then
      write(*,'(a)') '   FAIL: precision_abs attribute on rsl: expected in the rounded file only';  ok = .false.
   end if

   ! ---- (3) restart: lossless with the table's precision and the switch off ----------
   call vilma_restart_write(a, a%time, filename="test_io_prec_rst.nc", folder="libvilma/include", init=.true.)
   call check_restart_file(F_RST, a, ok)

   ! exact restore and bit-for-bit continuation
   b%par = p; call solid_earth_init(b, z_bed_eq, h_ice_eq)
   call vilma_restart_read(b, F_RST)
   call check_state("restore", b, a, ok)
   c%par = p; call solid_earth_init(c, z_bed_eq, h_ice_eq)
   call vilma_restart_read(c, F_RST)
   do step = 1, K2
      call solid_earth_update(a, h_ice, dt_yr)
      call solid_earth_update(c, h_ice, dt_yr)
   end do
   call check_state("continuation", c, a, ok)
   call solid_earth_finalize(a);  call solid_earth_finalize(b);  call solid_earth_finalize(c)

   ! ---- (4) &ctl write_lossless; sizes -------------------------------------------------
   call check_ctl(ok)
   call report_sizes(ok)

   call sht_grid_destroy(sht)
   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: diagnostics rounded within their table precision, lossless switches and'
      write(*,'(a)') '       restarts exact, rounded files smaller'
   else
      write(*,'(a)') ' FAIL: output precision did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine model_field(se, name, f)
      !! The model's field behind a diagnostic name.
      type(solid_earth), intent(in)  :: se
      character(len=*),  intent(in)  :: name
      real(wp),          intent(out) :: f(:,:)
      select case (name)
      case ("h_ice");   f = se%gg%h_ice
      case ("rsl");     f = se%gg%rsl
      case ("z_bed");   f = se%gg%z_bed
      case ("dz_bed");  f = se%gg%dz_bed
      case ("dz_ss");   f = se%gg%dz_ss
      case ("C_ocean"); f = se%gg%C
      case default;     error stop 'test_io_prec: unknown field'
      end select
   end subroutine model_field

   subroutine table_prec(name, pabs, prel)
      character(len=*), intent(in)  :: name
      real(wp),         intent(out) :: pabs, prel
      type(var_io_type) :: v
      call find_var_io_in_table(v, name, tab, with_error=.true.)
      pabs = v%prec_abs;  prel = v%prec_rel
   end subroutine table_prec

   subroutine check_rounded(name, x, xref, ok, must_change)
      !! x within the table precision of xref (max |err| <= max(prec_abs, prec_rel |xref|)),
      !! and x /= xref somewhere (the rounding is real).
      character(len=*), intent(in)    :: name
      real(wp),         intent(in)    :: x(:,:), xref(:,:)
      logical,          intent(inout) :: ok
      logical, optional, intent(in)   :: must_change
      real(wp) :: pabs, prel, tol(size(x,1),size(x,2)), emax, ratio
      call table_prec(name, pabs, prel)
      tol = max(pabs, prel*abs(xref))
      emax = maxval(abs(x - xref))
      ratio = 0.0_wp
      if (pabs > 0.0_wp .or. prel > 0.0_wp) ratio = maxval(abs(x - xref)/max(tol, tiny(1.0_wp)))
      write(*,'(a,a12,a,es9.2,a,es9.2,a,f6.3)') '   (1) ', name, ' prec_abs =', pabs, '  max|err| =', emax, &
           '  err/prec =', ratio
      if (pabs <= 0.0_wp .and. prel <= 0.0_wp) then
         write(*,'(a)') '   FAIL: '//name//' has no precision in the table';  ok = .false.
      else if (any(abs(x - xref) > tol)) then
         write(*,'(a)') '   FAIL: '//name//' exceeds its table precision';  ok = .false.
      else if (all(x == xref) .and. (.not. present(must_change) .or. must_change)) then
         write(*,'(a)') '   FAIL: '//name//' was not rounded';  ok = .false.
      end if
   end subroutine check_rounded

   subroutine check_exact(name, x, xref, ok)
      character(len=*), intent(in)    :: name
      real(wp),         intent(in)    :: x(:,:), xref(:,:)
      logical,          intent(inout) :: ok
      if (any(x /= xref)) then
         write(*,'(a,es9.2)') '   FAIL: '//name//' not exact, max|err| =', maxval(abs(x - xref));  ok = .false.
      end if
   end subroutine check_exact

   subroutine check_scalar(file, name, val, n, ok)
      character(len=*), intent(in)    :: file, name
      real(wp),         intent(in)    :: val
      integer,          intent(in)    :: n
      logical,          intent(inout) :: ok
      real(wp) :: x
      call nc_read(file, name, x, start=[n], count=[1])
      if (x /= val) then
         write(*,'(a,es9.2)') '   FAIL: scalar '//name//' not exact, err =', abs(x - val);  ok = .false.
      end if
   end subroutine check_scalar

   subroutine check_restart_file(file, se, ok)
      !! Every variable in the restart file equals the model state, bit for bit,
      !! and carries no precision attribute, though the table rounds the diagnostics.
      character(len=*),  intent(in)    :: file
      type(solid_earth), intent(in)    :: se
      logical,           intent(inout) :: ok
      integer :: i
      real(wp), allocatable :: x(:,:)
      allocate(x(np,nl))
      do i = 1, size(ROUNDED)
         call model_field(se, trim(ROUNDED(i)), ref)
         call nc_read(file, trim(ROUNDED(i)), x)
         call check_exact("restart "//trim(ROUNDED(i)), x, ref, ok)
         if (nc_exists_attr(file, trim(ROUNDED(i)), "precision_abs")) then
            write(*,'(a)') '   FAIL: restart '//trim(ROUNDED(i))//' carries a precision attribute';  ok = .false.
         end if
      end do
      call check_scalar(file, "bsl",  se%bsl,  1, ok)
      call check_scalar(file, "dphi", se%dphi, 1, ok)
      write(*,'(a)') '   (3) restart file: diagnostics exact, no precision attributes'
   end subroutine check_restart_file

   subroutine check_state(label, x, y, ok)
      !! The restored / continued state x equals y exactly: memory, diagnostics, clock.
      character(len=*),  intent(in)    :: label
      type(solid_earth), intent(in)    :: x, y
      logical,           intent(inout) :: ok
      real(wp) :: d
      d = max(maxval(abs(x%resp%Are - y%resp%Are)), maxval(abs(x%resp%Aim - y%resp%Aim)), &
              maxval(abs(x%resp%Bre - y%resp%Bre)), maxval(abs(x%resp%Bim - y%resp%Bim)), &
              maxval(abs(x%resp%Cre - y%resp%Cre)), maxval(abs(x%resp%Cim - y%resp%Cim)), &
              maxval(abs(x%gg%rsl - y%gg%rsl)), maxval(abs(x%gg%z_bed - y%gg%z_bed)), &
              maxval(abs(x%gg%dz_bed - y%gg%dz_bed)), maxval(abs(x%gg%dz_ss - y%gg%dz_ss)), &
              abs(x%dphi - y%dphi), abs(x%time - y%time), abs(x%rotation%m - y%rotation%m))
      write(*,'(a,a,a,es9.2)') '   (3) ', label, ' max|x - y| =', d
      if (d /= 0.0_wp) then
         write(*,'(a)') '   FAIL: '//label//' is not bit-identical';  ok = .false.
      end if
   end subroutine check_state

   subroutine check_ctl(ok)
      !! &ctl write_lossless: .false. in the shipped vilma.nml, read when set .true.
      logical, intent(inout) :: ok
      type(vilma_ctl_class) :: cd, cl
      character(len=512) :: line
      integer :: u_in, u_out, ios
      call vilma_ctl_load(cd, "vilma.nml")
      open(newunit=u_in,  file="vilma.nml", status="old", action="read")
      open(newunit=u_out, file=F_NML, status="replace", action="write")
      do
         read(u_in, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, "write_lossless") > 0 .and. index(adjustl(line), "!") /= 1) &
            line = "    write_lossless = .true."
         write(u_out, '(a)') trim(line)
      end do
      close(u_in);  close(u_out)
      call vilma_ctl_load(cl, F_NML)
      if (cd%write_lossless .or. .not. cl%write_lossless) then
         write(*,'(a)') '   FAIL: &ctl write_lossless (.false. in vilma.nml, read when .true.)';  ok = .false.
      end if
   end subroutine check_ctl

   subroutine report_sizes(ok)
      !! File sizes of the rounded and the lossless diagnostics of the run.
      logical, intent(inout) :: ok
      integer(8) :: s_rnd, s_lsl, h_rnd, h_lsl
      inquire(file=F_RND,  size=s_rnd);  inquire(file=F_LSL,  size=s_lsl)
      inquire(file=F_HRND, size=h_rnd);  inquire(file=F_HLSL, size=h_lsl)
      write(*,'(a,i0,a,i0,a,f5.1,a)') '   (4) size, file_out:  rounded ', s_rnd, '  lossless ', s_lsl, &
           '  -> ', 100.0_wp*real(s_rnd,wp)/real(s_lsl,wp), ' %'
      write(*,'(a,i0,a,i0,a,f5.1,a)') '       size, file_hor:  rounded ', h_rnd, '  lossless ', h_lsl, &
           '  -> ', 100.0_wp*real(h_rnd,wp)/real(h_lsl,wp), ' %'
      if (s_rnd >= s_lsl .or. h_rnd >= h_lsl) then
         write(*,'(a)') '   FAIL: rounded output is not smaller than lossless';  ok = .false.
      end if
   end subroutine report_sizes

   subroutine make_fields(z_bed_eq, h_ice_eq, h_ice)
      !! Polar land cap (colat<50 deg, +500 m) over deep ocean (-4000 m); no
      !! reference ice; a 2 km grounded ice load on colat<30 deg, west half only.
      real(wp), intent(out) :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
            ! a half cap in longitude, so the horizontal displacement has an east part
            h_ice(i,j)    = merge(2000.0_wp, 0.0_wp,     thd < 30.0_wp .and. sht%lon(i) < pi)
         end do
      end do
      h_ice_eq = 0.0_wp
   end subroutine make_fields

end program test_io_prec
