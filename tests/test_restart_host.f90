program test_restart_host
   !! Restart with host-resolved terms (vilma_io): the reference terms are static and
   !! checked, the terms of the last update (the next ramp start) are restored. lmax 12,
   !! forward Euler, rotation on.
   !!
   !!   (1) state restore: the last snapshot read into a fresh hosted model gives the same
   !!       rsl and the same last terms, bit for bit;
   !!   (2) continuation: restoring the earlier snapshot and stepping on with the same ice
   !!       and terms reproduces the uninterrupted run; with the last terms reset to the
   !!       reference instead it does not (the ramp start is part of the state);
   !!   (3) mismatches stop: a hosted file into a model without host terms, a plain file
   !!       into a hosted model, and a hosted file into a model with another reference.
   !!       Each runs as a child process of this program and must exit with an error.
   use vilma_precision, only: wp
   use vilma_constants, only: pi
   use vilma_params,    only: vilma_param_class
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_sht,       only: sht_grid, sht_grid_init
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   use vilma_io,        only: vilma_restart_write, vilma_restart_read
   implicit none

   integer,  parameter :: LMAX = 12, K1 = 3, K2 = 3
   real(wp), parameter :: DEG = pi/180.0_wp
   character(len=*), parameter :: FHOST  = "libvilma/include/test_restart_host.nc"
   character(len=*), parameter :: FPLAIN = "libvilma/include/test_restart_host_plain.nc"
   type(sht_grid), target :: sht
   integer  :: np, nl
   real(wp), allocatable :: z_eq(:,:), h_eq(:,:), h1(:,:), wr(:,:), Ir(:,:), Cr(:,:)
   real(wp), allocatable :: w0(:,:), I0(:,:), C0(:,:), F0(:,:), w1(:,:), I1(:,:), C1(:,:), F1(:,:)
   character(len=256) :: mode, exe
   logical :: ok

   call sht_grid_init(sht, LMAX, nlat=2*LMAX+2, nphi=4*LMAX)
   np = sht%nphi;  nl = sht%nlat                  ! the grid may round nlat up
   allocate(z_eq(np,nl), h_eq(np,nl), h1(np,nl), wr(np,nl), Ir(np,nl), Cr(np,nl), &
            w0(np,nl), I0(np,nl), C0(np,nl), F0(np,nl), w1(np,nl), I1(np,nl), C1(np,nl), F1(np,nl))
   call scene()

   call get_command_argument(1, mode)
   if (len_trim(mode) > 0) then
      call child(trim(mode))                 ! a mismatch case, run as a child: must stop
      stop 0
   end if

   ok = .true.
   call roundtrip(ok)
   call get_command_argument(0, exe)
   call mismatches(trim(exe), ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: host terms round-trip through the restart, mismatches stop'
      call radial_fe_finalize()
   else
      write(*,'(a)') ' FAIL: restart with host terms'
      call radial_fe_finalize()
      error stop 1
   end if

contains

   subroutine scene()
      !! Land cap, shelf, deep ocean; ice-free reference; a host region (w = 1 / 1/2 by
      !! longitude, colat 20..70) whose coast and ice change between t0 and t1.
      integer  :: i, j
      real(wp) :: th, lo, ww, hh, fg, fr
      do j = 1, nl
         th = sht%colat(j)/DEG
         do i = 1, np
            lo = sht%lon(i)/DEG
            z_eq(i,j) = merge(500.0_wp, merge(-100.0_wp, -4000.0_wp, th < 62.0_wp), th < 50.0_wp)
            h1(i,j) = 0.0_wp
            if (th < 30.0_wp) h1(i,j) = 2000.0_wp
            if (th >= 50.0_wp .and. th < 62.0_wp) h1(i,j) = 300.0_wp
            ww = 0.0_wp
            if (th >= 20.0_wp .and. th < 70.0_wp .and. lo < 180.0_wp) ww = merge(1.0_wp, 0.5_wp, lo < 90.0_wp)
            hh = merge(3000.0_wp, 200.0_wp, th < 30.0_wp)
            fg = merge(1.0_wp, 0.0_wp, th < 54.0_wp)
            fr = merge(1.0_wp, 0.0_wp, th < 58.0_wp)
            wr(i,j) = ww;  Ir(i,j) = 0.0_wp;  Cr(i,j) = ww*(1.0_wp - fr)
            w1(i,j) = ww;  I1(i,j) = ww*hh*fg;  C1(i,j) = ww*(1.0_wp - fg)
            F1(i,j) = ww*merge(200.0_wp, -50.0_wp, fg > 0.5_wp)*((1.0_wp - fg) - (1.0_wp - fr))
            w0(i,j) = 0.5_wp*ww;  I0(i,j) = 0.5_wp*ww*0.25_wp*hh*fr;  C0(i,j) = 0.5_wp*ww*(1.0_wp - fr);  F0(i,j) = 0.0_wp
         end do
      end do
      h_eq = 0.0_wp
   end subroutine scene

   subroutine new_par(p)
      type(vilma_param_class), intent(out) :: p
      p%lmax = LMAX;  p%nlat = nl;  p%nphi = np
      p%rotation = .true.
   end subroutine new_par

   subroutine step(se, k)
      !! Update number k (1..K1+K2): the ice and the terms move from t0 values to t1 values.
      type(solid_earth), intent(inout) :: se
      integer,           intent(in)    :: k
      real(wp) :: f
      f = min(real(k,wp)/real(K1+K2,wp), 1.0_wp)
      call solid_earth_update(se, f*h1, 500.0_wp, w0 + f*(w1 - w0), I0 + f*(I1 - I0), C0 + f*(C1 - C0), F0 + f*(F1 - F0))
   end subroutine step

   subroutine roundtrip(ok)
      logical, intent(inout) :: ok
      type(vilma_param_class) :: p
      type(solid_earth) :: a, b, c, d
      integer  :: k
      real(wp) :: t1, d_res, d_last, d_cont, d_ctl
      call new_par(p)
      a%par = p;  b%par = p;  c%par = p;  d%par = p
      call solid_earth_init(a, z_eq, h_eq, w=wr, I_g=Ir, C=Cr)
      call solid_earth_update(a, h_eq, 0.0_wp, w0, I0, C0, F0)           ! seed: the ramp start is w0..
      do k = 1, K1
         call step(a, k)
      end do
      t1 = a%time
      call vilma_restart_write(a, t1, filename=FHOST, init=.true.)
      do k = K1+1, K1+K2
         call step(a, k)
      end do
      call vilma_restart_write(a, a%time, filename=FHOST, init=.false.)

      ! (1) the last snapshot into a fresh model
      call solid_earth_init(b, z_eq, h_eq, w=wr, I_g=Ir, C=Cr)
      call vilma_restart_read(b, FHOST)
      d_res  = max(maxval(abs(b%rsl - a%rsl)), maxval(abs(b%gg%C - a%gg%C)), abs(b%bsl - a%bsl), abs(b%bsl_ive - a%bsl_ive), &
                   abs(b%bsl_ive_host - a%bsl_ive_host))
      d_last = max(maxval(abs(b%gg%host%last%w - a%gg%host%last%w)), maxval(abs(b%gg%host%last%I_g - a%gg%host%last%I_g)), &
                   maxval(abs(b%gg%host%last%C - a%gg%host%last%C)), maxval(abs(b%gg%host%last%F - a%gg%host%last%F)))
      write(*,'(a)') ' (1) state restore'
      write(*,'(a,es10.2,a,es10.2)') '     max|rsl, C, bsl| =', d_res, '   max|last terms| =', d_last
      if (d_res /= 0.0_wp .or. d_last /= 0.0_wp .or. maxval(abs(a%gg%host%last%F)) < 1.0_wp) then
         write(*,'(a)') '     FAIL: restored state or last terms differ from the written ones';  ok = .false.
      end if

      ! (2) continuation from the earlier snapshot
      call solid_earth_init(c, z_eq, h_eq, w=wr, I_g=Ir, C=Cr)
      call vilma_restart_read(c, FHOST, time=t1)
      do k = K1+1, K1+K2
         call step(c, k)
      end do
      d_cont = max(maxval(abs(c%rsl - a%rsl)), maxval(abs(c%gg%dz_bed - a%gg%dz_bed)), abs(c%dphi - a%dphi), &
                   maxval(abs(c%gg%dIg - a%gg%dIg)))
      ! control: the same, but the ramp start reset to the reference
      call solid_earth_init(d, z_eq, h_eq, w=wr, I_g=Ir, C=Cr)
      call vilma_restart_read(d, FHOST, time=t1)
      d%gg%host%last%w = d%gg%host%ref%w;  d%gg%host%last%I_g = d%gg%host%ref%I_g
      d%gg%host%last%C = d%gg%host%ref%C;  d%gg%host%last%F = 0.0_wp
      do k = K1+1, K1+K2
         call step(d, k)
      end do
      d_ctl = maxval(abs(d%rsl - a%rsl))
      write(*,'(a)') ' (2) continuation from the earlier snapshot'
      write(*,'(a,es10.2,a,es10.2,a)') '     max|continued - uninterrupted| =', d_cont, &
           '   with the ramp start reset: ', d_ctl, ' m'
      if (d_cont > 1.0e-9_wp) then
         write(*,'(a)') '     FAIL: the restarted run does not continue the uninterrupted one';  ok = .false.
      end if
      if (d_ctl < 1.0e-3_wp) then
         write(*,'(a)') '     FAIL: the ramp start does not matter here, the check is blind';  ok = .false.
      end if
      call solid_earth_finalize(a);  call solid_earth_finalize(b)
      call solid_earth_finalize(c);  call solid_earth_finalize(d)

      ! a plain file for the mismatch cases
      call solid_earth_init(a, z_eq, h_eq)
      call solid_earth_update(a, h1, 500.0_wp)
      call vilma_restart_write(a, a%time, filename=FPLAIN, init=.true.)
      call solid_earth_finalize(a)
   end subroutine roundtrip

   subroutine mismatches(exe, ok)
      character(len=*), intent(in)    :: exe
      logical,          intent(inout) :: ok
      character(len=*), parameter :: M(3) = [character(len=24) :: 'host-into-plain', 'plain-into-host', 'other-reference']
      integer :: q, ex
      write(*,'(a)') ' (3) mismatches stop (child processes)'
      do q = 1, size(M)
         call execute_command_line(exe//' '//trim(M(q))//' > /dev/null 2>&1', exitstat=ex)
         write(*,'(a,a,a,i0)') '     ', trim(M(q)), ': exit status ', ex
         if (ex == 0) then
            write(*,'(a)') '     FAIL: the read did not stop';  ok = .false.
         end if
      end do
   end subroutine mismatches

   subroutine child(m)
      !! One mismatching read; reaching its end is the failure (exit status 0).
      character(len=*), intent(in) :: m
      type(vilma_param_class) :: p
      type(solid_earth) :: se
      call new_par(p)
      se%par = p
      select case (m)
      case ('host-into-plain')
         call solid_earth_init(se, z_eq, h_eq)
         call vilma_restart_read(se, FHOST)
      case ('plain-into-host')
         call solid_earth_init(se, z_eq, h_eq, w=wr, I_g=Ir, C=Cr)
         call vilma_restart_read(se, FPLAIN)
      case ('other-reference')
         call solid_earth_init(se, z_eq, h_eq, w=0.5_wp*wr, I_g=Ir, C=0.5_wp*Cr)
         call vilma_restart_read(se, FHOST)
      end select
      write(*,'(a)') ' NOT STOPPED'
   end subroutine child

end program test_restart_host
