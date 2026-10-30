program diag_restart_3d
   !! Diagnostic (not in `make check`): restart exactness with 3-D viscosity at production
   !! resolution. Runs the uninterrupted model K1+K2 updates and a fresh model restored after K1,
   !! and compares them bit for bit after every update. test_restart_3d is the small case in
   !! `make check`; this one takes the resolution, the viscosity file and the options from the
   !! command line, to run at lmax 170 on the FesmData Bagge2021 field:
   !!
   !!   diag_restart_3d.x lmax visc_file spin host k1 k2 dt_yr [scheme]
   !!
   !! lmax and visc_file: the grid is (2 lmax + 2) x (4 lmax), the field a FesmData product (structure 6);
   !! spin = 1: 1-D pre-spin-up of the uninterrupted model (the restored one never spins up);
   !! host = 1: host-resolved terms (the scene of test_restart_host); k1, k2: updates before and
   !! after the restore, of dt_yr years each; scheme: fe (default) | trap. Rotation is on.
   use vilma_precision, only: wp
   use vilma_constants, only: pi
   use vilma_params,    only: vilma_param_class
   use vilma_radial_fe, only: radial_fe_finalize
   use vilma_sht,       only: sht_grid, sht_grid_init
   use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_spinup, &
                              solid_earth_finalize
   use vilma_io,        only: vilma_restart_write, vilma_restart_read
   implicit none

   real(wp), parameter :: DEG = pi/180.0_wp
   character(len=512) :: arg, visc, file, scheme
   integer  :: lmax, np, nl, k1, k2, k, spin, host, nbad
   real(wp) :: dt, t1, f, dtry
   type(sht_grid), target :: sht
   type(vilma_param_class) :: p
   type(solid_earth) :: a, c
   real(wp), allocatable :: z_eq(:,:), h_eq(:,:), h1(:,:), wr(:,:), Ir(:,:), Cr(:,:)
   real(wp), allocatable :: w0(:,:), I0(:,:), C0(:,:), F0(:,:), w1(:,:), I1(:,:), C1(:,:), F1(:,:)
   real(wp), allocatable :: rsl(:,:,:), dzb(:,:,:), sc(:,:)
   logical :: ex, fixdt

   call get_command_argument(1, arg);  read(arg,*) lmax
   call get_command_argument(2, visc)
   call get_command_argument(3, arg);  read(arg,*) spin
   call get_command_argument(4, arg);  read(arg,*) host
   call get_command_argument(5, arg);  read(arg,*) k1
   call get_command_argument(6, arg);  read(arg,*) k2
   call get_command_argument(7, arg);  read(arg,*) dt
   scheme = "fe";  call get_command_argument(8, arg);  if (len_trim(arg) > 0) scheme = arg
   call get_command_argument(9, arg);  fixdt = (arg(1:1) == 'x')
   file = "diag_restart_3d.nc"

   call sht_grid_init(sht, lmax, nlat=2*lmax+2, nphi=4*lmax)
   np = sht%nphi;  nl = sht%nlat
   allocate(z_eq(np,nl), h_eq(np,nl), h1(np,nl), wr(np,nl), Ir(np,nl), Cr(np,nl), &
            w0(np,nl), I0(np,nl), C0(np,nl), F0(np,nl), w1(np,nl), I1(np,nl), C1(np,nl), F1(np,nl))
   call scene()
   allocate(rsl(np,nl,k1+k2), dzb(np,nl,k1+k2), sc(4,k1+k2))

   p%lmax = lmax;  p%nlat = nl;  p%nphi = np
   p%rotation = .true.;  p%scheme = trim(scheme)
   p%earth_response = "ve";  p%l_visc_3d = .true.;  p%visc_3d_file = trim(visc);  p%visc_3d_structure = 6
   p%pre_spinup_1d = (spin == 1);  p%equil_time_max = 0.0_wp

   a%par = p
   call init_model(a)
   if (spin == 1) call spin_model(a)
   do k = 1, k1+k2
      call step(a, k)
      rsl(:,:,k) = a%rsl;  dzb(:,:,k) = a%dz_bed;  sc(:,k) = [a%bsl, a%dphi, a%fill_int, a%worst_mass_resid]
      write(*,'(a,i0,a,es10.3,a,es10.3)') ' A step ', k, ': bsl ', a%bsl, '  max|rsl| ', maxval(abs(a%rsl))
      if (k == k1) then
         t1 = a%time;  dtry = a%stepper%dt_try
         call vilma_restart_write(a, t1, filename=trim(file), init=.true.)
      end if
   end do

   c%par = p
   call init_model(c)
   call vilma_restart_read(c, trim(file), time=t1)
   write(*,'(a,2es24.16,l2)') ' dt_try A, C:', dtry, c%stepper%dt_try, dtry == c%stepper%dt_try
   if (fixdt) c%stepper%dt_try = dtry
   nbad = 0
   do k = k1+1, k1+k2
      call step(c, k)
      write(*,'(a,i0,a,es10.2,a,es10.2,a,es10.2,a,es10.2)') ' C step ', k, ': max|rsl - A| ', &
           maxval(abs(c%rsl - rsl(:,:,k))), '  max|dz_bed - A| ', maxval(abs(c%dz_bed - dzb(:,:,k))), &
           '  |bsl - A| ', abs(c%bsl - sc(1,k)), '  |fill_int - A| ', abs(c%fill_int - sc(3,k))
      if (.not. (all(c%rsl == rsl(:,:,k)) .and. all(c%dz_bed == dzb(:,:,k)) .and. all([c%bsl, c%dphi, c%fill_int, c%worst_mass_resid] == sc(:,k)))) nbad = nbad + 1
   end do
   write(*,'(a,i0,a,i0)') ' updates differing after the restore: ', nbad, ' of ', k2
   inquire(file=trim(file), exist=ex)
   call radial_fe_finalize()
   if (nbad /= 0) error stop 'diag_restart_3d: the restarted run is not bit-identical'
   write(*,'(a)') ' BIT-IDENTICAL'

contains

   subroutine init_model(se)
      type(solid_earth), intent(inout) :: se
      if (host == 1) then
         call solid_earth_init(se, z_eq, h_eq, w=wr, I_g=Ir, C=Cr)
      else
         call solid_earth_init(se, z_eq, h_eq)
      end if
   end subroutine init_model

   subroutine spin_model(se)
      type(solid_earth), intent(inout) :: se
      if (host == 1) then
         call solid_earth_spinup(se, 0.5_wp*h1, w0, I0, C0, F0)
      else
         call solid_earth_spinup(se, 0.5_wp*h1)
      end if
   end subroutine spin_model

   subroutine step(se, k)
      type(solid_earth), intent(inout) :: se
      integer,           intent(in)    :: k
      real(wp) :: f
      f = min(0.5_wp + real(k,wp)/real(2*(k1+k2),wp), 1.0_wp)
      if (host == 1) then
         call solid_earth_update(se, f*h1, dt, w0 + f*(w1 - w0), I0 + f*(I1 - I0), C0 + f*(C1 - C0), F0 + f*(F1 - F0))
      else
         call solid_earth_update(se, f*h1, dt)
      end if
   end subroutine step

   subroutine scene()
      !! test_restart_host's scene: land cap, shelf, deep ocean; ice-free reference; a host region
      !! (w = 1 / 1/2 by longitude, colat 20..70) whose coast and ice change between t0 and t1.
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
            if (host == 1 .and. th >= 20.0_wp .and. th < 70.0_wp .and. lo < 180.0_wp) ww = merge(1.0_wp, 0.5_wp, lo < 90.0_wp)
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

end program diag_restart_3d
