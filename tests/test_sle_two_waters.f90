program test_sle_two_waters
   !! The sea-level equation with two waters (vilma_sle header): melt water rho_water
   !! (fresh, the volume a melt adds to the ocean) and sea water rho_seawater (the
   !! flotation test, the ocean load), with the dilution term Δρ·C·D.
   !!
   !! Geometry: northern hemisphere land (+2000 m) carrying a grounded polar cap
   !! (10 degrees, 1000 m), the southern hemisphere ocean with a shelf (−100 m,
   !! colatitude 90-120 degrees) and a deep basin (−4000 m). One elastic solve each
   !! with two waters (1000, 1028) and with one (1000, 1000). Checks:
   !!   (1) volume: the ocean mean of C·rsl is the melt's fresh-water volume,
   !!       −(ρ_i/ρ_w)∫I dΩ / ∫C dΩ, in both (the subgrid fill is zero: no coast moves);
   !!   (2) mass: the converged load integrates to zero in both, to round-off;
   !!   (3) the dilution term matters: with two waters the load without it,
   !!       ρ_i I + ρ_sw C·rsl, integrates to (ρ_sw − ρ_w) times the melt volume;
   !!   (4) flotation is sea water's: an ice column with ρ_w d < ρ_i I < ρ_sw d floats
   !!       with two waters and grounds with one.
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, rho_ice
   use vilma_earth_structure, only: earth_model, build_M3L70V01
   use vilma_response,        only: response, response_init_elastic, response_destroy
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_synthesis, &
                                    sht_grid_surface_integral
   use vilma_sle,             only: sle_solve, sle_solver, sle_result, ocean_function
   use vilma_field,           only: spherical_cap
   implicit none

   integer,  parameter :: LMAX = 32
   real(wp), parameter :: RW = 1000.0_wp, RSW = 1028.0_wp, DEG = pi/180.0_wp
   type(sht_grid)    :: sht
   type(earth_model) :: em
   type(response)    :: resp
   real(wp), allocatable :: topo0(:,:), ice(:,:)
   real(wp) :: th, vol, mass, nodil
   real(wp) :: tf(1,1), hf(1,1), cf(1,1)
   integer  :: i, j
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   em = build_M3L70V01()
   call response_init_elastic(resp, em, LMAX)
   allocate(topo0(sht%nphi,sht%nlat), ice(sht%nphi,sht%nlat))
   do j = 1, sht%nlat
      th = sht%colat(j)/DEG
      do i = 1, sht%nphi
         topo0(i,j) = merge(2000.0_wp, merge(-100.0_wp, -4000.0_wp, th <= 120.0_wp), th <= 90.0_wp)
      end do
   end do
   call spherical_cap(sht, 0.0_wp, 0.0_wp, 10.0_wp*DEG, 1000.0_wp, ice)

   write(*,'(a)') ' SLE with two waters: polar cap on land, shelf and deep ocean'
   call solve_case(RSW, 'two waters (1000, 1028)', vol, mass, nodil)
   ! (3) with two waters
   if (abs(nodil - (RSW - RW)*vol)/abs((RSW - RW)*vol) > 1.0e-9_wp) then
      write(*,'(a,2es14.6)') ' FAIL (3): load without the dilution, expected', nodil, (RSW - RW)*vol;  ok = .false.
   end if
   call solve_case(RW, 'one water (1000, 1000)', vol, mass, nodil)

   ! (4) flotation: depth 4000 m, I with ρ_i I / d = 1014
   tf = -4000.0_wp;  hf = 4000.0_wp*1014.0_wp/rho_ice
   call ocean_function(tf, hf, cf, rho_ice, RSW)
   if (cf(1,1) /= 1.0_wp) then;  write(*,'(a)') ' FAIL (4): should float in sea water';  ok = .false.;  end if
   call ocean_function(tf, hf, cf, rho_ice, RW)
   if (cf(1,1) /= 0.0_wp) then;  write(*,'(a)') ' FAIL (4): should ground in fresh water';  ok = .false.;  end if

   call response_destroy(resp);  call sht_grid_destroy(sht)
   if (.not. ok) error stop 1
   write(*,'(a)') ' PASS: volume, mass, the dilution term and flotation with two waters'

contains

   subroutine solve_case(rho_sw, label, vol, mass, nodil)
      real(wp),         intent(in)  :: rho_sw
      character(len=*), intent(in)  :: label
      real(wp),         intent(out) :: vol, mass, nodil
      type(sle_solver) :: sle
      type(sle_result) :: res
      real(wp), allocatable :: rsl(:,:), C(:,:), load(:,:)
      complex(wp), allocatable :: sig(:)
      real(wp) :: cint, ref
      allocate(rsl(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat), load(sht%nphi,sht%nlat), sig(sht%nlm))
      sle%rho_ice = rho_ice;  sle%rho_water = RW;  sle%rho_seawater = rho_sw
      call sle_solve(sle, sht, resp, ice, ice, topo0, rsl, C, res, sigma_lm=sig)
      call sht_grid_synthesis(sht, sig, load)
      vol  = -(rho_ice/RW)*sht_grid_surface_integral(sht, ice)       ! the melt's fresh volume (a growth: < 0)
      cint = sht_grid_surface_integral(sht, C)
      ref  = rho_ice*sht_grid_surface_integral(sht, ice)
      mass = sht_grid_surface_integral(sht, load)
      nodil = sht_grid_surface_integral(sht, rho_ice*ice*(1.0_wp - C) + rho_sw*C*rsl)
      write(*,'(3a,f10.5,a,f10.5,a,es10.2)') '   ', label, ': ocean mean rsl', &
           sht_grid_surface_integral(sht, C*rsl)/cint, ' m, fresh volume', vol/cint, ' m, |mass|/ice', abs(mass)/ref
      if (abs(sht_grid_surface_integral(sht, C*rsl) - vol) > 1.0e-9_wp*abs(vol)) then
         write(*,'(a)') ' FAIL (1): the ocean mean is not the fresh-water volume';  ok = .false.
      end if
      ! degree 0 survives the analysis and synthesis exactly (Gauss quadrature)
      if (abs(mass) > 1.0e-9_wp*ref) then
         write(*,'(a)') ' FAIL (2): the load does not conserve mass';  ok = .false.
      end if
   end subroutine solve_case

end program test_sle_two_waters
