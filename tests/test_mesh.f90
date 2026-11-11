program test_mesh
   !! Validate the radial mesh (radial_mesh_build). Every mesh must span
   !! [0, r_earth] exactly (Martinec 2000 meshes ⟨0,a⟩ through the centre), be
   !! strictly ascending, place a node on every layer interface — including the
   !! CMB (no element straddles a discontinuity) — and keep every element within
   !! the target size for ITS OWN depth (default 5/10/40 km above 70 km / to
   !! 670 km / below), however thick the layer that holds it.
   !!
   !!   (1) M3-L70-V01: 217 elements, uniform within each layer (the layers lie
   !!       within one depth zone each).
   !!   (2) the 11-layer PREM table: 219 elements; its 80 km lithosphere crosses
   !!       the 70 km zone boundary and gets 14 x 5 km + 1 x 10 km.
   !!   (3) one layer from the surface to the CMB: resolved by depth like the
   !!       layered models, 10 km elements across 100-200 km depth (a rule that
   !!       sizes by the layer's mid-depth gives 40 km there).
   !!   (4) a layer interface 1 km below a zone boundary leaves no sliver element.
   !!   (5) the mesh table carried by the earth model (namelist mesh_dr,
   !!       mesh_depth, mesh_dr_scale): half the sizes double the element count,
   !!       and a one-zone table gives a uniform mesh.
   use vilma_precision,       only: wp
   use vilma_earth_structure, only: earth_model, earth_layer, build_M3L70V01, &
                                    RHEOL_MAXWELL, RHEOL_FLUID
   use vilma_radial_fe,       only: radial_mesh_build, radial_mesh
   implicit none

   real(wp), parameter :: km = 1.0e3_wp
   real(wp), parameter :: tolnode = 1.0_wp     ! 1 m
   real(wp), parameter :: A = 6371.0_wp*km
   real(wp), parameter :: DR_DEF(3) = [5.0_wp, 10.0_wp, 40.0_wp]*km
   real(wp), parameter :: DEPTH_DEF(2) = [70.0_wp, 670.0_wp]*km
   ! The PREM table of input/vilma_defaults.nml (layer tops, surface first).
   real(wp), parameter :: PREM_TOP(11) = [6371.0_wp, 6291.0_wp, 6151.0_wp, 5971.0_wp, &
        5771.0_wp, 5701.0_wp, 5600.0_wp, 4943.0_wp, 4287.0_wp, 3630.0_wp, 3480.0_wp]*km

   type(earth_model) :: em
   type(radial_mesh) :: mesh
   logical  :: ok
   integer  :: i, ne_m3
   real(wp) :: dr, depth

   ok = .true.

   ! --- (1) M3-L70-V01 ----------------------------------------------------------
   em = build_M3L70V01()
   call radial_mesh_build(mesh, em)
   call check_mesh('M3-L70-V01', em, mesh, DR_DEF, DEPTH_DEF)
   ne_m3 = mesh%ne
   call expect_int('M3-L70-V01 ne', mesh%ne, 217)
   do i = 1, mesh%ne
      if (abs((mesh%r(i+1) - mesh%r(i)) - layer_uniform_dr(em, mesh, i)) > 1.0e-6_wp) then
         print '(a,i0)', ' FAIL: M3-L70-V01 is not uniform within its layers at element ', i
         ok = .false.;  exit
      end if
   end do

   ! --- (2) PREM, 11 layers -----------------------------------------------------
   em = stack(PREM_TOP)
   call radial_mesh_build(mesh, em)
   call check_mesh('PREM-11', em, mesh, DR_DEF, DEPTH_DEF)
   call expect_int('PREM-11 ne', mesh%ne, 219)
   ! lithosphere: 14 elements of 5 km above 70 km depth, one of 10 km below
   call expect_int('PREM-11 lithosphere elements', count(mesh%elem_layer == 1), 15)
   dr = mesh%r(mesh%nr) - mesh%r(mesh%nr-1)
   if (abs(dr - 5.0_wp*km) > tolnode) then
      print '(a,f8.3,a)', ' FAIL: PREM-11 surface element is ', dr/km, ' km, not 5'
      ok = .false.
   end if
   dr = mesh%r(mesh%nr-14) - mesh%r(mesh%nr-15)
   if (abs(dr - 10.0_wp*km) > tolnode) then
      print '(a,f8.3,a)', ' FAIL: PREM-11 element at 70-80 km depth is ', dr/km, ' km, not 10'
      ok = .false.
   end if

   ! --- (3) one mantle layer, surface to CMB ------------------------------------
   em = stack([A, 3480.0_wp*km])
   call radial_mesh_build(mesh, em)
   call check_mesh('one-layer mantle', em, mesh, DR_DEF, DEPTH_DEF)
   ! 14 + 60 + 2221/40 = 129.5 -> 130 mantle elements, as a 5/10/40 km stack has
   call expect_int('one-layer mantle elements', count(mesh%elem_layer == 1), 130)
   do i = 1, mesh%ne
      depth = A - 0.5_wp*(mesh%r(i) + mesh%r(i+1))
      if (depth > 100.0_wp*km .and. depth < 200.0_wp*km .and. &
          mesh%r(i+1) - mesh%r(i) > 10.0_wp*km + tolnode) then
         print '(a,f7.2,a,f7.1,a)', ' FAIL: one-layer mantle: element of ', &
              (mesh%r(i+1) - mesh%r(i))/km, ' km at ', depth/km, ' km depth (target 10 km)'
         ok = .false.
      end if
   end do

   ! --- (4) no sliver where an interface sits just below a zone boundary --------
   em = stack([A, A - 71.0_wp*km, 3480.0_wp*km])
   call radial_mesh_build(mesh, em)
   call check_mesh('71 km lithosphere', em, mesh, DR_DEF, DEPTH_DEF)
   call expect_int('71 km lithosphere elements', count(mesh%elem_layer == 1), 15)
   do i = 1, mesh%ne
      if (mesh%elem_layer(i) == 1 .and. mesh%r(i+1) - mesh%r(i) < 4.0_wp*km) then
         print '(a,f7.3,a)', ' FAIL: sliver element of ', (mesh%r(i+1) - mesh%r(i))/km, ' km'
         ok = .false.
      end if
   end do

   ! --- (5) the mesh table of the earth model -----------------------------------
   em = build_M3L70V01()
   em%mesh_dr = 0.5_wp*DR_DEF;  em%mesh_depth = DEPTH_DEF
   call radial_mesh_build(mesh, em)
   call check_mesh('M3-L70-V01, half size', em, mesh, em%mesh_dr, em%mesh_depth)
   call expect_int('M3-L70-V01 ne at half size', mesh%ne, 2*ne_m3)

   em%mesh_dr = [20.0_wp*km];  em%mesh_depth = [real(wp) ::]
   call radial_mesh_build(mesh, em)
   call check_mesh('M3-L70-V01, one 20 km zone', em, mesh, em%mesh_dr, em%mesh_depth)
   ! ceiling of 70, 350, 250, 2221, 3480 km over 20 km
   call expect_int('M3-L70-V01 ne, one 20 km zone', mesh%ne, 4 + 18 + 13 + 112 + 174)

   if (ok) then
      print '(a)', ' PASS: radial mesh spans the sphere, captures interfaces, sizes elements by their depth'
   else
      print '(a)', ' FAIL: radial mesh consistency check'
      error stop 1
   end if

contains

   function stack(r_top) result(e)
      !! A layered model with the given layer tops (surface first), the last
      !! layer reaching the centre. Only the geometry matters to the mesh.
      real(wp), intent(in) :: r_top(:)
      type(earth_model) :: e
      integer  :: k, n
      real(wp) :: r_bot
      n = size(r_top)
      e%name = "stack";  e%r_earth = r_top(1);  e%r_core = r_top(n)
      allocate(e%layers(n))
      do k = 1, n
         r_bot = 0.0_wp
         if (k < n) r_bot = r_top(k+1)
         e%layers(k) = earth_layer(r_bot, r_top(k), 4000.0_wp, 1.0e11_wp, 1.0e21_wp, &
                                   merge(RHEOL_FLUID, RHEOL_MAXWELL, k == n))
         if (k == n) e%layers(k)%mu = 0.0_wp
      end do
   end function stack

   pure function target_at(depth, dr_zone, depth_zone) result(h)
      !! Target element size at a depth, from the zone table.
      real(wp), intent(in) :: depth, dr_zone(:), depth_zone(:)
      real(wp) :: h
      integer :: k
      h = dr_zone(size(dr_zone))
      do k = size(depth_zone), 1, -1
         if (depth <= depth_zone(k)) h = dr_zone(k)
      end do
   end function target_at

   pure function layer_uniform_dr(e, m, i) result(h)
      !! Element size if the layer of element i were divided uniformly.
      type(earth_model), intent(in) :: e
      type(radial_mesh), intent(in) :: m
      integer,           intent(in) :: i
      real(wp) :: h
      h = (e%layers(m%elem_layer(i))%r_top - e%layers(m%elem_layer(i))%r_bot) &
          / real(count(m%elem_layer == m%elem_layer(i)), wp)
   end function layer_uniform_dr

   subroutine check_mesh(name, e, m, dr_zone, depth_zone)
      !! The properties every mesh has.
      character(*),      intent(in) :: name
      type(earth_model), intent(in) :: e
      type(radial_mesh), intent(in) :: m
      real(wp),          intent(in) :: dr_zone(:), depth_zone(:)
      integer  :: i, k
      logical  :: found, good
      real(wp) :: dr, rmid, target, hit

      good = .true.
      ! node count, bounds, strictly ascending
      if (m%nr /= m%ne + 1 .or. size(m%r) /= m%nr .or. size(m%elem_layer) /= m%ne) good = .false.
      if (m%r(1) /= 0.0_wp) good = .false.
      if (abs(m%r(m%nr) - e%r_earth) > tolnode) good = .false.
      do i = 2, m%nr
         if (m%r(i) <= m%r(i-1)) good = .false.
      end do
      ! every layer interface is a node
      do k = 2, size(e%layers)
         found = .false.
         do i = 1, m%nr
            if (abs(m%r(i) - e%layers(k)%r_top) <= tolnode) found = .true.
         end do
         if (.not. found) then
            print '(3a,f8.1,a)', ' ', name, ': MISSING interface node at ', e%layers(k)%r_top/km, ' km'
            good = .false.
         end if
      end do
      ! element sizes: within the target at the element's shallow or deep end
      ! (an element across a zone boundary lies between the two targets), and
      ! the layer assignment contains the midpoint
      hit = 0.0_wp
      do i = 1, m%ne
         dr     = m%r(i+1) - m%r(i)
         rmid   = 0.5_wp*(m%r(i+1) + m%r(i))
         target = max(target_at(e%r_earth - m%r(i+1), dr_zone, depth_zone), &
                      target_at(e%r_earth - m%r(i),   dr_zone, depth_zone))
         if (dr > target + tolnode) then
            print '(3a,i0,a,f7.2,a,f7.2,a)', ' ', name, ': element ', i, ' size ', dr/km, &
                 ' km exceeds target ', target/km, ' km'
            good = .false.
         end if
         if (rmid < e%layers(m%elem_layer(i))%r_bot .or. &
             rmid > e%layers(m%elem_layer(i))%r_top) good = .false.
         hit = max(hit, dr)
      end do
      print '(1x,a,a,i0,a,f7.2,a)', name, ': ', m%ne, ' elements, largest ', hit/km, ' km'
      if (.not. good) then
         print '(3a)', ' FAIL: ', name, ': mesh consistency'
         ok = .false.
      end if
   end subroutine check_mesh

   subroutine expect_int(what, got, want)
      character(*), intent(in) :: what
      integer,      intent(in) :: got, want
      if (got /= want) then
         print '(3a,i0,a,i0)', ' FAIL: ', what, ' = ', got, ', expected ', want
         ok = .false.
      end if
   end subroutine expect_int

end program test_mesh
