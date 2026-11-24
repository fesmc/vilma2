module vilma_constants
   !! Physical and reference constants for VILMA2.
   !!
   !! Values follow the GIA-community conventions used by the Spada et al. (2011)
   !! and Martinec et al. (2018) benchmarks, so that validation runs can be
   !! compared directly. The ice density is the exception: its default is the
   !! standard 917; the benchmarks' 931 is set by the benchmark tests themselves.
   !! Densities must be applied consistently (see doc/design.md, "conserve mass
   !! not volume").
   use vilma_precision, only: wp
   implicit none
   public

   real(wp), parameter :: pi      = acos(-1.0_wp)
   real(wp), parameter :: deg2rad = pi/180.0_wp
   real(wp), parameter :: rad2deg = 180.0_wp/pi

   ! --- Gravitation -----------------------------------------------------------
   real(wp), parameter :: grav_G = 6.67430e-11_wp   !! Newtonian G [m^3 kg^-1 s^-2]

   ! --- Reference Earth (PREM-consistent, as in the GIA benchmarks) -----------
   real(wp), parameter :: r_earth = 6.371e6_wp      !! mean radius [m]
   real(wp), parameter :: m_earth = 5.972e24_wp     !! mass [kg]
   real(wp), parameter :: g_surf  = 9.81_wp         !! surface gravity [m s^-2]
   real(wp), parameter :: rho_earth_mean = 5511.0_wp!! mean density [kg m^-3]

   ! --- Material densities ------------------------------------------------------
   ! The defaults of the &vilma parameters rho_ice, rho_water and rho_seawater that
   ! the model actually uses (a host with its own densities sets those). Ice: the
   ! standard 917, not the GIA benchmarks' 931. Two waters (vilma_sle): melt water,
   ! fresh, whose volume a melt adds to the ocean (1000), and sea water, of the
   ! flotation test and the ocean columns that move (1028). The GIA benchmarks have
   ! one water of 1000 (rho_seawater = rho_water).
   real(wp), parameter :: rho_ice      = 917.0_wp   !! ice  [kg m^-3]
   real(wp), parameter :: rho_water    = 1000.0_wp  !! melt water (fresh) [kg m^-3]
   real(wp), parameter :: rho_seawater = 1028.0_wp  !! sea water [kg m^-3]

   ! --- Rotation (for the rotational-feedback / TPW module) -------------------
   real(wp), parameter :: omega_earth = 7.292115e-5_wp  !! mean rotation rate [rad s^-1]

   ! --- Time --------------------------------------------------------------------
   real(wp), parameter :: sec_per_year = 31556926.0_wp  !! seconds in a Julian year
   real(wp), parameter :: kyr          = 1000.0_wp*sec_per_year

end module vilma_constants
