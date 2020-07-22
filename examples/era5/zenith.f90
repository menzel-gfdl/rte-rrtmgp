module zenith_mod
use mo_rte_kind, only: wp
implicit none


public :: gauss


contains


!> @brief Finds sunrise time as a fraction of a twelve hour period
!!         and then finds zeniths at three times in between sunrise and miday
!!         the times chosen, so gaussian integration can be used
!!         Maths of phys +chemistry Margenau+ Murphy page 462
!!
!!         to find the zenith agles the approx from Paltrige +platt p63 is used
!!         cosZ = sin(decl)sin(lat)+cos(decl)cos(lat)cos(Hangle)
subroutine gauss(day, lat_in, weight, z)

  integer, intent(in) :: day !< Julian day (1 = January 1st, and so on).
  real(kind=wp), intent(in) :: lat_in !< Latitude [degrees].
  real(kind=wp), dimension(3), intent(out) :: weight !< Gaussian weights.
  real(kind=wp), dimension(3), intent(out) :: z !< Cosine of zenith angles for gaussian weights.

  real(kind=wp) :: ca
  real(kind=wp) :: d
  real(kind=wp) :: doya
  real(kind=wp), dimension(3):: g
  real(kind=wp) :: h
  integer :: ifrac
  real(kind=wp) :: lat
  real(kind=wp), parameter :: pi = asin(1.0)*2.0
  real(kind=wp) :: sfrac

  !Gaussian integrations fractions are:
  g(1) = 0.11270166
  g(2) = 0.5
  g(3) = 0.88729833

  !Convert latitude to radians.
  lat = pi*lat_in/180.0

  doya = 2.*pi*float(day)/365.0
  d = 0.006918 - 0.399912*cos(doya) + 0.070257*sin(doya) &
      - 0.006758*cos(2.*doya) + 0.000907*sin(2.*doya) &
      - 0.002697*cos(3.*doya) + 0.001480*sin(3.*doya)

  !Find hangle at sunrise (cosZ = 0).
  if (d .ge. 0.0) then
    !Northern summer,
    if (abs(lat) .lt. ((pi/2.0) - d)) then
      h = acos(-(sin(d)*sin(lat))/(cos(d)*cos(lat)))
    elseif (lat .gt. 0.0) then
      !Polar day,
      h = pi
    else
      !Polar night,
      h = 0.0
    endif
  else
    !Southern summer.
    if (abs(lat) .lt. abs(-(pi/2.0) - d)) then
      h = acos(-(sin(d)*sin(lat))/(cos(d)*cos(lat)))
    elseif (lat .lt. 0.0) then
      !Polar day.
      h = pi
    else
      !Polar night.
      h = 0.0
    endif
  endif

  !Fraction of half day with sun up (or whole day).
  sfrac = h/pi

  !Compute zenith at gaussian fractions of time between sunrise and midday.
  do ifrac = 1, 3
    !Find hour angle at each fraction.
    ca = pi*(sfrac - g(ifrac)*sfrac)
    h = ca
    z(ifrac) = sin(d)*sin(lat) + cos(d)*cos(lat)*cos(h)
  enddo

  !Compute weights
  !equation of time comes in here though as the number
  !of seconds in each day is alterd by Eq
  weight(1) = sfrac*5.0/18.0
  weight(2) = sfrac*4.0/9.0
  weight(3) = sfrac*5.0/18.0
  weight(:) = weight(:)*(1.0 + 0.0167*cos(2.0*pi*(real(day) - 3)/365.0))**2
end subroutine gauss


end module zenith_mod
