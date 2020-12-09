module func_subs

    private

public :: pol_cart, pol_angles, rel_energy, rel_angle,  rotate, rand_shuffle, sort, energy, dalitz, rms, random_normal, randomEG

interface pol_cart
  module procedure :: pol_cart
end interface pol_cart

interface pol_angles
  module procedure :: pol_angles
end interface pol_angles

interface rel_energy
  module procedure :: rel_energy
end interface rel_energy

interface rel_angle
  module procedure :: rel_angle
end interface rel_angle

interface rotate
  module procedure :: rotate
end interface rotate

interface rand_shuffle
  module procedure :: rand_shuffle
end interface rand_shuffle

interface sort
  module procedure :: sort
end interface sort

interface energy
  module procedure :: energy
end interface energy

interface dalitz
  module procedure :: dalitz
end interface dalitz

interface rms
  module procedure :: rms
end interface rms

interface random_normal
  module procedure :: random_normal
end interface random_normal

interface randomEG
  module procedure :: randomEG
end interface randomEG




contains   !All the required functions
   
   
!! Polar to Cartesian Conversion
function pol_cart(theta,phi)
real :: pol_cart(3)
real,intent(in) :: theta, phi
real :: d2r, r2d ,theta1, phi1
        
d2r=0.017453292
r2d=1./d2r
        
theta1 = theta*d2r
phi1 = phi*d2r
        
pol_cart(1) = sin(theta1)*cos(phi1)
pol_cart(2) = sin(theta1)*sin(phi1)
pol_cart(3) = cos(theta1)
return 
end function






!!Find Theta and Phi from vector Components
subroutine pol_angles(vec,theta,phi)
!real :: pol_angles(2)
use, intrinsic :: ieee_arithmetic
real,intent(in) :: vec(3)
real,intent(out) :: theta, phi
real :: d2r, r2d, r 

d2r=0.017453292
r2d=1./d2r
r = sqrt(vec(1)**2. + vec(2)**2. + vec(3)**2.)

theta = acos(vec(3)/r)*r2d

if (vec(1) .gt. 0.0 .and. vec(2) .gt. 0.0) then
phi = atan(vec(2)/vec(1))*r2d
!phi = acos(vec(1)/(r*sin(acos(vec(3)/r))))*r2d
elseif (vec(1) .lt. 0.0 .and. vec(2) .gt. 0.0) then
phi = atan(vec(2)/vec(1))*r2d + 180.
elseif (vec(1) .lt. 0.0 .and. vec(2) .lt. 0.0) then
phi = atan(vec(2)/vec(1))*r2d + 180.
else
phi = 360. + atan(vec(2)/vec(1))*r2d
endif

if (ieee_is_nan(phi)) then
print*,vec,r
endif
return
end subroutine

!Function for Calculating Relative Energy
function rel_energy(m1, m2, E1, E2, gamma)
    real :: rel_energy
    real, intent(in) :: m1, m2, E1, E2, gamma
    real :: d2r, r2d

    d2r=0.017453292
    r2d=1./d2r

    rel_energy = (m2*E1 + m1*E2 - 2.0*sqrt(m1*m2*E1*E2)*cos(gamma*d2r))/(m1 + m2)

    return
   end function rel_energy



   !Function to Caclculate Relative Angle Between two Vectors
   function rel_angle(theta1, phi1, theta2, phi2)
    !use, intrinsic :: ieee_arithmetic
    real :: rel_angle
    real, intent(in) :: theta1, phi1, theta2, phi2
    real :: d2r, r2d, th1, th2, ph1, ph2

    d2r=0.017453292
    r2d=1./d2r

    th1 = theta1*d2r
    th2 = theta2*d2r
    ph1 = phi1*d2r
    ph2 = phi2*d2r

    rel_angle = acos(cos(th1)*cos(th2) + sin(th1)*sin(th2)*cos(ph1-ph2))*r2d
    !if (ieee_is_nan(rel_angle)) then
    !print*, theta1,theta2,phi1,phi2
    !endif
    return
    end function rel_angle




!!Rotation Matrix
function rotate(vec,axis,angle)
real :: rotate(3)
real,intent(in) :: vec(3)
real,intent(in) :: angle
integer,intent(in) :: axis
real :: r2d, d2r, angle1

d2r=0.017453292
r2d=1./d2r

angle1 = angle*d2r

if (axis == 1) then       !# For Rotation around X axis
rotate(1) = vec(1)
rotate(2) = vec(2)*cos(angle1) - vec(3)*sin(angle1)
rotate(3) = vec(2)*sin(angle1) + vec(3)*cos(angle1)
elseif (axis == 2) then   !# For Rotation around Y axis
rotate(1) = vec(1)*cos(angle1) + vec(3)*sin(angle1)
rotate(2) = vec(2)
rotate(3) = -vec(1)*sin(angle1) + vec(3)*cos(angle1)
elseif (axis == 3) then   !# For Rotation around Z axis
rotate(1) = vec(1)*cos(angle1) - vec(2)*sin(angle1)
rotate(2) = vec(1)*sin(angle1) + vec(2)*cos(angle1)
rotate(3) = vec(3)
endif

return
end function




!!Find Energy of a particle given its velocity vector
function energy(mass, vec)
real :: energy
real,intent(in) :: vec(3), mass

energy = 0.5*mass*(vec(1)**2. + vec(2)**2. + vec(3)**2.)

return
end function


!Randomly Shuffle an array of 3 elements
function rand_shuffle(arr)
implicit none
real :: rand_shuffle(9), r, a(3), b(3), c(3)
real, intent(in) :: arr(9)
integer :: i

a = arr(1:3)
b = arr(4:6)
c = arr(7:9)

!call randomEG(r)
call random_number(r)

i = int(r*6)
!print*,i

select case(i)
    case(0)
        rand_shuffle = (/a,b,c/)
    case(1)
        rand_shuffle = (/a,c,b/)
    case(2)
        rand_shuffle = (/b,a,c/)
    case(3)
        rand_shuffle = (/c,a,b/)
    case(4)
        rand_shuffle = (/b,c,a/)
    case(5)
        rand_shuffle = (/c,b,a/)
end select

return
end function rand_shuffle


!Sorting Subroutine
recursive subroutine sort(a, first, last)
  implicit none
  real  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
!print*,(first+last) / 2
  i = first
  j = last
  do
     do while (a(i) > x)
        i=i+1
     end do
     do while (x > a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call sort(a, first, i-1)
  if (j+1 < last)  call sort(a, j+1, last)
end subroutine sort




!!Dalitz Plot Coordinate Transformation
function dalitz(e)
real :: dalitz(2)
real, intent(in) :: e(3)
real :: e12, e23, e31

e12 = e(1)
e23 = e(2)
e31 = e(3)

dalitz(1) = sqrt(3.0)*(e12 - e23)/2.0
dalitz(2) = (2.0*e31 - e12 - e23)/2.0
return
end function


function rms(vec)
  implicit none
  
  real :: rms, mean, sqmean, variance
  real, intent(in) :: vec(3)
  
  mean = sum(vec)/3.0
  sqmean = (vec(1)**2.0 + vec(2)**2.0 + vec(3)**2.0)/3.0
  !variance = (mean-vec(1))**2.0 + (mean-vec(2))**2.0 + (mean-vec(3))**2.0
  rms = sqrt(sqmean - mean**2.0)
  return
 end function rms



!!Normal Random Number Generator
SUBROUTINE random_normal(r) 
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.   This version uses the default
!  uniform random number generator which is in your fortran library.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

!  Fortran 90 version by Alan Miller (alan @ mel.dms.csiro.au)

IMPLICIT NONE
REAL :: ran_norm, r

!     Local variables
REAL, PARAMETER :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,   &
                   half = 0.5, r1 = 0.27597, r2 = 0.27846
REAL            :: u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
ran_norm = v/u
r = ran_norm
RETURN

END SUBROUTINE random_normal






!Random Seed Generator
subroutine randomEG()
implicit none

integer :: k, i, n=10
real :: r
integer, dimension(8) :: values 
! Declare an assumed shape, dynamic array
integer, dimension(:), allocatable :: seed

! gfortran subroutine to return date and time information 
! from the real time system clock. Works down to milliseconds 
! and stores the eight return values in array values.
call date_and_time(VALUES=values)
! restart the state of the pseudorandom number generator
! k = minimum size of seed (12 on my system)
call random_seed(size=k)
! allocate memory to seed
allocate(seed(k))

! assign information in values to seed
seed(:) = values(:)
! seed the random number generator
call random_seed(put=seed)

!do i=1,n
!    call random_number(r)
!    print *, r
!end do
return
end subroutine randomEG


end module func_subs
