!Decay Simulation for 16O+16O at beam energy 95 MeV
program simulation
implicit none

real :: m1, m2, m3, m4
real :: E1, Ex, Eth, Ecm1, Ecm2, Ecm3, E3, E4, theta3, theta4, phi3, phi4, d2r, r2d, a , b , d , r , s, t, u
real :: P1(3), P2(3), P3(3), P4(3), En1(5000), En2(5000), En3(5000), EnSc(5000), EnO(5000)
real :: Ealp1, Ealp2, Ealp3, E_Be, E_C
real :: v1(3), v_1(3), v2(3), v_2(3), v3(3), v_3(3), v4(3), v_4(3), v_Be(3), v__Be(3), v_C(3), v__C(3), v_O(3), v(3) 
real :: theta_1, phi_1, theta_2, phi_2, theta_3, phi_3, theta, phi, psi, kai, jai
real :: v11,v22,v33, factor, theta_phi3(2), e12, e23, e31, e(3)

integer :: i,j,N,N_hit
real :: theta_min, theta_max, phi_min, phi_max

d2r=0.017453292
r2d=1./d2r

E1 = 95.0 !Beam energy
N = 5000

theta_min = 25.
theta_max = 45.

phi_min = -5.0
phi_max = 5.0

open(10, file='theta_phi_seq.txt')
open(11, file='theta_phi_DDL.txt')
open(12, file='theta_phi_DDE.txt')
open(13, file='theta_phi_DDphi.txt')

!System Info
m1 = 16.0 !In amu
m2 = 16.0 !In amu
m3 = 16.0 !In amu
m4 = 16.0 !In amu

Eth = 14.44 !4-alpha breakup threshold

!N_hit = 1

call randomEG(r)  !Initialising the random generator

i = 1
do while (i <= N)   !Starting of the main Monte Carlo Loop

!call random_normal(r)
Ex = 15.1 !+ r*0.02    !Normally Distributed Excitation Energy with mean Exc. energy 7.654 MeV & sd = 200 keV

call random_number(r)
call random_number(s)
!write(*,*)r,s
!theta3 = 41.01  !Recoil 12C Polar Angle
!phi3 = 0.0   !Recoil 12C Azimuthal Angle 
theta4 = theta_min + (theta_max - theta_min)*r  !scattered 12C Polar Angle
phi4 = phi_min + (phi_max - phi_min)*s   !scattered 12C Azimuthal Angle 

!Kinematics : Where particle 4 is scattered and particle 3 is recoil and We are finding out the energy of the scattered particle
a=(1.0+m4/m3)
b=-2.0*sqrt(E1*(m1*m4/m3**2.0))*cos(theta4*d2r)
d=(Ex - (1-m1/m3)*E1)

   P1 = (/0.0,0.0,sqrt(2.*m1*E1)/)
   P2 = (/0.0,0.0,0.0/)

if ((b**2. - 4.*a*d) .ge. 0.) then
    !Write(*,*)N_hit
    !N_hit = N_hit + 1
    E4 = ((-b + sqrt(b**2. - 4*a*d))/(2.0*a))**2.   !scattered 12C KE
    P4 = sqrt(2.*m4*E4)*pol_cart(theta4,phi4)
    E3 = E1 - E4 - Ex
    P3 = P1+P2-P4
    v_O = (P3/m3)  !Velocity of recoil Excited 16O nucleus

EnSc(i) = energy(m4,P4/m4)
EnO(i) = energy(m3,v_O)
i = i+1
else
write(*,*)(b**2. - 4.*a*d)
endif

enddo


write(*,*)"min scattered 16O", minval(EnSc), "max scattered 16O", maxval(EnSc)
write(*,*)"min excited 16O", minval(EnO), "max excited 16O", maxval(EnO)
!Write(*,*)N_hit


contains  !Functions and subroutines

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
function pol_angles(vec)
real :: pol_angles(2)
real,intent(in) :: vec(3)
real :: theta, phi, d2r, r2d, r 

d2r=0.017453292
r2d=1./d2r
r = sqrt(vec(1)**2. + vec(2)**2. + vec(3)**2.)

pol_angles(1) = acos(vec(3)/r)*r2d

if (vec(1) .gt. 0.0 .and. vec(2) .gt. 0.0) then
pol_angles(2) = atan(vec(2)/vec(1))*r2d
!phi = acos(vec(1)/(r*sin(acos(vec(3)/r))))*r2d
elseif (vec(1) .lt. 0.0 .and. vec(2) .gt. 0.0) then
pol_angles(2) = atan(vec(2)/vec(1))*r2d + 180.
elseif (vec(1) .lt. 0.0 .and. vec(2) .lt. 0.0) then
pol_angles(2) = atan(vec(2)/vec(1))*r2d + 180.
else
pol_angles(2) = 360. + atan(vec(2)/vec(1))*r2d
endif

return
end function





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





!!Dalitz Plot Coordinate Transformation
function dalitz(e12, e23, e31)
real :: dalitz(2)
real :: e12, e23, e31

dalitz(1) = sqrt(3.)*(e12 - e23)
dalitz(2) = (2.*e31 - e12 - e23)
return
end function








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
subroutine randomEG(r)
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

end program

