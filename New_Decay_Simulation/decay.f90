program simulation
 use mod_functional
 use gplot
implicit none

real :: E1, Ex, Eth, Ecm1, Ecm2, Ecm3, E3, E4, theta3, theta4, phi3, phi4
real :: P1(3), P2(3), P3(3), P4(3)
real :: d2r, r2d, a , b , d , r , s, t, u
real :: Ealp1, Ealp2, Ealp3, E_Be, v1(3), v_1(3), v2(3), v_2(3), v3(3), v_3(3), v_Be(3), v__Be(3), v_C(3), v(3) 
real :: theta_1, phi_1, theta_2, phi_2, theta_3, phi_3, theta, phi, psi, kai, jai
real :: v11,v22,v33, factor, theta_phi3(2), e12, e23, e31, e(3)

integer :: i,j,N,N_hit
real :: theta_min, theta_max, phi_min, phi_max

d2r=0.017453292
r2d=1./d2r

E1=75.0 !Beam energy
N = 1000
!write(*,*)'ENTER MIN. & MAX. THETA OF SCATERRED 12C:'
!read(*,*)theta_min, theta_max
theta_min = 20.0
theta_max = 40.0
!write(*,*)'ENTER MIN. & MAX. PHI OF THE DETECTOR:'
!read(*,*)phi_min, phi_max
phi_min = -10.0
phi_max = 10.0

open(14,file='dalitz_seq.txt')
open(15,file='dalitz_ddl.txt')
open(16,file='dalitz_dde.txt')
open(17,file='dalitz_ddphi.txt')


Eth = 7.27

call randomEG(r)



do i=1,N
call random_normal(r)
Ex = 7.654 + r*0.0086   !Normally Distributed Excitation Energy with mean Exc. energy 7.654 MeV & sd = 200 keV


call random_number(r)
call random_number(s)
theta3 = theta_min + (theta_max - theta_min)*r  !Recoil 12C Polar Angle
phi3 = phi_min + (phi_max - phi_min)*s   !Recoil 12C Azimuthal Angle 

!Kinematics
a=1.
b=-sqrt(E1)*cos(theta3*d2r)
d=Ex/2

P1 = (/0.0,0.0,sqrt(2.*12.0*E1)/)
P2 = (/0.0,0.0,0.0/)

if ((b**2. - 4*a*d) .gt. 0. .or. (b**2. - 4*a*d) .eq. 0.) then
E3 = (0.5*(-b + sqrt(b**2. - 4*a*d)))**2.   !recoil 12C KE
P3 = sqrt(2.*12.*E3)*pol_cart(theta3,phi3)
!write(*,*)P3
!E4 = E1 - E3 - Ex
!P4 = sqrt(2.*12.*E4)
P4 = P1+P2-P3
v_C = (P4/12.)!*pol_cart(theta4, phi4) !Velocity of Excited Target 12C nucleus

!endif

!!!!!!!!!!!!!!!########## DECAY of 12C in DIFFERENT MODES ##############!!!!!!!!!!!!!!!!!!!


!## SEQUENTIALl DECAY OF 12C INTO 8Be+4He ##

Ecm2 = Ex - 7.36
Ealp1 = (2./3.)*Ecm2
E_Be = Ecm2 - Ealp1

!Randomly distribute Polar and Azimuthal Angles in CMF
call random_number(r)
call random_number(s)
!write(*,*)r,s
theta_1 = 180.*r   !Polar Angle of emitted alpha-1
phi_1 = 360.*s     !Azimuthal Angle of emitted alpha-1

v_1 = sqrt(2*Ealp1/4.)*pol_cart(theta_1,phi_1)
v__Be = - 0.5*v_1

!Convert CMF to LF
v1 = v_1 + v_C
v_Be = v__Be + v_C

!write(10,*)pol_angles(v1), pol_angles(v_Be)


!Breakup of 8Be in CMF
Ecm3 = (Ex - Eth) - Ecm2

Ealp2 = 0.5*Ecm3
Ealp3 = Ecm3 - Ealp2

!Randomly distribute Polar and Azimuthal Angles in CMF
call random_number(r)
call random_number(s)
!write(*,*)r,s
theta_2 = 180.*r   !Polar Angle of emitted alpha-1
phi_2 = 360.*s     !Azimuthal Angle of emitted alpha-1

v_2 = sqrt(2*Ealp2/4.)*pol_cart(theta_2,phi_2)
v_3 = - v_2

!Convert CMF to LF
v2 = v_2 + v_Be
v3 = v_3 + v_Be


!write(*,*)4*(v1+v2+v3), 12*v_C   !to check if momentum conservation was followed

!For Dalitz Plot
e12 = energy(2., (v1-v2))
e23 = energy(2., (v2-v3))
e31 = energy(2., (v3-v1))
e = (/e12,e23,e31/)
!print*,e
e = reverse(sort(e))
!print*,e
!e = rand_shuffle(e)
!write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.
write(14,*)dalitz(e)

endif
enddo 



do j=1,5000

call random_normal(r)
Ex = 7.654 + (-1.+2.*r)*0.0043    !Normally Distributed Excitation Energy with mean Exc. energy 7.654 MeV & sd = 200 keV


!if (Ex .gt. Eth) then


call random_number(r)
call random_number(s)
!write(*,*)r,s
!theta3 = (18.+20.*r)  !Recoil 12C Polar Angle
!phi3 = -16. + 32.*s   !Recoil 12C Azimuthal Angle 
theta3 = theta_min + (theta_max - theta_min)*r  !Recoil 12C Polar Angle
phi3 = phi_min + (phi_max - phi_min)*s   !Recoil 12C Azimuthal Angle 

!Kinematics
a=1.
b=-sqrt(E1)*cos(theta3*d2r)
d=Ex/2

P1 = (/0.0,0.0,sqrt(2.*12.0*E1)/)
P2 = (/0.0,0.0,0.0/)

if ((b**2. - 4*a*d) .gt. 0. .or. (b**2. - 4*a*d) .eq. 0.) then
E3 = (0.5*(-b + sqrt(b**2. - 4*a*d)))**2.   !recoil 12C KE
P3 = sqrt(2.*12.*E3)*pol_cart(theta3,phi3)
!write(*,*)P3
!E4 = E1 - E3 - Ex
!P4 = sqrt(2.*12.*E4)
P4 = P1+P2-P3
!phi4 = 180. + phi3  !No Momentum perpendicular to the plane so phi3-phi4=180
!call pol_cart(theta3, phi3)
v_C = (P4/12.)!*pol_cart(theta4, phi4) !Velocity of Excited Target 12C nucleus
!write(*,*)12.*(sqrt(v_C(1)**2.+v_C(2)**2.+v_C(3)**2.))
!endif



!## DECAY IN DDL MODE ##


Ecm1 = Ex - Eth
Ealp1 = Ecm1/2.
Ealp2 = 0.
Ealp3 = Ecm1/2.

call random_number(r)
call random_number(s)
!write(*,*)r,s
theta_1 = 180.*r   !Polar Angle of emitted alpha-1
phi_1 = 360.*s     !Azimuthal Angle of emitted alpha-1

v_1 = sqrt(2*Ealp1/4.)*pol_cart(theta_1,phi_1)
v_3 = -v_1
v_2 = - (4.*v_1 + 4.*v_3)/4.  !Baiscally (0,0,0) in CM, but following momentum conservation in centre of mass

!Convert CMF to LF
v1 = v_1 + v_C
v2 = v_2 + v_C
v3 = v_3 + v_C



!For Dalitz Plot
e12 = energy(2., (v1-v2))
e23 = energy(2., (v2-v3))
e31 = energy(2., (v3-v1))
e = (/e12,e23,e31/)
!e = sort(e)
!e = rand_shuffle(e)
!write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.
write(15,*)dalitz(e)

!## DECAY IN DDE MODE ##

call random_number(r)
call random_number(s)

Ecm1 = Ex - Eth
Ealp1 = Ecm1/3.0 + (-1.+2.*r)*0.0043
Ealp2 = Ecm1/3.0 + (-1.+2.*r)*0.0043
Ealp3 = Ecm1 - Ealp1 - Ealp2

!In CMF emission angles
call random_number(r)
call random_number(s)
call random_number(t)

psi = 360.0*r
kai = 360.0*s
jai = 360.0*t


theta_1 = 90.0
phi_1 = jai + 0.0

theta_2 = 90.0
phi_2 = phi_1 + 120.0

theta_3 = 90.0
phi_3 = phi_2 + 120.0 

!theta = 120.0
!phi = 120.0


v_1 = sqrt(2*Ealp1/4.)*pol_cart(theta_1,phi_1)
v_2 = sqrt(2*Ealp2/4.)*pol_cart(theta_2,phi_2)
v_3 = sqrt(2*Ealp3/4.)*pol_cart(theta_3,phi_3)
!print*,v_1+v_2+v_3


v_1 = rotate(v_1,1,psi)
v_2 = rotate(v_2,1,psi)
v_3 = rotate(v_3,1,psi)

v_1 = rotate(v_1,2,kai)
v_2 = rotate(v_2,2,kai)
v_3 = rotate(v_3,2,kai)


!v1 = rotate(v1,3,jai)
!v2 = rotate(v2,3,jai)
!v3 = rotate(v3,3,jai)

!print*,v_1+v_2+v_3

!write(*,*)energy(4., v_1) , energy(4., v_2) , energy(4., v_3)

!write(*,*)(v_1 + v_2 + v_3) !To check if the rotation was given properly
!write(*,*)energy(4., v_1) + energy(4., v_2) + energy(4., v_3)

!Convert CMF to LF
v1 = v_1 + v_C
v2 = v_2 + v_C
v3 = v_3 + v_C

!write(*,*)4*(v1+v2+v3), 12*v_C   !to check if momentum conservation was followed


!For Dalitz Plot
e12 = energy(2., (v1-v2))
e23 = energy(2., (v2-v3))
e31 = energy(2., (v3-v1))
e = (/e12,e23,e31/)
!e = sort(e)
e = rand_shuffle(e)
!write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.
write(16,*)dalitz(e)



!## DECAY IN DDphi MODE ##

Ecm1 = Ex - Eth

!Distributing energy between 3 alphas'
call random_number(r)
call random_number(s)
v11 = 0.5*r
v22 = (1-v11)*s
v33 = 1.0 - v22 - v11

if (v11+v22 .ge. v33 .and. v22+v33 .ge. v11 .and. v33+v11 .ge. v22) then
!In CMF emission angles
call random_number(r)
call random_number(s)
call random_number(t)

psi = 360.0*r
kai = 360.0*s
jai = 360.0*t

theta_1 = 90.0   !Polar Angle of emitted alpha-1
phi_1 = jai     !Azimuthal Angle of emitted alpha-1

theta_2 = 90.0
phi_2 = phi_1 + acos((v33**2.0 - v11**2.0 - v22**2.0)/(2.*v11*v22))*r2d

theta_3 = 90.0
phi_3 = phi_2 + acos((v11**2.0 - v22**2.0 - v33**2.0)/(2.*v33*v22))*r2d

!write(*,*)phi_1,phi_2,phi_3,phi_3+acos((v22**2.0 - v11**2.0 - v33**2.0)/(2.*v33*v11))*r2d
!write(*,*)v11*cos(phi_1*d2r)+v22*cos(phi_2*d2r)+v33*cos(phi_3*d2r)

v_1 = v11*pol_cart(theta_1,phi_1)
v_2 = v22*pol_cart(theta_2,phi_2)
v_3 = v33*pol_cart(theta_3,phi_3)

!print*,'Sum', v_1+v_2+v_3

v_1 = rotate(v_1,1,psi)
v_2 = rotate(v_2,1,psi)
v_3 = rotate(v_3,1,psi)

v_1 = rotate(v_1,2,kai)
v_2 = rotate(v_2,2,kai)
v_3 = rotate(v_3,2,kai)


!v1 = rotate(v1,3,jai)
!v2 = rotate(v2,3,jai)
!v3 = rotate(v3,3,jai)
!print*,'Sum', v_1+v_2+v_3

!write(*,*)(energy(4., v_1) + energy(4., v_2) + energy(4., v_3))


factor = sqrt(Ecm1/(energy(4., v_1) + energy(4., v_2) + energy(4., v_3)))

v_1 = factor*v_1
v_2 = factor*v_2
v_3 = factor*v_3

!write(*,*)(v_1 + v_2 + v_3)


!write(*,*)(energy(4., v_1) + energy(4., v_2) + energy(4., v_3))



!Convert CMF to LF
v1 = v_1 + v_C
v2 = v_2 + v_C
v3 = v_3 + v_C

!write(*,*)4*(v1+v2+v3), 12*v_C   !to check if momentum conservation was followed



!For Dalitz Plot
e12 = energy(2., (v1-v2))
e23 = energy(2., (v2-v3))
e31 = energy(2., (v3-v1))
e = (/e12,e23,e31/)
!e = sort(e)
e = rand_shuffle(e)
!write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.
write(17,*)dalitz(e)


endif

endif
enddo

 close(14)
 close(15)
 close(16)
 close(17)
!Plot theta and phi various decay modes in gnuplot
call execute_command_line('gnuplot -p dalitz_seq.plt')
call execute_command_line('gnuplot -p dalitz_ddl.plt')
!call execute_command_line('gnuplot -p theta_phi_DDE.plt')
!call execute_command_line('gnuplot -p theta_phi_DDphi.plt')






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
function pol_angles(vec)
real :: pol_angles(2)
real,intent(in) :: vec(3)
real :: theta, phi, d2r, r2d, r 

d2r=0.017453292
r2d=1./d2r
r = sqrt(vec(1)**2. + vec(2)**2. + vec(3)**2.)

pol_angles(1) = acos(vec(3)/r)*r2d

if (vec(2) .ge. 0.0) then
!pol_angles(2) = atan(vec(2)/vec(1))*r2d
pol_angles(2) = acos(vec(1)/(r*sin(acos(vec(3)/r))))*r2d
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


!Randomly Shuffle an array of 3 elements
function rand_shuffle(arr)
implicit none
real :: rand_shuffle(3), r, a, b, c
real, intent(in) :: arr(3)
integer :: i

a = arr(1)
b = arr(2)
c = arr(3)

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

