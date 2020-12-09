program simulation
implicit none

real :: E1, Ex, Eth, Ecm1, Ecm2, Ecm3, E3, E4, theta3, theta4, phi3, phi4, P3, P4, d2r, r2d, a , b , d , r , s, t, u
real :: Ealp1, Ealp2, Ealp3, E_Be, v1(3), v_1(3), v2(3), v_2(3), v3(3), v_3(3), v_Be(3), v__Be(3), v_C(3), v(3) 
real :: theta_1, phi_1, theta_2, phi_2, theta_3, phi_3, theta, phi, psi, jai
real :: th_phi1(2), th_phi2(2), th_phi3(2)
real :: e12, e23, e31, xdalitz, ydalitz
real :: v11,v22,v33, factor, theta_phi3(2),x,y
!For Inputs
integer :: i, j, N, N_hit_seq, N_hit_DDL, N_hit_DDE
real :: theta_min, theta_max, phi_min, phi_max

d2r=0.017453292
r2d=1./d2r

E1=75.0 !Beam energy
N_hit_seq = 0
N_hit_DDL = 0
N_hit_DDE = 0

write(*,*)'ENTER THE NUMBER OF EVENTS:'
!read(*,*)N
N=100000
write(*,*)'ENTER MIN. & MAX. THETA OF SCATERRED 12C:'
!read(*,*)theta_min, theta_max
theta_min = 20.0
theta_max = 40.0
write(*,*)'ENTER MIN. & MAX. PHI OF THE DETECTOR:'
!read(*,*)phi_min, phi_max
phi_min = -16.0
phi_max = 16.0


open(10, file='theta_phi_seq.txt')
open(11, file='theta_phi_DDL.txt')
open(12, file='theta_phi_DDE.txt')
open(13, file='theta_phi_DDphi.txt')
open(14,file='dalitz_DDL.txt')


!write(*,*)'Scattered Energy', 'Recoil Energy', 'Recoil Angle', 'Exitation Energy', 'Scattered Angle'

Ex = 7.65      !excitation energy of the recoil 12C
Eth = 7.27

do i=1,N
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

if ((b**2. - 4*a*d) .gt. 0. .or. (b**2. - 4*a*d) .eq. 0.) then
E3 = (0.5*(-b + sqrt(b**2. - 4*a*d)))**2.   !recoil 12C KE
P3 = sqrt(2.*12.*E3)
!write(*,*)P3
E4 = E1 - E3 - Ex  !Scattered 12C Energy
P4 = sqrt(2.*12.*E4)
theta4 = (asin(P3*sin(theta3*d2r)/P4))*r2d
phi4 = 180. + phi3  !No Momentum perpendicular to the plane so phi3-phi4=180
!call pol_cart(theta3, phi3)
v_C = (P4/12.)*pol_cart(theta4, phi4) !Velocity of Excited Target 12C nucleus
!write(*,*)12.*(sqrt(v_C(1)**2.+v_C(2)**2.+v_C(3)**2.))
endif





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

!Breakup of 8Be in CMF
Ecm3 = 0.092
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

!Write into file
write(10,*)pol_angles(v1), energy(4.,v1)
write(10,*)pol_angles(v2), energy(4.,v2)
write(10,*)pol_angles(v3), energy(4.,v3)

!write(*,*)4*(v1+v2+v3), 12*v_C   !to check if momentum conservation was followed
!write(*,*)pol_angles(v1)

!write(*,*)energy(2., (v3-v2))














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

!Write into file
write(11,*)pol_angles(v1), energy(4.,v1)
write(11,*)pol_angles(v2), energy(4.,v2)
write(11,*)pol_angles(v3), energy(4.,v3)


!write(*,*)energy(2., (v1-v2))














!## DECAY IN DDE MODE ##

Ecm1 = Ex - Eth
Ealp1 = Ecm1/3.
Ealp2 = Ecm1/3.
Ealp3 = Ecm1 - Ealp1 - Ealp2

!In CMF emission angles
call random_number(r)
call random_number(s)
call random_number(t)
!write(*,*)r,s,t
!First taking the z=0 plane
theta_1 = 90.   !Polar Angle of emitted alpha-1
phi_1 = 360.*r*d2r     !Azimuthal Angle of emitted alpha-1

theta_2 = 90.
phi_2 = phi_1 + 120.*d2r

theta_3 = 90.
phi_3 = phi_2 + 120.*d2r

!For rotating the z=0 plane around x axis by psi and around y axis by jai
psi = 360.*s*d2r
jai = 360.*t*d2r

!v_1 = (/cos(phi_1), sin(phi_1), 0./)
!v_2 = (/cos(phi_2), sin(phi_2), 0./)
!v_3 = (/cos(phi_3), sin(phi_3), 0./)

v_1 = sqrt(2.*Ealp1/4.)*(/(cos(jai)*cos(phi_1)+sin(jai)*sin(psi)*sin(phi_1)), (cos(psi)*sin(phi_1)), &
(sin(jai)*cos(phi_1)- cos(jai)*sin(psi)*sin(phi_1))/)

v_2 = sqrt(2.*Ealp2/4.)*(/(cos(jai)*cos(phi_2)+sin(jai)*sin(psi)*sin(phi_2)), (cos(psi)*sin(phi_2)), &
(sin(jai)*cos(phi_2)- cos(jai)*sin(psi)*sin(phi_2))/)

v_3 = sqrt(2.*Ealp3/4.)*(/(cos(jai)*cos(phi_3)+sin(jai)*sin(psi)*sin(phi_3)), (cos(psi)*sin(phi_3)), &
(sin(jai)*cos(phi_3)- cos(jai)*sin(psi)*sin(phi_3))/)

!write(*,*)(v_1 + v_2 + v_3) !To check if the rotation was given properly

!Convert CMF to LF
v1 = v_1 + v_C
v2 = v_2 + v_C
v3 = v_3 + v_C

!write(*,*)4*(v1+v2+v3), 12*v_C   !to check if momentum conservation was followed

!Write into file
write(12,*)pol_angles(v1), energy(4.,v1)
write(12,*)pol_angles(v2), energy(4.,v2)
write(12,*)pol_angles(v3), energy(4.,v3)


!write(*,*)energy(2., (v1-v2))











!## DECAY IN DDphi MODE ##

Ecm1 = Ex - Eth

!Distributing energy between 3 alphas'
call random_number(r)
call random_number(s)
!x = -0.5 + 1.0*r
!y = -0.5 + 1.0*s
v11 = 0.5*r
v22 = (1-v11)*s
v33 = 1.0 - v22 - v11
!if (sqrt(x*x+y*y) .lt. 0.5) then
!write(*,*)x,y
!v11 = 1.0/3.0 - (1.0/3.0)*y - (1/sqrt(3.0))*x
!v22 = 1.0/3.0 - (1.0/3.0)*y + (1/sqrt(3.0))*x
!v33 = 1.0 - v22 - v11
!write(*,*)v11,v22,v33
if (v11+v22 .ge. v33 .and. v22+v33 .ge. v11 .and. v33+v11 .ge. v22) then
!In CMF emission angles
call random_number(r)
theta_1 = 90.0   !Polar Angle of emitted alpha-1
phi_1 = 360.*r     !Azimuthal Angle of emitted alpha-1

theta_2 = 90.0
phi_2 = phi_1 + acos((v33**2.0 - v11**2.0 - v22**2.0)/(2.*v11*v22))*r2d

theta_3 = 90.0
phi_3 = phi_2 + acos((v11**2.0 - v22**2.0 - v33**2.0)/(2.*v33*v22))*r2d

!write(*,*)phi_1,phi_2,phi_3,phi_3+acos((v22**2.0 - v11**2.0 - v33**2.0)/(2.*v33*v11))*r2d
!write(*,*)v11*cos(phi_1*d2r)+v22*cos(phi_2*d2r)+v33*cos(phi_3*d2r)

v_1 = v11*pol_cart(theta_1,phi_1)
v_2 = v22*pol_cart(theta_2,phi_2)
v_3 = v33*pol_cart(theta_3,phi_3)

!v_1 = v11*(/cos(phi_1), sin(phi_1), 0./)
!v_2 = v22*(/cos(phi_2), sin(phi_2), 0./)
!v_3 = v33*(/cos(phi_3), sin(phi_3), 0./)
!v_3 = - v_1 - v_2
!write(*,*)(energy(2., v_1) + energy(2., v_2) + energy(2., v_3))
!write(*,*)acos((v33**2.0 - v11**2.0 - v22**2.0)/(2.*v11*v22))*r2d,&
!acos((v22**2.0 - v11**2.0 - v33**2.0)/(2.*v33*v11))*r2d,&
!acos((v11**2.0 - v22**2.0 - v33**2.0)/(2.*v33*v22))*r2d
!theta_phi3 = pol_angles(v_3)
!theta_3 = theta_phi3(1)
!write(*,*)theta_3
!phi_3 = theta_phi3(2)

!write(*,*)(v_1 + v_2 + v_3)
!write(*,*)(energy(4., v_1) + energy(4., v_2) + energy(4., v_3))


factor = sqrt(Ecm1/(energy(4., v_1) + energy(4., v_2) + energy(4., v_3)))

v_1 = factor*v_1
v_2 = factor*v_2
v_3 = factor*v_3

!write(*,*)(v_1 + v_2 + v_3)

v11 = sqrt(v_1(1)**2.0 + v_1(2)**2.0 + v_1(3)**2.0)
v22 = sqrt(v_2(1)**2.0 + v_2(2)**2.0 + v_2(3)**2.0)
v33 = sqrt(v_3(1)**2.0 + v_3(2)**2.0 + v_3(3)**2.0)

!write(*,*)(energy(4., v_1) + energy(4., v_2) + energy(4., v_3))


!For rotating the z=0 plane around x axis by psi and around y axis by jai
call random_number(s)
call random_number(t)
psi = 360.*s*d2r
jai = 360.*t*d2r


v_1 = (/(cos(jai)*v_1(1)+sin(jai)*sin(psi)*v_1(2)), (cos(psi)*v_1(2)), &
(sin(jai)*v_1(1)- cos(jai)*sin(psi)*v_1(2))/)

v_2 = (/(cos(jai)*v_2(1)+sin(jai)*sin(psi)*v_2(2)), (cos(psi)*v_2(2)), &
(sin(jai)*v_2(1)- cos(jai)*sin(psi)*v_2(2))/)

v_3 = (/(cos(jai)*v_3(1)+sin(jai)*sin(psi)*v_3(2)), (cos(psi)*v_3(2)), &
(sin(jai)*v_3(1)- cos(jai)*sin(psi)*v_3(2))/)

!write(*,*)(v_1 + v_2 + v_3) !To check if the rotation was given properly
!write(*,*)(energy(4., v_1) + energy(4., v_2) + energy(4., v_3))

!Convert CMF to LF
v1 = v_1 + v_C
v2 = v_2 + v_C
v3 = v_3 + v_C

!write(*,*)4*(v1+v2+v3), 12*v_C   !to check if momentum conservation was followed

!Write into file
write(13,*)pol_angles(v1), energy(4.,v1)
write(13,*)pol_angles(v2), energy(4.,v1)
write(13,*)pol_angles(v3), energy(4.,v1)

!For Dalitz Plot
e12 = energy(2., (v1-v2))
e23 = energy(2., (v2-v3))
e31 = energy(2., (v3-v1))

!write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.

endif



!End Of Do loop
enddo


 close(10)
 close(11)
 close(12)
 close(13)
 close(14)


!Plot theta and phi various decay modes in gnuplot
!call execute_command_line('gnuplot -p theta_phi_seq.plt')
!call execute_command_line('gnuplot -p theta_phi_DDL.plt')
!call execute_command_line('gnuplot -p theta_phi_DDE.plt')
call execute_command_line('gnuplot -p theta_phi_DDphi.plt')
!call execute_command_line('./eff_seq') !Calculate Efficiency also









contains   !All the required functions

!! Polar to Cartesian Conversion
function pol_cart(theta,phi)
real :: pol_cart(3)
real :: theta, phi, d2r, r2d 

d2r=0.017453292
r2d=1./d2r

theta = theta*d2r
phi = phi*d2r

pol_cart(1) = sin(theta)*cos(phi)
pol_cart(2) = sin(theta)*sin(phi)
pol_cart(3) = cos(theta)
return 
end function






!!Find Theta and Phi from vector Components
function pol_angles(vec)
real :: pol_angles(2)
real :: vec(3)
real :: theta, phi, d2r, r2d, r 

d2r=0.017453292
r2d=1./d2r
r = sqrt(vec(1)**2. + vec(2)**2. + vec(3)**2.)

pol_angles(1) = acos(vec(3)/r)*r2d
pol_angles(2) = atan(vec(2)/vec(1))*r2d

return
end function





!!Find Energy of a particle given its velocity vector
function energy(mass, vec)
real :: energy
real :: vec(3)
real mass
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



end program

