program simulation
    !use mod_functional
    !use gplot
    use func_subs
   implicit none
   
   real :: E1, Ex, Eth, Ecm1, Ecm2, Ecm3, E3, E4, theta3, theta4, phi3, phi4
   real :: P1(3), P2(3), P3(3), P4(3)
   real :: d2r, r2d, a , b , d , r , s, t, u
   real :: Ealp1, Ealp2, Ealp3, E_Be, v1(3), v_1(3), v2(3), v_2(3), v3(3), v_3(3), v_Be(3), v__Be(3), v_C(3), v(3) 
   real :: theta_1, phi_1, theta_2, phi_2, theta_3, phi_3, theta, phi, psi, kai, jai
   real :: v11,v22,v33, factor, theta_phi3(2), e12, e23, e31, e(3), En1, En2, gamma, th1 ,th2, ph1, ph2
   
   integer :: i,j,N,N_hit
   real :: theta_min, theta_max, phi_min, phi_max
   
   call randomEG()

   d2r=0.017453292
   r2d=1./d2r
   
   call random_normal(r)
   E1=75.0 + 0.002*r!Beam energy
   N = 100
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

   
   
   
   do i=1,N
   call random_normal(r)
   Ex = 7.654 + r*0.0000086   !Normally Distributed Excitation Energy with mean Exc. energy 7.654 MeV & sd = 200 keV
   
   
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
   !e12 = energy(2., (v1-v2))
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v1) + 0.04*r
   En2 = energy(4.0, v2) + 0.04*s
   call pol_angles(v1, th1, ph1)
   call pol_angles(v2, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e12 = rel_energy(4.0, 4.0, En1, En2 , gamma)

   !e23 = energy(2., (v2-v3))
   !print*, e23
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v2) + 0.04*r
   En2 = energy(4.0, v3) + 0.04*s
   call pol_angles(v2, th1, ph1)
   call pol_angles(v3, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e23 = rel_energy(4.0, 4.0, En1, En2 , gamma)
   !print*, e23
   !e31 = energy(2., (v3-v1))
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v3) + 0.04*r
   En2 = energy(4.0, v1) + 0.04*s
   call pol_angles(v3, th1, ph1)
   call pol_angles(v1, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e31 = rel_energy(4.0, 4.0, En1, En2 , gamma)

   e = (/e12,e23,e31/)
   !print*,e
   !e = reverse(sort(e))
   !print*,minval(e)
   e = rand_shuffle(e)
   !write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.
   write(14,*)dalitz(e),minval(e),rms(e)
   
   endif
   enddo 
   
   



   !Direct Decay Components
   do j=1,1000
   
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
   !e12 = energy(2., (v1-v2))
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v1) + 0.04*r
   En2 = energy(4.0, v2) + 0.04*s
   call pol_angles(v1, th1, ph1)
   call pol_angles(v2, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e12 = rel_energy(4.0, 4.0, En1, En2 , gamma)

   !e23 = energy(2., (v2-v3))
   !print*, e23
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v2) + 0.04*r
   En2 = energy(4.0, v3) + 0.04*s
   call pol_angles(v2, th1, ph1)
   call pol_angles(v3, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e23 = rel_energy(4.0, 4.0, En1, En2 , gamma)
   !print*, e23
   !e31 = energy(2., (v3-v1))
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v3) + 0.04*r
   En2 = energy(4.0, v1) + 0.04*s
   call pol_angles(v3, th1, ph1)
   call pol_angles(v1, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e31 = rel_energy(4.0, 4.0, En1, En2 , gamma)

   e = (/e12,e23,e31/)
   !print*,e
   !e = reverse(sort(e))
   !print*,e
   e = rand_shuffle(e)
   !write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.
   write(15,*)dalitz(e),minval(e),rms(e)
   




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
   !e12 = energy(2., (v1-v2))
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v1) + 0.04*r
   En2 = energy(4.0, v2) + 0.04*s
   call pol_angles(v1, th1, ph1)
   call pol_angles(v2, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e12 = rel_energy(4.0, 4.0, En1, En2 , gamma)

   !e23 = energy(2., (v2-v3))
   !print*, e23
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v2) + 0.04*r
   En2 = energy(4.0, v3) + 0.04*s
   call pol_angles(v2, th1, ph1)
   call pol_angles(v3, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e23 = rel_energy(4.0, 4.0, En1, En2 , gamma)
   !print*, e23
   !e31 = energy(2., (v3-v1))
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v3) + 0.04*r
   En2 = energy(4.0, v1) + 0.04*s
   call pol_angles(v3, th1, ph1)
   call pol_angles(v1, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e31 = rel_energy(4.0, 4.0, En1, En2 , gamma)

   e = (/e12,e23,e31/)
   !print*,e
   !e = reverse(sort(e))
   !print*,e
   e = rand_shuffle(e)
   !write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.
   write(16,*)dalitz(e),minval(e),rms(e)
   
   
   
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
   !e12 = energy(2., (v1-v2))
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v1) + 0.04*r
   En2 = energy(4.0, v2) + 0.04*s
   call pol_angles(v1, th1, ph1)
   call pol_angles(v2, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e12 = rel_energy(4.0, 4.0, En1, En2 , gamma)

   !e23 = energy(2., (v2-v3))
   !print*, e23
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v2) + 0.04*r
   En2 = energy(4.0, v3) + 0.04*s
   call pol_angles(v2, th1, ph1)
   call pol_angles(v3, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e23 = rel_energy(4.0, 4.0, En1, En2 , gamma)
   !print*, e23
   !e31 = energy(2., (v3-v1))
   call random_normal(r)
   call random_normal(s)
   En1 = energy(4.0, v3) + 0.04*r
   En2 = energy(4.0, v1) + 0.04*s
   call pol_angles(v3, th1, ph1)
   call pol_angles(v1, th2, ph2)
   call random_normal(r)
   call random_normal(s)
   call random_normal(t)
   call random_normal(u)
   th1 = th1 + 1.0*r
   ph1 = ph1 + 1.0*r
   th2 = th2 + 1.0*r
   ph2 = ph2 + 1.0*r 
   gamma = rel_angle(th1,ph1,th2,ph2)
   e31 = rel_energy(4.0, 4.0, En1, En2 , gamma)

   e = (/e12,e23,e31/)
   !print*,e
   !e = reverse(sort(e))
   !print*,e
   e = rand_shuffle(e)
   !write(14,*)(sqrt(3.)*(e12 - e23)/2.), (2.*e31 - e12 - e23)/2.
   write(17,*)dalitz(e),minval(e),rms(e)
   
   
   endif
   
   endif
   enddo
   
    close(14)
    close(15)
    close(16)
    close(17)
   !Plot theta and phi various decay modes in gnuplot
   !call execute_command_line('gnuplot -p dalitz_seq.plt')
   !call execute_command_line('gnuplot -p dalitz_ddl.plt')
   !call execute_command_line('gnuplot -p dalitz_dde.plt')
   !call execute_command_line('gnuplot -p dalitz_ddphi.plt')

   
    
   !contains

   end program
   
   