program test
implicit none

real :: v(5),v1(3),v2(3)
real :: th,ph

v = (/0.2,7.1,6.,150.,5./)
v1 = (/-0.1,0.1,0.01/)
v2 = (/0.1,-0.1,-0.01/)

!print*,sqrt(energy(2.,v))
call pol_angles(v,th,ph)
!print*,sqrt(energy(2.,v))*pol_cart(th,ph)
call sort(v,1,5)
print*,v

contains

!!Find Energy of a particle given its velocity vector
function energy(mass, vec)
    real :: energy
    real,intent(in) :: vec(3), mass
    
    energy = 0.5*mass*(vec(1)**2. + vec(2)**2. + vec(3)**2.)
    
    return
end function

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



!Function for sorting an array of 3 elements
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
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
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


!!Find Theta and Phi from vector Components
subroutine pol_angles(vec,theta,phi)
!real :: pol_angles(2)
real,intent(in) :: vec(3)
real,intent(out) :: theta, phi
real :: d2r, r2d, r 

d2r=0.017453292
r2d=1./d2r
r = sqrt(vec(1)**2. + vec(2)**2. + vec(3)**2.)

theta = acos(vec(3)/r)*r2d

if (vec(2) .ge. 0.0) then
!pol_angles(2) = atan(vec(2)/vec(1))*r2d
phi = acos(vec(1)/(r*sin(acos(vec(3)/r))))*r2d
elseif (vec(1) .lt. 0.0 .and. vec(2) .lt. 0.0) then
phi = atan(vec(2)/vec(1))*r2d + 180.
else
phi = 360. + atan(vec(2)/vec(1))*r2d
endif

return
end subroutine

end program
