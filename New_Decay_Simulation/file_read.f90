program file_read
use, INTRINSIC :: ieee_arithmetic
implicit none

real :: a
integer :: i

open(10,file='dalitz_seq.txt')
i = 0

do while(i .le. 10000)
read(10,*)a
write(*,*)a,i
i = i+1
end do

 close(10)
end program
