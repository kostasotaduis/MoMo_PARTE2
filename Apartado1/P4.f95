program P4

use modP4

implicit none

real*8, dimension(:,:), allocatable :: r, v, rnew, vnew !variables array
real*8 :: rho, Temp, dt!parametros real
real*8 :: t !variables real
integer*4 :: N, ktotal !parametros integer
integer*4 :: i, j, k !auxiliar integer
real*8 :: rnd, n1, n2 !auxiliar real
real*8 :: x

!La random seed_____________________________
integer, dimension(:), allocatable :: state !
integer :: state_size                       !
integer :: cput                             !
                                            !
call random_seed(size=state_size)           !
allocate(state(state_size))                 !
call system_clock(cput)                     !
state = cput                                !
call random_seed(put=state)                 !
!Acabo random seed__________________________!

!setear parametros
N = 125
m = 1.0d0
sigma = 1.0d0
rho = 0.7d0
L = N**(1.0d0/3.0d0) * sigma * rho**(-1.0d0/3.0d0)
!L = 10.0
!L = 10.0d0
rc = 4.0d0 * sigma
Temp = 100.0d0

dt = 0.001d0
ktotal = 100000

!open(unit=0, file="prueba.dat", status="replace", action="write")

!colocar particulas en square lattice
call latsc(N, rho, sigma, r)

!allocatear resto de arrays
allocate(rnew(N,3), v(N,3), vnew(N,3))
v = 0.0d0
!iniciar velocidades que correspondan a una temperatura

do i = 1, N
	do j = 1, 3
		call random_number(rnd)
		!call BoxMuller(((kB * Temp)**0.5), n1, n2)
		v(i,j) = 10.0d0 * (floor(rnd*2.0d0)*2-1.0d0)
	end do
end do

open (unit=2, file="vel0.dat", status="replace", action="write")

do i = 1, N
	write(unit=2, fmt=*) norm2(v(i,:))
end do

print*, L

!iniciar simulacion
open (unit=0, file="p4_Eulerdt0001.dat", status="replace", action="write")
!open (unit=1, file="p4b.xyz", status="replace", action="write")
rnew = r
vnew = v
t = 0.0
write(unit=0, fmt=*) t, Ek(v), ELenardpbc(r)
do k = 1, ktotal
	!call vVerlet(r, v, dt, rnew, vnew)
	call Euler(r, v, dt, rnew, vnew)
	r = rnew
	v = vnew
	do i = 1, N
		call pbc(r(i,:))
	end do
	t = t + dt
	if ((k/1000)*1000 == k) then
		print*, real(k)/real(ktotal)*100, "% completado"
		write(unit=0, fmt=*) t, Ek(v), ELenardpbc(r), norm2(ptot(v))
		!write(unit=1, fmt="(i0)") N
		!write(unit=1, fmt=*)
		!do i = 1, N
		!	write(unit=1, fmt="(a, f0.4, a, f0.4, a, f0.4)") "H ", r(i,1), " ", r(i,2), " ", r(i,3)
		!end do
	end if
end do

open (unit=2, file="p4b_dt000001.dat", status="replace", action="write")

do i = 1, N
	write(unit=2, fmt=*) norm2(v(i,:))
end do

close(unit=0)
close(unit=2)

end program P4