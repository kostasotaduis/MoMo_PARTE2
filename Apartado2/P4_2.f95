program P4

use modP4

implicit none

real*8, dimension(:,:), allocatable :: r, v, rnew, vnew !variables array
real*8 :: rho, Temp, dt, nu!parametros real
real*8 :: t, avEk, avEp, p, avp !variables real
integer*4 :: N, ktotal, keq !parametros integer
integer*4 :: i, j, k !auxiliar integer
real*8 :: rnd, n1, n2 !auxiliar real

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
N      = 125
m      = 1.0d0 !realmente no la vamos a usar
sigma  = 1.0d0
rho    = 0.8d0
L      = N**(1.0d0/3.0d0) * sigma * rho**(-1.0d0/3.0d0)
rc     = 4.0d0 * sigma
nu     = 10.0d0 !potencia del termostato

dt = 0.00001d0

!open(unit=0, file="prueba.dat", status="replace", action="write")

!colocar particulas en square lattice
call latsc(N, rho, sigma, r)

!allocatear resto de arrays
allocate(rnew(N,3), v(N,3), vnew(N,3))
v = 0.0d0

!colocar simulacion en estado fluido (termostato a 50kbe)
ktotal = 20000
Temp = 50.0d0
do k = 1, ktotal
	call vVerlet(r, v, dt, rnew, vnew)
	r = rnew
	v = vnew
	do i = 1, N
		call pbc(r(i,:))
	end do
	call Ander(v, Temp, nu*dt)
end do
print*, "todo seteado"
!iniciar simulacion
!abrimos archivo de datos
open (unit=0, file="p4_2_rho08.dat", status="replace", action="write")
open (unit=1, file="tray_rho08.xyz", status="replace", action="write")

!preparamos los parametros para la simulacion
ktotal = 200000
Temp = 1.2d0
rnew = r
vnew = v
t = 0.0d0
keq = ktotal - 10000 !momento en el que consideramos equilibrio, para calcular las medias
avEp = 0.0d0
avEk = 0.0d0
avp  = 0.0d0
write(unit=0, fmt=*) t, Ek(v), ELenardpbc(r), norm2(ptot(v))
do k = 1, ktotal
	call vVerletP(r, v, dt, rnew, vnew, p)
	r = rnew
	v = vnew
	do i = 1, N
		call pbc(r(i,:))
	end do
	call Ander(v, Temp, nu*dt)
	t = t + dt
	if ((k/100)*100 == k) then
		print*, real(k)/real(ktotal)*100, "% completado"
		!escribir datos de energias para grafico
		write(unit=0, fmt=*) t, Ek(v), ELenardpbc(r), norm2(ptot(v))
		!escribir datos para la visualización
		write(unit=1, fmt="(i0)") N
		write(unit=1, fmt=*)
		do i = 1, N
		write(unit=1, fmt="(a, f0.4, a, f0.4, a, f0.4)") "H ", r(i,1), " ", r(i,2), " ", r(i,3)
		end do
	end if
	if (k > keq) then
		avEp = (1.0d0/(k-keq)) * ((k-keq-1) * avEp + ELenardpbc(r))
		avEk = (1.0d0/(k-keq)) * ((k-keq-1) * avEk + Ek(v))
		avp  = (1.0d0/(k-keq)) * ((k-keq-1) * avp + p)
	end if
end do

open (unit=2, file="p4b_2_rho08.dat", status="replace", action="write")

avp = avp + kB * Temp * rho
write(unit=2, fmt=*) N, rho, avEp, avEk, avp

close(unit=0)
close(unit=1)
close(unit=2)

end program P4