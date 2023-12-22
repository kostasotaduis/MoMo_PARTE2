module modP4

real*8, public :: L, sigma, m
real*8, public :: rc
real*8, public, parameter :: kB = 1.0d0

public :: pbc, dpbc, latsc, Ander, vVerlet, FLenard, Ek, ELenardpbc

contains

subroutine pbc(r)

real*8, dimension(:), intent(inout) :: r
integer*4 :: i

do i = 1, 3
	if (r(i) >= L) then
		r(i) = r(i) - L
	else if (r(i) < 0.0d0) then
		r(i) = r(i) + L
	end if
end do

end subroutine pbc


function dpbc(r1, r2) result(distpbc)

real*8, dimension(:), intent(in) :: r1, r2
real*8 :: distpbc
real*8, dimension(3) :: r, zeroL
integer*4 :: i

zeroL = (/L, L, L/)

do i = 1, 3
	r(i) = (r2(i)-r1(i)+0.5*L) - L * sign(1.0d0, r2(i)-r1(i)+0.5*L) &
     * nint(abs((r2(i)-r1(i)+0.5*L)-0.5*L)/L) - 0.5*L
end do

distpbc = norm2(r)

end function dpbc

subroutine BoxMuller(sigma, n1, n2)

real*8, intent(in) :: sigma
real*8 :: x1, x2
real*8, intent(out) :: n1, n2

call random_number(x1)
call random_number(x2)

n1 = sigma * (-2.0 * log(1.0-x1))**0.5 * cos(2.0 * pi * x2)
n2 = sigma * (-2.0 * log(1.0-x2))**0.5 * sin(2.0 * pi * x1)

end subroutine BoxMuller

subroutine latsc(N, rho, sigma, r)

real*8, intent(in) :: rho, sigma
real*8, dimension(:,:), allocatable, intent(out) :: r
integer*4 :: i, j, k, N, Ln, m
real*8 :: a

Ln = nint(N**(1.0d0/3.0d0))

a = sigma * ((1.0d0 / rho) ** (1.0d0/3.0d0))
allocate(r(N,3))
m = 1
do i = 0, Ln-1
	do j = 0, Ln-1
		do k = 0, Ln-1
			r(m,1) = 0.0 + i * a
			r(m,2) = 0.0 + j * a
			r(m,3) = 0.0 + k * a
			m = m + 1
		end do
	end do
end do
			
end subroutine latsc

!Andersen thermostat

subroutine Ander(v, T)

real*8, dimension(:,:), intent(inout) :: v
real*8, intent(in) :: T
real*8 :: r, n1, n2
integer*4 :: i, j

call random_number(r)

i = nint(r*size(v,1) + 0.5)

do j = 1, 3
	call BoxMuller(((kB * T)**0.5), n1, n2)
	v(i,j) = n1
end do

end subroutine Ander

subroutine vVerlet(r, v, dt, rnew, vnew)

real*8, dimension(:,:), intent(in) :: v, r
real*8, dimension(:,:), allocatable, intent(out) :: vnew, rnew
real*8, intent(in) :: dt
real*8, dimension(3) :: f, fnew
integer*4 :: N, i, j
real*8 :: dist

N = size(r, 1)
allocate(rnew(N,3), vnew(N,3))

do i = 1, N
	f = 0.0d0
	do j = 1, N
		if (i /= j) then
			f = f + FLenard(r(i,:), r(j,:))
		end if
	end do
	rnew(i,:) = r(i,:) + v(i,:)*dt + 0.5d0 * f * dt**2.0d0
end do

do i = 1, N
	f = 0.0d0
	fnew = 0.0d0
	do j = 1, N
		if (i /= j) then
			f = f + FLenard(r(i,:), r(j,:))
			fnew = fnew + FLenard(rnew(i,:), rnew(j,:))
		end if
	end do
	vnew(i,:) = v(i,:) + 0.5d0 * (f + fnew) * dt
end do


end subroutine vVerlet

subroutine Euler(r, v, dt, rnew, vnew)

real*8, dimension(:,:), intent(in) :: v, r
real*8, dimension(:,:), allocatable, intent(out) :: vnew, rnew
real*8, intent(in) :: dt
real*8, dimension(3) :: f
integer*4 :: N, i, j

N = size(r, 1)
allocate(rnew(N,3), vnew(N,3))

do i = 1, N
	f = 0.0d0
	do j = 1, N
		if (i /= j) then
			f = f + FLenard(r(i,:), r(j,:))
		end if
	end do
	vnew(i,:) = v(i,:) + f * dt
	rnew(i,:) = r(i,:) + v(i,:) * dt
end do

end subroutine Euler

function FLenard(r1, r2) result(FL)

real*8, dimension(:), intent(in) :: r1, r2
real*8, dimension(3) :: rr, ur, FL
real*8 :: dist
integer*4 :: i


rr = r1 - r2



do i = 1, 3
	if (abs(rr(i)) > (L/2.0d0)) then
		rr(i) = rr(i) - L * (rr(i))/(abs(rr(i)))
	end if
end do

dist = norm2(rr)

!if (dist < sigma * 0.1) then
!	print*, "la liaste", dist
!end if

ur = rr / dist
FL = 0.0d0
if (dist <= rc) then
	FL = 4.0d0 * (12.0d0 * (sigma**12.0d0)/(dist**14.0d0) - &
	6.0d0 * (sigma**6.0d0)/(dist**8.0d0)) * rr
end if
!if (norm2(FL) > 10000.0) then
!	print*, "cagaste", dist
!end if

end function Flenard

function Ek(v) result(Ekin)

real*8, dimension(:,:), intent(in) :: v
real*8 :: Ekin
integer*4 :: i

Ekin = 0.0d0
do i = 1, size(v, 1)
	Ekin = Ekin + 0.5d0 * norm2(v(i,:))**2.0d0
end do

end function Ek

function ELenardpbc(r) result(ELen)

real*8, dimension(:,:), intent(in) :: r
real*8, dimension(3) :: rr
real*8 :: ELen
real*8 :: d2, Vij
integer*4 :: N, i, j, k

N = size(r,1)

ELen = 0.0d0

do i = 1, N
	do j = i+1, N
		if (j == i) then
			cycle
		end if
		rr = r(j,:)-r(i,:)
		do k = 1, 3
			if (abs(rr(k)) > (0.5d0*L)) then
				rr(k) = rr(k) - L * (rr(k))/(abs(rr(k)))
			end if
		end do
		d2 = (norm2(rr))**2.0
		if ((d2 <= rc**2.0d0)) then
			Vij = 4.0d0 * ((1.0d0/(d2**6.0d0)) - (1.0d0/(d2**3.0d0))) - &
			 4.0d0 * ((1.0d0/(rc**12.0d0)) - (1.0d0/(rc**6.0d0)))
			ELen = ELen + Vij
		end if
	end do
end do

end function ELenardpbc

function ptot(v) result(ptt)

real*8, dimension(:,:), intent(in) :: v
real*8, dimension(3) :: ptt
integer*4 :: i, N

N = size(v,1)

ptt = 0.0d0
do i = 1, N
	ptt = ptt + v(i,:)
end do

end function ptot

end module modP4