! A is supposed to be row full rank (implies m<=n)
!-----------------------------------------------------------------------
subroutine initqp ( n, G, c, m, A, b, x, y, lambda, info )
!-----------------------------------------------------------------------
implicit none
! inputs
integer,                          intent(in)  :: n
double precision, dimension(n,n), intent(in)  :: G
double precision, dimension(n),   intent(in)  :: c
integer,                          intent(in)  :: m
double precision, dimension(m,n), intent(in)  :: A
double precision, dimension(m),   intent(in)  :: b


! outputs
double precision, dimension(n),   intent(out) :: x
double precision, dimension(m),   intent(out) :: y
double precision, dimension(m),   intent(out) :: lambda
integer,                          intent(out) :: info

! local variable
integer i
external solvesysqp

! allocatable
double precision, dimension(:), allocatable :: rd
double precision, dimension(:), allocatable :: rp
double precision, dimension(:), allocatable :: rxs
double precision, dimension(:), allocatable :: deltax
double precision, dimension(:), allocatable :: deltay
double precision, dimension(:), allocatable :: deltalambda


! allocation
allocate(rd(n))
allocate(rp(m))
allocate(rxs(m))
allocate(deltax(n))
allocate(deltay(m))
allocate(deltalambda(m))

!init
info = 0
x      = 1.0d0/dble(n)
y      = 1.0d0
lambda = 1.0d0


!rd
rd = matmul(G, x) + c - matmul(transpose(A), lambda)

!rp
rp = matmul(A, x) - y - b

!rxs
rxs = - y * lambda

call solvesysqp ( m, n, A, G, x, y, lambda, rd, rp, rxs, &
                        deltax, deltay, deltalambda, info )

if (info .NE. 0) then
  return
endif

#ifdef DEBUG
print*, "init"
print*, "yy", (y(i), i=1,m)
print*, "lamb", (lambda(i), i=1,m)
#endif

y = max(1.0d0, abs(y+deltay))
lambda =max(1.0d0, abs(lambda+deltalambda))

#ifdef DEBUG
print*, "sort"
print*, "yy", (y(i), i=1,m)
print*, "lamb", (lambda(i), i=1,m)
#endif

return
end


