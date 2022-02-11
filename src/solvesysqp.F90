!
!-----------------------------------------------------------------------
subroutine solvesysqp ( m, n, A, G, x, y, lambda, rd, rp, rxs, deltax, deltay, deltalambda, info )
!-----------------------------------------------------------------------
implicit none
! inputs
integer,                          intent(in)  :: m
integer,                          intent(in)  :: n
double precision, dimension(m,n), intent(in)  :: A
double precision, dimension(n,n), intent(in)  :: G
double precision, dimension(n),   intent(in)  :: x
double precision, dimension(m),   intent(in)  :: y
double precision, dimension(m),   intent(in)  :: lambda
double precision, dimension(n),   intent(in)  :: rd
double precision, dimension(m),   intent(in)  :: rp
double precision, dimension(m),   intent(in)  :: rxs

! outputs
double precision, dimension(n),   intent(out) :: deltax
double precision, dimension(m),   intent(out) :: deltay
double precision, dimension(m),   intent(out) :: deltalambda
integer,                          intent(out) :: info

! local variables
integer i, j, nrhs, mn
double precision prec, cho
external dpotrs

! allocatable variables
double precision, dimension(:),   allocatable :: Adeltax
double precision, dimension(:,:), allocatable :: AD2
double precision, dimension(:,:), allocatable :: AD2A
double precision, dimension(:),   allocatable :: xsl
double precision, dimension(:),   allocatable :: cro
double precision, dimension(:,:), allocatable :: T1
double precision, dimension(:),   allocatable :: rd1
double precision, dimension(:),   allocatable :: rp1
double precision, dimension(:),   allocatable :: rxs1

parameter(prec = 1.0d-9, cho = 0)

! init
info = 0
mn = max(m,n)

! allocations
allocate(Adeltax(m))
allocate(AD2(m,n))
allocate(AD2A(n,n))
allocate(xsl(m))
allocate(cro(m))
allocate(T1(n,n))
allocate(rd1(n))
allocate(rp1(m))
allocate(rxs1(m))



#ifdef DEBUG
print*, "A=[[", ((A(i,j),j=1,n),char(10), i=1,m), "]]"
print*, "G=[[", ((G(i,j),j=1,n),char(10), i=1,n), "]]"
print*, "lambda=[", (lambda(i),i=1,m), "]"
print*, "y=[", (y(i),i=1,m), "]"
print*, "x=[", (x(i),i=1,n), "]"
print*, "rd=[", (rd(i),i=1,n), "]'"
print*, "rp=[", (rp(i),i=1,m), "]'"
print*, "rxs=[", (rxs(i),i=1,m), "]'"
#endif

! calculates (G+A'Y^-1 ΛA)
xsl = lambda / y   ! y^-1 lambda
do j=1,n
  AD2(:,j) = A(:,j)*xsl(j)
enddo
AD2A = matmul(transpose(A), AD2)
T1 = G + AD2A

#ifdef DEBUG
print*, "xsl=", (xsl(i), i=1,m)
print*, "AD2=[[", ((AD2(i, j),j=1,n),char(10), i=1,m), "]]"
print*, "AD2A=[[", ((AD2A(i,j),j=1,n),char(10), i=1,n), "]]"
print*, "T1=[",  ((T1(i,j),j=1,n),char(10), i=1,n), "]'"
#endif


! calculates deltax(val provisoire) = -rd+A'Y^{-1}Λ[-rp-rxs Λ^{-1}] =  -rd+A2'[-rp-rxs Λ^{-1}]
cro = -rp - rxs/lambda
deltax = matmul(transpose(AD2), cro) - rd

#ifdef DEBUG
print*, "termb=[",  (deltax(i),i=1,n), "]'"
#endif

! Solve (G+A'Y^-1 ΛA)deltax = -rd+A'Y^{-1}Λ[-rp-y+sigma mu Λ^{-1}e]
! cholesky factorization
if (cho == 0) then
  !  lapack cholesky
  call dpotrf('L', n, T1, n, info)
else ! cho==1
  ! modified cholesky
  call modchol2( m, T1, T1, info )
endif
if (info .NE. 0) then
  return
endif
#ifdef DEBUG
print*, "chol=", ((T1(i,j), j=1,n), char(10), i=1,n)
#endif

if (info .NE. 0) THEN
#ifdef DEBUG
  print*, "info1=", info
#endif
  return
endif
nrhs = 1
call dpotrs('U', n, nrhs, T1, n, deltax, n, info) ! dpotrs
if (info .NE. 0) THEN
#ifdef DEBUG
  print*, "info2=", info
#endif
  return
endif
#ifdef DEBUG
print*, "sol=[",  (deltax(i),i=1,n), "]'"
#endif


! calculates deltay = rp + AΔx
deltay = rp + matmul(A, deltax)

! calculates deltalambda = -Y^-1(rxs+ΛΔx)
deltalambda = -rxs/y -lambda*deltay/y

#ifdef DEBUG
print*, "deltax=[", (deltax(i),i=1,n), "]'"
print*, "deltalambda=[", (deltalambda(i),i=1,m), "]'"
print*, "deltay=[", (deltay(i),i=1,m), "]'"
#endif

!check the system of equations
rd1 = matmul(G, deltax) - matmul(transpose(A), deltalambda)
if (sum(abs(rd1+rd)).GT. prec) THEN
#ifdef DEBUG
  print*, "rd=", (rd(i), i=1,n)
  print*, "rd1=", (rd1(i), i=1,n)
  print*, "prob 1:", sum(abs(rd1+rd))
#endif
  info = -17
endif

rp1 = matmul(A, deltax) - deltay
if (sum(abs(rp1+rp)).GT. prec) THEN
#ifdef DEBUG
  print*, "rp=", (rp(i), i=1,m)
  print*, "rp1=", (rp1(i), i=1,m)
  print*, "prob 2:", sum(abs(rp1+rp))
#endif
  info = -18
endif

rxs1 = lambda*deltay+y*deltalambda
if (sum(abs(rxs1+rxs)).GT. prec) THEN
#ifdef DEBUG
  print*, "rxs=", (rxs(i), i=1,m)
  print*, "rxs1=", (rxs1(i), i=1,m)
  print*, "prob 3:", sum(abs(rxs1+rxs))
#endif
  info = -19
endif


RETURN
END


