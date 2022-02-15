!> Quadratic programming
!> from Nocedal Wright
!> @param n the number of control variables
!> @param G the matrix objective (n,n)
!> @param c the vector objective (n)
!> @param m the number of linear constraints
!> @param A the matrix constraint (m,n)
!> @param b the vector constraint (m)
!> @param isx0 1 if an primal initial value is provided, 0 if not.
!> @param x0 the primal initial value (n)
!> @param itermax maximal number of iterations
!> @param mutol value of mu (duality parameter) stopping criterion.
!> @param x the primal solution (n)
!> @param y the slack solution (m)
!> @param lambda the dual solution (m)
!> @param fobj the objective solution
!> @param iter the number of iterations realized
!> @param info information code
!-----------------------------------------------------------------------
subroutine quadprosimp(n, G, c, m, A, b, isx0, x0, itermax, mutol, x, y, lambda, fobj, iter, info)

implicit none
! inputs
integer,                          intent(in)  :: n
double precision, dimension(n,n), intent(in)  :: G
double precision, dimension(n),   intent(in)  :: c
integer,                          intent(in)  :: m
double precision, dimension(m,n), intent(in)  :: A
double precision, dimension(m),   intent(in)  :: b
integer,                          intent(in)  :: isx0
double precision, dimension(n),   intent(in)  :: x0
integer,                          intent(in)  :: itermax
double precision,                 intent(in)  :: mutol


! outputs
double precision, dimension(n),   intent(out) :: x
double precision, dimension(m),   intent(out) :: y
double precision, dimension(m),   intent(out) :: lambda
double precision,                 intent(out) :: fobj

integer,                          intent(out) :: iter
integer,                          intent(out) :: info


! locall variables
double precision alpha_aff, alpha_aff1, alpha_aff2, alpha, alpha1, alpha2, sigma
double precision ps, mu_aff, tau, mu
#ifdef DEBUG
integer i,j
#endif

EXTERNAL initqp, solvesysqp, calcalpha


! allocatable
double precision, dimension(:), allocatable :: rd
double precision, dimension(:), allocatable :: rp
double precision, dimension(:), allocatable :: rxs
double precision, dimension(:), allocatable :: deltax
double precision, dimension(:), allocatable :: deltay
double precision, dimension(:), allocatable :: deltalambda
double precision, dimension(:), allocatable :: yaff
double precision, dimension(:), allocatable :: lambdaaff


! init
info = 0
iter = 0

#ifdef DEBUG
print*, "G=[[", ((G(i,j),j=1,n),char(10), i=1,n), "]]"
print*, "c=[", (c(i),i=1,n), "]"
print*, "A=[[", ((A(i,j),j=1,n),char(10), i=1,m), "]]"
print*, "b=[", (b(i),i=1,m), "]"
#endif


! allocation
allocate(rd(n))
allocate(rp(m))
allocate(rxs(m))
allocate(deltax(n))
allocate(deltay(m))
allocate(deltalambda(m))
allocate(yaff(m))
allocate(lambdaaff(m))


! init
call initqp ( n, G, c, m, A, b, isx0, x0, x, y, lambda, info )
if (info /= 0) then
  return
endif
#ifdef DEBUG
  print*, "resp cont=", matmul(A,x)-b
#endif


#ifdef DEBUG
  print*, "x0=", (x(i), i=1,n)
  print*, "y0=", (y(i), i=1,m)
  print*, "lambda0=", (lambda(i), i=1,m)
#endif

! duality measure
ps = dot_product(y, lambda)
mu = ps/dble(m)

! algorithm
iter = 0
do while ((mu > mutol).AND.(iter < itermax))
    iter = iter + 1

    !rd
    rd = matmul(G, x) + c - matmul(transpose(A), lambda)

    !rp
    rp = matmul(A, x) - y - b

    !rxs
    rxs = y * lambda

    ! affine scaling step
    call solvesysqp ( m, n, A, G, y, lambda, rd, rp, rxs, &
                         deltax, deltay, deltalambda, info )
    if (info /= 0) then
      return
    endif

!   duality measure
    ps = dot_product(lambda, y)
    mu = ps/dble(m)
#ifdef DEBUG
    print*, "mu=", mu
#endif


! length step
    call calcalpha(m, lambda, deltalambda, 1.0d0, alpha_aff1, info)
    if (info < 0) then
      return
    endif
    call calcalpha(m, y, deltay, 1.0d0, alpha_aff2, info)
    if (info < 0) then
      return
    endif
    alpha_aff = min(alpha_aff1, alpha_aff2)

#ifdef DEBUG
    print*, "alpha_aff=", alpha_aff
#endif


    yaff = y+alpha_aff*deltay
    lambdaaff = lambda+alpha_aff*deltalambda
    ps = dot_product(yaff, lambdaaff)
    mu_aff = ps / dble(m)
#ifdef DEBUG
    print*, "mu_aff=", mu_aff
#endif

    ! centering
    sigma = (mu_aff/mu)**3
#ifdef DEBUG
    print*, "sigma=", sigma
#endif

    ! corrector/centering step
    rxs = rxs + deltalambda*deltay - sigma*mu
    call solvesysqp ( m, n, A, G, y, lambda, rd, rp, rxs, deltax, deltay, deltalambda, info )
    if (info /= 0) then
      return
    endif

    ! step length
    tau = 1.0d0 - sqrt(0.1d0/dble(iter)) !0.9d0
    call calcalpha(m, lambda, deltalambda, tau, alpha1, info)
    if (info < 0) then
      return
    endif
    call calcalpha(m, y, deltay, tau, alpha2, info)
    if (info < 0) then
      return
    endif
    alpha = min(alpha1, alpha2)

    ! update
    x = x + alpha*deltax
    lambda = lambda + alpha*deltalambda
    y  = y + alpha*deltay
    fobj = dot_product(c,x) + 0.5*dot_product(x, matmul(G,x))

#ifdef DEBUG
    print*, "~~~~"
    print*, "alpha=", alpha
    print*, "x=", (x(i), i=1,n)
    print*, "lambda=", (lambda(i), i=1,m)
    print*, "y=", (y(i), i=1,m)
    print*, "mu=", mu
    print*, "~~~~"
#endif



enddo

return
end subroutine
