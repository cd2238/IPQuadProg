! -----------------------------------------------------------------------*/
subroutine calcalpha(m, lambda, deltalambda, tau, alpha, info)
!-------------------------------------------------------------------------
implicit none
integer,                        intent(in)  :: m
double precision, dimension(m), intent(in)  :: lambda
double precision, dimension(m), intent(in)  :: deltalambda
double precision,               intent(in)  :: tau

double precision,               intent(out) :: alpha
integer,                        intent(out) :: info

! local variables
integer i
double precision epsi
parameter (epsi=1.0d-15)

info = 0
alpha = 0.99d0
do i=1, m
  if (deltalambda(i) .GT. epsi) then
    if ( - tau*lambda(i)/deltalambda(i)  .GT. 0.99d0 ) then
      info = -1
      return
    endif
  elseif(deltalambda(i) .LT. -epsi) then
    alpha = min(alpha, - tau*lambda(i)/deltalambda(i) )
  endif
enddo

alpha = max(0.0d0, alpha)



return
end subroutine
