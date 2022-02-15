!> Main program in order to execute the quadratic program (quadprosimp)
!> Example of the command line : ./quadpro n10c04
!> the maximal number of iterations itermax, and the duality measure tolerance mutol has to be modified in the code below.
!-----------------------------------------------------------------------
program main

     implicit none
     integer :: num_args, allocock, i, isx0
     character(len=12), dimension(:), allocatable   :: args
     integer                                        :: nb_control_variables
     double precision, dimension(:,:), allocatable  :: quadratic_objective
     double precision, dimension(:), allocatable    :: linear_objective
     integer                                        :: nb_inequality_constraint
     double precision, dimension(:,:), allocatable  :: constraint_matrix
     double precision, dimension(:), allocatable    :: constraint_vector
     integer                                        :: itermax
     double precision                               :: mutol
     double precision, dimension(:), allocatable    :: x0


     double precision, dimension(:), allocatable    :: x
     double precision, dimension(:), allocatable    :: xref

     double precision, dimension(:), allocatable    :: lambda   
     double precision, dimension(:), allocatable    :: y
     double precision                               :: fobj
     double precision                               :: fref

     integer                                        :: iter
     integer                                        :: info

     integer                                        :: iosxref
     integer                                        :: iosfref


     ! parameters
     parameter(itermax = 1000, mutol=1.0d-10)

     ! init
     iosfref = 0
     iosxref = 0

     ! read data
     num_args = command_argument_count()
     allocate(args(num_args))  ! 
     if (num_args >= 1) then
         call get_command_argument(1,args(1))
         print*, "dataset=", args(1)
     else
         print*, "1 argument is needed (name of the data folder)."
         return
     end if
     open(12, file="../data/"//trim(args(1))//"/nb_control_variables.txt")
     read(12,*) nb_control_variables
     !print*, nb_control_variables
     close(12)
     open(12, file="../data/"//trim(args(1))//"/nb_inequality_constraint.txt")
     read(12,*) nb_inequality_constraint
     !print*, nb_inequality_constraint
     close(12) 
     
     ! allocate arrays
     allocate (quadratic_objective(nb_control_variables, nb_control_variables), stat = allocock)
     if (allocock /= 0) return   
     allocate (linear_objective(nb_control_variables), stat = allocock)
     if (allocock /= 0) return
     allocate (constraint_matrix(nb_inequality_constraint, nb_control_variables), stat = allocock)
     if (allocock /= 0) return    
     allocate (constraint_vector(nb_inequality_constraint), stat = allocock)
     if (allocock /= 0) return
     allocate (xref(nb_control_variables), stat = allocock)
     if (allocock /= 0) return

     open(12, file="../data/"//trim(args(1))//"/quadratic_objective.txt")
     read(12,*) quadratic_objective
     !print*, linear_objective
     close(12)
     open(12, file="../data/"//trim(args(1))//"/linear_objective.txt")
     read(12,*) linear_objective
     !print*, linear_objective
     close(12)
     open(12, file="../data/"//trim(args(1))//"/constraint_matrix.txt")
     read(12,*) constraint_matrix
     close(12) 
     open(12, file="../data/"//trim(args(1))//"/constraint_vector.txt")
     read(12,*) constraint_vector
     !print*, constraint_vector
     close(12) 
     open(12, file="../data/"//trim(args(1))//"/xref.txt")
     read(12,*, iostat= iosxref) xref
     close(12)
     open(12, file="../data/"//trim(args(1))//"/fref.txt")
     read(12,*, iostat= iosfref) fref
     close(12)

     ! outputs
     allocate (x(nb_control_variables), stat = allocock)
     if (allocock /= 0) return 
     allocate (lambda(nb_inequality_constraint), stat = allocock)
     if (allocock /= 0) return 
     allocate (y(nb_inequality_constraint), stat = allocock)
     if (allocock /= 0) return


     
     allocate(x0(nb_control_variables))
     isx0 = 0
     x0 = 0.0d0
 
     ! perform optimization
     call quadprosimp ( nb_control_variables, quadratic_objective, linear_objective, &
                       nb_inequality_constraint, constraint_matrix, constraint_vector, &
                       isx0, x0, &
                       itermax, mutol, &
                       x, y, lambda, fobj, iter, info )

     
#ifdef DEBUG
     ! write results
     print*, (x(i), i=1, nb_control_variables)
     print*, (lambda(i), i=1, nb_inequality_constraint)     
     print*, (y(i), i=1, nb_inequality_constraint)
#endif
     print*, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
     print*, "Results :"
     print*, "xobj =", (x(i), i=1,nb_control_variables)
     if (iosxref == 0) then
       print*, "xref =", (xref(i), i=1,nb_control_variables)
     endif
     print*, "fobj =", fobj
     if (iosfref == 0) then
       print*, "fref =", fref
     endif
     print*, "constraints violated ? (violation if > 0)"
     do i = 1, nb_inequality_constraint
        print*, "const ", i, ":", constraint_vector(i) - dot_product(constraint_matrix(i,:), x)
     enddo
     print*, "iter/itermax :", iter, "/", itermax
     print*, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

     deallocate(x)
     deallocate(xref)
     deallocate(y)
     deallocate(x0)
     deallocate(quadratic_objective)
     deallocate(lambda)
     deallocate(args)
     deallocate(linear_objective)
     deallocate(constraint_matrix)
     deallocate(constraint_vector)

end program
