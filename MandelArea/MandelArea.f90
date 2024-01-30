program mandelarea

implicit none

real, dimension(:), allocatable :: x !x-axis
real, dimension(:), allocatable :: inside    !notinmset counter

real    :: r, i !real & imaginary to create a random complex
complex :: z    !the random complex to be initialized
integer :: b, j ! loop stuff

integer :: currentiter = 0
integer, parameter :: startiter = 1000            !start at iteration startiter
integer, parameter :: iterstep = 0        !add(or mult by) iterstep at each batch
integer, parameter :: batch = 10               !number of point on x-axis
integer, parameter :: loopmax = 10000000          !montecarlo loop

ALLOCATE(x(batch)) 
ALLOCATE(inside(batch)) 

x = 0
inside = 0
currentiter =  startiter

DO b = 1, batch
    x(b) = b - 1
end do

DO b = 1, batch
    !x(b) = real(currentiter)
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(inside,b,currentiter) private(z,r,i,x)
    DO j = 1, loopmax
        call random_number(r)
        call random_number(i)    
        z = CMPLX(r*4.0-2.0, i*4.0-2.0)    

        if(isInMSet(z, currentiter)) then
            !$omp atomic update
            inside(b) = inside(b) + 1
        end if
    end do
    !$OMP END PARALLEL DO
    print *, b, currentiter, inside(b), (inside(b) / (loopmax - inside(b)) * 100.0)
    currentiter = currentiter + iterstep
end do

print *, "average : ", sum((inside / (loopmax - inside) * 100.0) / batch)


DEALLOCATE(x)
DEALLOCATE(inside)

CONTAINS 

pure function isInMSet(c, maxiter)

    complex, intent(in) :: c
    integer, intent(in) :: maxiter
    integer :: n
    complex :: z
    logical :: isInMSet
    
    z = CMPLX(0,0)
    n = 0
    IF(((ABS(c - CMPLX(-1,0))) < 0.25) .OR. ((ABS(1.0 - SQRT(1-(4*c)))) < 1.0)) THEN
        isInMSet = .TRUE.
        return  ! this if is kind of pointless if i forget the return ^^
    END IF

    DO WHILE(ABS(z) < 4 .AND. (n < maxiter))
        z = z*z+c
        n = n+1
    END do

    if(n >= maxiter) then
        isInMSet = .TRUE.
    else
        isInMset = .FALSE.
    end if    

end function isInMSet

end program mandelarea