module utilities
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    contains
        function stringfromint(x)
            character(len=:), allocatable :: stringfromint
            integer, intent(in) :: x
            integer :: numdigits

            numdigits = int(log10(real(x))) + 1
            allocate(character(len=numdigits) :: stringfromint)
            write(unit=stringfromint, fmt="(i0)") x
        end function

        function linspace(a,b,N)
            integer, intent(in) :: N
            real(kind=dp), intent(in) :: a,b
            real(kind=dp), dimension(:), allocatable :: linspace
            real(kind=dp) :: dx
            integer :: i

            dx = (b - a)/(N - 1)
            linspace = [ (a + i*dx, i=0, N-1) ]
        end function

        subroutine linfit(x, y, slope, const)
            real(kind=dp), dimension(:), intent(in) :: x, y
            real(kind=dp), intent(in out) :: slope, const
            real(kind=dp), dimension(:,:), allocatable :: A
            real(kind=dp), dimension(:), allocatable :: b, work
            integer :: lda, ldb, lwork, info, m

            m = size(x)
            allocate(A(m,2))
            allocate(b(m))
            allocate(work(1))

            A(:,1) = x(:)
            A(:,2) = 1
            b(:) = y(:)

            lda = m
            ldb = m

            lwork = -1
            call dgels("N",m,2,1,A,lda,b,ldb,work,lwork,info)
            lwork = int(work(1))
            deallocate(work)
            allocate(work(lwork))
            call dgels("N",m,2,1,A,lda,b,ldb,work,lwork,info)

            slope = b(1)
            const = b(2)

            if(info /= 0) then
                error stop "Problems with least squares"
            endif
        end subroutine
end module
