module utilities
    !! Useful procedures when doing percolation.
    implicit none
    integer, parameter :: dp = kind(1.0d0)
        !! Kind used for all **real** variables.
    contains
        function stringfromint(x)
            !! Make a string of "correct" length from a positive integer.
            character(len=:), allocatable :: stringfromint
                !! String containing the given integer, without spaces.
            integer, intent(in) :: x
                !! Positive integer to be converted.
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

        function find_intersection(array1, array2, num_labels) result(intersect_label)
            integer, dimension(:), intent(in) :: array1, array2
            integer :: intersect_label
            integer, intent(in) :: num_labels
            integer :: L, i
            logical, dimension(:), allocatable :: label_found

            L = size(array1)

            allocate(label_found(0:num_labels))
            !/intersectsnippetstart/!
            label_found(0:num_labels) = .false.

            do i=1,L
                label_found(array1(i)) = .true.
            end do

            do i=1,L
                if(array2(i) /= 0 .and. label_found(array2(i))) then
                    intersect_label = array2(i)
                    return
                end if
            end do
            !/intersectsnippetend/!

            intersect_label = -1
        end function
end module
