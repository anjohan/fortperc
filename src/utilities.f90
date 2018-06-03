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
            !! Create an array of **N** linearly spaced *values* (not intervals)
            !! from **a** to **b**.
            !! Similar to [`numpy.linspace(a, b, N)`](https://docs.scipy.org/doc/numpy/reference/generated/numpy.linspace.html).
            integer, intent(in) :: N
                !! Number of values (not intervals).
            real(kind=dp), intent(in) :: a
                !! Lower endpoint.
            real(kind=dp), intent(in) :: b
                !! Upper endpoint.
            real(kind=dp), dimension(:), allocatable :: linspace
                !! Array of linearly spaced values.
            real(kind=dp) :: dx
            integer :: i

            dx = (b - a)/(N - 1)
            linspace = [ (a + i*dx, i=0, N-1) ]
        end function

        subroutine linfit(x, y, slope, const)
            !! Compute a linear fit for the given data, return
            !! the slope and the constant term.
            !! `dgels` from LAPACK solves the linear least squares problem.
            real(kind=dp), dimension(:), intent(in) :: x, y
                !! Values to be fitted.
            real(kind=dp), intent(inout) :: slope
                !! \\(a\\) in \\(y=ax+b\\).
            real(kind=dp), intent(inout) :: const
                !! \\(b\\) in \\(y=ax+b\\).

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
            !! Find the common element in two arrays, given a total of
            !! **num_labels** unique elements. Used by [[find_spanning_cluster]].
            integer, dimension(:), intent(in) :: array1, array2
                !! Array to analyse.
            integer, intent(in) :: num_labels
                !! The known number of unique non-zero elements.
            integer :: intersect_label
                !! The first non-zero common element.
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
