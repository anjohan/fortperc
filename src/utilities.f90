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

        subroutine find_random_point(matrix, i, j)
            !! Find a random position on **matrix** such that
            !! **matrix(i, j)** is `.true.`.

            logical, dimension(:,:), intent(in) :: matrix
                !! Matrix whose `.true.` values are allowed positions for the
                !! random walker.

            integer, intent(out) :: i, j
                !! Returned random point on **matrix**.

            integer :: L
            real(kind=dp) :: x0_real, y0_real

            L = size(matrix, 1)

            ! Find random starting point.
            do
                call random_number(x0_real)
                call random_number(y0_real)

                x0_real = L*x0_real + 1
                y0_real = L*y0_real + 1

                i = int(x0_real)
                j = int(y0_real)

                if(matrix(i,j)) then
                    return
                end if
            end do
        end subroutine

        subroutine periodic_wraparound(x, L)
            integer, intent(inout) :: x
            integer, intent(in) :: L

            if(x == 0) then
                !write(*,*) "before:", x
                x = L
                !write(*,*) "after:", x
            else if(x == L+1) then
                !write(*,*) "before:", x
                x = 1
                !write(*,*) "after:", x
            end if
        end subroutine

        function mark_percolating_with_periodic(label_matrix, &
                                                percolating_label, &
                                                num_clusters) &
                                                result(matrix)
            !! Return a `logical` matrix where the `.true.` values are the sites
            !! belonging to the percolating cluster, when periodic boundary
            !! conditions are considered.

            integer, dimension(:,:), intent(in) :: label_matrix
                !! Labelled matrix from [[label]] or [[hoshen_kopelman]].
            integer, intent(in) :: percolating_label
                !! Known label of the percolating cluster.
            integer, intent(in) :: num_clusters
                !! The known number of clusters.
            logical, dimension(:,:), allocatable :: matrix
                !! Matrix where the `.true.` values are the sites belonging to
                !! the percolating cluster with periodic boundary conditions.

            integer :: L, i, j
            logical, dimension(:), allocatable :: connected_to_percolating

            L = size(label_matrix, 1)

            allocate(matrix(L,L))
            allocate(connected_to_percolating(0:num_clusters))

            connected_to_percolating(:) = .false.
            connected_to_percolating(percolating_label) = .true.

            do i = 1, L
                if(label_matrix(i,1) == percolating_label) then
                    connected_to_percolating(label_matrix(i, L)) = .true.
                end if

                if(label_matrix(i,L) == percolating_label) then
                    connected_to_percolating(label_matrix(i, 1)) = .true.
                end if

                if(label_matrix(1,i) == percolating_label) then
                    connected_to_percolating(label_matrix(L, i)) = .true.
                end if

                if(label_matrix(L, i) == percolating_label) then
                    connected_to_percolating(label_matrix(1, i)) = .true.
                end if
            end do

            connected_to_percolating(0) = .false.

            !$omp parallel do
            do j = 1, L
                do i = 1, L
                    matrix(i,j) = connected_to_percolating(label_matrix(i,j))
                end do
            end do
            !$omp end parallel do
        end function
end module
