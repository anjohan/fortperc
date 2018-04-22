module clusterlabelling
    use utilities
    implicit none
    contains
        function create_binary_matrix(p, L) result(binary_matrix)
            logical, dimension(:,:), allocatable :: binary_matrix
            real(kind=dp), intent(in) :: p
            integer, intent(in) :: L

            real(kind=dp), dimension(:,:), allocatable :: p_matrix
            allocate(p_matrix(L,L))

            call random_seed()
            call random_number(p_matrix)

            binary_matrix = p_matrix < p

        end function

        !/labelsubroutinestart/!
        subroutine label(matrix, labelled_matrix, number_of_labels)
            logical, dimension(:,:), allocatable, intent(in) :: matrix
            integer, dimension(:,:), allocatable, intent(inout) :: labelled_matrix
            integer, intent(inout) :: number_of_labels
            integer :: L, i, j

            L = size(matrix,1)
            number_of_labels = 0

            if(.not. allocated(labelled_matrix)) then
                allocate(labelled_matrix(L,L))
            end if

            labelled_matrix = 0

            do j=1,L
                do i=1,L
                    if(matrix(i,j) .and. labelled_matrix(i,j) == 0) then
                        number_of_labels = number_of_labels + 1
                        call growcluster(matrix,labelled_matrix,i,j,number_of_labels)
                    endif
                enddo
            enddo
        end subroutine
        !/labelsubroutineend/!

        !/growclustersubroutinestart/!
        recursive subroutine growcluster(matrix, labelled_matrix, i, j, label)
            logical, dimension(:,:), allocatable, intent(in) :: matrix
            integer, dimension(:,:), allocatable, intent(inout) :: labelled_matrix
            integer, intent(in) :: i, j, label
            integer :: L

            L = size(matrix,1)

            labelled_matrix(i,j) = label

            if(i<L) then
                if(matrix(i+1,j) .and. labelled_matrix(i+1,j)==0) then
                    call growcluster(matrix, labelled_matrix, i+1, j, label)
                endif
            endif
            if(j<L) then
                if(matrix(i,j+1) .and. labelled_matrix(i,j+1)==0) then
                    call growcluster(matrix, labelled_matrix, i, j+1, label)
                endif
            endif
            if(i>1) then
                if(matrix(i-1,j) .and. labelled_matrix(i-1,j)==0) then
                    call growcluster(matrix, labelled_matrix, i-1, j, label)
                endif
            endif
            if(j>1) then
                if(matrix(i,j-1) .and. labelled_matrix(i,j-1)==0) then
                    call growcluster(matrix, labelled_matrix, i, j-1, label)
                endif
            endif
        end subroutine
        !/growclustersubroutineend/!

        function find_spanning_cluster(labelled_matrix, number_of_labels) result(spanning_label)
            integer :: spanning_label
            integer, dimension(:,:), allocatable, intent(in) :: labelled_matrix
            integer, intent(in) :: number_of_labels
            integer :: L

            L = size(labelled_matrix,1)

            spanning_label = find_intersection(labelled_matrix(:,1), labelled_matrix(:,L), number_of_labels)

            if (spanning_label == -1) then
                spanning_label = find_intersection(labelled_matrix(1,:), labelled_matrix(L,:), number_of_labels)
            end if
        end function

        function find_intersection(array1, array2, number_of_labels) result(intersect_label)
            integer, dimension(:), intent(in) :: array1, array2
            integer :: intersect_label
            integer, intent(in) :: number_of_labels
            integer :: L, i
            logical, dimension(:), allocatable :: label_found

            L = size(array1)

            allocate(label_found(0:number_of_labels))
            !/intersectsnippetstart/!
            label_found(0:number_of_labels) = .false.

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

        function spanning_density_one_sample(p, L) result(spanning_density)
            real(kind=dp) :: spanning_density
            real(kind=dp), intent(in) :: p
            integer, intent(in) :: L

            logical, dimension(:,:), allocatable :: binary_matrix
            integer, dimension(:,:), allocatable :: labelled_matrix

            integer :: number_of_labels, spanning_label

            binary_matrix = create_binary_matrix(p,L)
            call label(binary_matrix, labelled_matrix, number_of_labels)

            spanning_label = find_spanning_cluster(labelled_matrix, number_of_labels)

            if (spanning_label == -1) then
                spanning_density = 0
                return
            end if

            spanning_density = count(labelled_matrix == spanning_label)/real(L**2,kind=dp)
        end function

        function spanning_density(p, L, number_of_samples)
            real(kind=dp) :: spanning_density
            real(kind=dp), intent(in) :: p
            integer, intent(in) :: L, number_of_samples

            integer :: i
            real(kind=dp), dimension(:), allocatable :: results

            allocate(results(number_of_samples))

            !$omp parallel do
            do i=1,number_of_samples
                results(i) = spanning_density_one_sample(p, L)
            end do
            !$omp end parallel do

            spanning_density = sum(results)/number_of_samples
        end function

        function has_spanning_cluster_one_sample(p, L) result(has_spanning)
            logical :: has_spanning
            real(kind=dp), intent(in) :: p
            integer, intent(in) :: L

            logical, dimension(:,:), allocatable :: binary_matrix
            integer, dimension(:,:), allocatable :: labelled_matrix

            integer :: number_of_labels, spanning_label

            binary_matrix = create_binary_matrix(p,L)
            call label(binary_matrix, labelled_matrix, number_of_labels)

            spanning_label = find_spanning_cluster(labelled_matrix, number_of_labels)

            has_spanning = spanning_label /= -1
        end function

        function spanning_probability(p, L, number_of_samples)
            real(kind=dp) :: spanning_probability
            real(kind=dp), intent(in) :: p
            integer, intent(in) :: L, number_of_samples

            integer :: i
            logical, dimension(:), allocatable :: results

            allocate(results(number_of_samples))

            !$omp parallel do
            do i=1,number_of_samples
                results(i) = has_spanning_cluster_one_sample(p, L)
            end do
            !$omp end parallel do

            spanning_probability = count(results)/real(number_of_samples,kind=dp)
        end function
end module clusterlabelling
