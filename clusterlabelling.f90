module clusterlabelling
    use utilities
    implicit none
    real(kind=dp), parameter :: pc = 0.592746
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

        function sizes(labelled_matrix, number_of_labels)
            integer, dimension(:), allocatable :: sizes
            integer, dimension(:,:), intent(in) :: labelled_matrix
            integer, intent(in) :: number_of_labels
            integer :: L,i,j

            L = size(labelled_matrix,1)
            allocate(sizes(number_of_labels))
            sizes(:) = 0

            do j=1,L
                do i=1,L
                    associate (label => labelled_matrix(i,j))
                        if(label /= 0) then
                            sizes(label) = sizes(label) + 1
                        endif
                    end associate
                enddo
            enddo
        end function

        subroutine cluster_number_density(p, L, num_samples, final_bin_mids, results)
            integer, intent(in) :: L, num_samples
            real(kind=dp), intent(in) :: p
            real(kind=dp), dimension(:), intent(inout), allocatable :: final_bin_mids, results

            integer :: num_bins, i, j, num_labels, sizeindex, maxindex, spanning_label
            real(kind=dp) :: a, loga

            logical, dimension(:,:), allocatable :: binary_matrix
            integer, dimension(:,:), allocatable :: label_matrix
            integer, dimension(:), allocatable :: clustersizes, histogram
            real(kind=dp), dimension(:), allocatable :: bin_edges, bin_sizes, bin_mids

            a = 1.5d0
            loga = log(a)

            num_bins = ceiling(log(1.0d0*L*L)/loga)
            bin_edges = a**[(i, i=0,num_bins)]
            bin_mids = 0.5*(bin_edges(1:num_bins) + bin_edges(2:num_bins+1))
            bin_sizes = bin_edges(2:num_bins+1) - bin_edges(1:num_bins)

            allocate(histogram(1:num_bins))

            do i = 1, num_samples
                binary_matrix = create_binary_matrix(p, L)
                call label(binary_matrix, label_matrix, num_labels)
                clustersizes = find_sizes(label_matrix, num_labels)
                spanning_label = find_spanning_cluster(label_matrix, num_labels)

                do j = 1, num_labels
                    if(j /= spanning_label) then
                        sizeindex = floor(log(1.0d0*clustersizes(j))/loga) + 1
                        histogram(sizeindex) = histogram(sizeindex) + 1
                    end if
                end do
            end do

            do maxindex = num_bins, 1, -1
                if(histogram(maxindex) /= 0) exit
            end do
            final_bin_mids = bin_mids(:maxindex)
            results = histogram(:maxindex)/(bin_sizes(:maxindex)*num_samples)


        end subroutine
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

        subroutine label_wrongly(matrix, labelled_matrix, number_of_labels)
            logical, dimension(:,:), allocatable, intent(in) :: matrix
            integer, dimension(:,:), allocatable, intent(inout) :: labelled_matrix
            integer, intent(inout) :: number_of_labels
            integer, dimension(:), allocatable :: label_map, reduced_label_map
            integer :: L, i, j, highest_label
            integer :: up, left


            L = size(matrix,1)
            number_of_labels = 0

            if(.not. allocated(labelled_matrix)) then
                allocate(labelled_matrix(L,L))
            end if

            labelled_matrix = 0
            highest_label = 0

            label_map = [ (0, i=1,L*L) ]
            ! write(*,*) "Label map:", label_map

            do j=1,L
                do i=1,L
                    if(matrix(i,j)) then
                        if(i==1) then
                            up = 0
                        else
                            up = labelled_matrix(i-1, j)
                        endif
                        if(j==1) then
                            left = 0
                        else
                            left = labelled_matrix(i, j-1)
                        endif

                        if(left==0) then
                            if(up==0) then
                                highest_label = highest_label + 1
                                labelled_matrix(i,j) = highest_label
                                label_map(highest_label) = highest_label
                            else
                                labelled_matrix(i,j) = up
                            endif
                        else
                            if(up==0 .or. left==up) then
                                labelled_matrix(i,j) = left
                            else
                                labelled_matrix(i,j) = min(left,up)
                                label_map(max(left,up)) = min(left,up)
                            endif
                        endif
                    endif
                enddo
            enddo

            ! write(*,*) "Labels before reduction:"
            ! do i = 1, highest_label
            !     write(*,*) i, label_map(i)
            ! enddo

            do i = highest_label, 1, -1
                j = i
                do while(label_map(j) /= j)
                    j = label_map(j)
                end do
                label_map(i) = j
            enddo

            ! write(*,*) "Labels after reduction:"
            ! do i = 1, highest_label
            !     write(*,*) i, label_map(i)
            ! enddo

            allocate(reduced_label_map(highest_label))
            reduced_label_map = 0
            number_of_labels = 0
            do i = 1, highest_label
                if(label_map(i) == i) then
                    number_of_labels = number_of_labels + 1
                    reduced_label_map(i) = number_of_labels
                endif
            enddo

            ! write(*,*) "Numer of  labels:", number_of_labels

            ! write(*,*) "Overwriting labels."

            do concurrent(j=1:L, i=1:L, labelled_matrix(i,j)/=0)
                labelled_matrix(i,j) = reduced_label_map(label_map(labelled_matrix(i,j)))
            enddo
        end subroutine
end module clusterlabelling
