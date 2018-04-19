module clusterlabeling
    implicit none
    contains
        function create_binary_matrix(p, L) result(binary_matrix)
            logical, dimension(:,:), allocatable :: binary_matrix
            real, intent(in) :: p
            integer, intent(in) :: L

            real, dimension(:,:), allocatable :: p_matrix
            allocate(p_matrix(L,L))

            call random_seed()
            call random_number(p_matrix)

            binary_matrix = p_matrix < p

        end function

        subroutine label(matrix, labeled_matrix, number_of_labels)
            logical, dimension(:,:), allocatable, intent(in) :: matrix
            integer, dimension(:,:), allocatable, intent(inout) :: labeled_matrix
            integer, intent(inout) :: number_of_labels
            integer :: L, i, j

            L = size(matrix,1)
            number_of_labels = 0

            if(.not. allocated(labeled_matrix)) then
                allocate(labeled_matrix(L,L))
            end if

            labeled_matrix = 0

            do j=1,L
                do i=1,L
                    print *, "i = ", i, "j = ", j
                    if(matrix(i,j) .and. labeled_matrix(i,j) == 0) then
                        number_of_labels = number_of_labels + 1
                        print *,number_of_labels
                        call growcluster(matrix,labeled_matrix,i,j,number_of_labels)
                    endif
                enddo
            enddo
        end subroutine

        recursive subroutine growcluster(matrix, labeled_matrix, i, j, label)
            logical, dimension(:,:), allocatable, intent(in) :: matrix
            integer, dimension(:,:), allocatable, intent(inout) :: labeled_matrix
            integer, intent(in) :: i, j, label
            integer :: L

            L = size(matrix,1)

            labeled_matrix(i,j) = label

            if(i<L .and. matrix(i+1,j) .and. labeled_matrix(i+1,j)==0) then
                call growcluster(matrix, labeled_matrix, i+1, j, label)
            endif
            if(j<L .and. matrix(i,j+1) .and. labeled_matrix(i,j+1)==0) then
                call growcluster(matrix, labeled_matrix, i, j+1, label)
            endif
            if(i>1 .and. matrix(i-1,j) .and. labeled_matrix(i-1,j)==0) then
                call growcluster(matrix, labeled_matrix, i-1, j, label)
            endif
            if(j>1 .and. matrix(i,j-1) .and. labeled_matrix(i,j-1)==0) then
                call growcluster(matrix, labeled_matrix, i, j-1, label)
            endif
       end subroutine

       function find_spanning_cluster(labeled_matrix, number_of_labels) result(spanning_label)
           integer :: spanning_label
           integer, dimension(:,:), allocatable, intent(in) :: labeled_matrix
           integer, intent(in) :: number_of_labels
           integer :: L

           L = size(labeled_matrix,1)

           spanning_label = find_intersection(labeled_matrix(:,1), labeled_matrix(:,L), number_of_labels)

           if (spanning_label == -1) then
               spanning_label = find_intersection(labeled_matrix(1,:), labeled_matrix(L,:), number_of_labels)
           end if
       end function

       function find_intersection(array1, array2, number_of_labels) result(intersect_label)
           integer, dimension(:), intent(in) :: array1, array2
           integer :: intersect_label
           integer, intent(in) :: number_of_labels
           integer :: L, i
           logical, dimension(:), allocatable :: label_found

           L = size(array1)

           allocate(label_found(number_of_labels))
           label_found = .false.

           do i=1,L
               if(array1(i) /= 0) then
                   label_found(array1(i)) = .true.
               end if
           end do

           do i=1,L
               if(array2(i) /= 0 .and. label_found(array2(i))) then
                   intersect_label = array2(i)
                   return
                end if
            end do

            intersect_label = -1

        end function

end module clusterlabeling
