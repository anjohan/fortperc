module clusterlabeling
    implicit none
    contains
        subroutine label(matrix, labeled_matrix, number_of_labels)
            logical, dimension(:,:), allocatable, intent(in) :: matrix
            integer, dimension(:,:), allocatable, intent(inout) :: labeled_matrix
            integer, intent(inout) :: number_of_labels
            integer :: L, i, j

            L = size(matrix,1)
            number_of_labels = 0

            do j=1,L
                do i=1,L
                    ! print *, "i = ", i, "j = ", j
                    if(matrix(i,j) .and. labeled_matrix(i,j) == 0) then
                        number_of_labels = number_of_labels + 1
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

end module clusterlabeling
