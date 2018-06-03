module hk
    !! A module containing a Fortran implementation of the
    !! Hoshen-Kopelman algorithm, as well as the necessary
    !! union-find procedures.

    implicit none
    contains

        function hoshen_kopelman(matrix) result(num_clusters)
            !! Hoshen-Kopelman algorithm for labelling clusters on
            !! a binary matrix.
            !!
            !! The function takes in a binary matrix
            !! (zeros and positive numbers), labels the clusters of
            !! positive numbers and returns the total number of clusters
            !! found. The argument matrix is overwritten with the labels.
            integer, intent(inout), dimension(:,:) :: matrix
                !! The matrix to be labelled.
            integer :: num_clusters
                !! The total number of disjoint clusters labelled.
            integer :: m, n, i, j, up, left, label
            integer, dimension(:), allocatable :: labels, new_labels

            m = size(matrix, 1)
            n = size(matrix, 2)
            allocate(labels(m*n/2+1))
            labels = 0
            num_clusters = 0

            do j = 1, n
                do i = 1, m
                    if(matrix(i,j) > 0) then
                        if(i == 1) then
                            up = 0
                        else
                            up = matrix(i-1, j)
                        end if
                        if(j == 1) then
                            left = 0
                        else
                            left = matrix(i, j-1)
                        end if

                        ! New cluster
                        if(up==0 .and. left==0) then
                            num_clusters = num_clusters + 1
                            labels(num_clusters) = num_clusters
                            matrix(i,j) = num_clusters

                        ! Site binds clusters
                        else if(up > 0 .and. left > 0) then
                            matrix(i,j) = uf_union(up, left, labels)

                        ! Only one neighbour
                        else
                            matrix(i,j) = max(up, left)
                        end if
                    end if
                end do
            end do

            allocate(new_labels(m*n/2+1))
            new_labels = 0
            num_clusters = 0

            do j = 1, n
                do i = 1, m
                    if(matrix(i,j) > 0) then
                        label = uf_find(matrix(i,j), labels)
                        if(new_labels(label) == 0) then
                            num_clusters = num_clusters + 1
                            new_labels(label) = num_clusters
                        end if
                        matrix(i,j) = new_labels(label)
                    end if
                end do
            end do
        end function

        function uf_find(x, labels) result(y)
            !! Union-Find find algorithm:
            !! Find the lowest corresponding label.
            !! Relabelling is done when necessary.
            integer, intent(in) :: x
                !! Label for which to find the lowest corresponding label.
            integer, dimension(:), intent(inout) :: labels
                !! List of labels. **labels(i)** points to the lowest
                !! corresponding label of label **i**.
            integer :: y, z, tmp

            y = x
            do while(labels(y) /= y)
                y = labels(y)
            end do

            tmp = x
            do while(labels(tmp) /= tmp)
                z = labels(tmp)
                labels(tmp) = y
                tmp = z
            end do
        end function

        function uf_union(x, y, labels) result(canonical_label)
            !! Union-Find union algorithm:
            !! Merge two labels, return the result.
            integer, intent(in) :: x, y
                !! Labels to merge.
            integer, dimension(:), intent(inout) :: labels
                !! List of labels.
            integer :: canonical_label

            canonical_label = uf_find(y,labels)
            labels(uf_find(x,labels)) = canonical_label
        end function
end module hk
