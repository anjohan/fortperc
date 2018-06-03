module percolation
    !! A module with procedures for calculating various properties which
    !! are interesting when doing numerical percolation experiments.
    use utilities
    use hk
    implicit none

    private :: has_spanning_cluster_one_sample, spanning_density_one_sample

    real(kind=dp), parameter :: pc = 0.592746
        !! Known critical probability for a two-dimensional site-percolating
        !! system.
    contains
        function create_binary_matrix(p, L) result(binary_matrix)
            !! Create a random, binary (`logical`) matrix,
            !! which can be used in percolation experiments.

            logical, dimension(:,:), allocatable :: binary_matrix
                !! Randomly created matrix, where each element is `.true.`
                !! or `.false.` if a randomly generated number is smaller
                !! or greater than p.
            real(kind=dp), intent(in) :: p
                !! Probability for each matrix element to be `.true.`.
            integer, intent(in) :: L

            real(kind=dp), dimension(:,:), allocatable :: p_matrix
            allocate(p_matrix(L,L))

            call random_seed()
            call random_number(p_matrix)

            binary_matrix = p_matrix < p
        end function

        subroutine label(matrix, labelled_matrix, num_labels)
            !! Alternative interface to the Hoshen-Kopelman algorithm
            !! from the hk module, which uses a binary matrix created
            !! by [[create_binary_matrix]].
            logical, dimension(:,:), intent(in) :: matrix
                !! Binary matrix where clusters will be identified.
            integer, dimension(:,:), allocatable, intent(inout) :: labelled_matrix
                !! Integer matrix which will store the labels of each cluster.
                !! Reallocated if necessary.
            integer, intent(out) :: num_labels
                !! Overwritten with the total number of disjoint clusters.
            integer :: L

            L = size(matrix,1)

            if(.not. allocated(labelled_matrix)) then
                allocate(labelled_matrix(L,L))
            else if(size(labelled_matrix,1) /= L) then
                deallocate(labelled_matrix)
                allocate(labelled_matrix(L,L))
            endif
            where(matrix)
                labelled_matrix = 1
            else where
                labelled_matrix = 0
            end where
            num_labels = hoshen_kopelman(labelled_matrix)
        end subroutine

        function find_sizes(labelled_matrix, num_labels) result(sizes)
            !! Count the number of sites belonging to each cluster
            !! in the labelled matrix.
            integer, dimension(:), allocatable :: sizes
                !! Array of cluster sizes. The i'th element is the size of
                !! the cluster with label i.
            integer, dimension(:,:), intent(in) :: labelled_matrix
                !! Integer matrix with labelled clusters, resulting from
                !! the Hoshen-Kopelman algorithm.
            integer, intent(in) :: num_labels
                !! The known number of clusters.
            integer :: L,i,j

            L = size(labelled_matrix,1)
            allocate(sizes(num_labels))
            sizes(:) = 0

            do j=1,L
                do i=1,L
                    associate (label => labelled_matrix(i,j))
                        if(label /= 0) then
                            sizes(label) = sizes(label) + 1
                        end if
                    end associate
                end do
            end do
        end function

        subroutine cluster_number_density(p, L, num_samples, bin_mids, results, binsize_base)
            !! Calculate the number density of clusters with different sizes.
            !! The cluster number density is defined as
            !! \\begin{equation}
            !!      n(s,p) = \sum_{\text{MC-samples}}
            !!                   \frac{\text{number of clusters with size }s}
            !!                        {L^2\cdot\text{number of MC-samples}}.
            !! \\end{equation}
            !! Direct calculations will usually give bad results, as there
            !! will be very few
            !! clusters with large sizes compared to the numbers of clusters
            !! with small sizes. This is circumvented by doing logarithmic
            !! binning and averaging, i.e.
            !! \\begin{equation}
            !!      n\left([s+\Delta s),p\right) = \sum_{\text{MC-samples}}
            !!                   \frac{\text{number of clusters with size }s
            !!                         \in[s,s+\Delta s)}
            !!                        {\Delta s\cdot L^2\cdot\text{number of MC-samples}},
            !!      \label{eq:nsp}
            !! \\end{equation}
            !! where \\(\Delta s\\) are logarithmically distributed. After execution,
            !! **bin_mids** will contain the centres of the bins.

            integer, intent(in) :: L
                !! Percolating systems will be \\(L\times L\\).
            integer, intent(in) :: num_samples
                !! Results will be averaged over this number of Monte Carlo-samples.
                !! Sampling is parallelised.
            real(kind=dp), intent(in) :: p
                !! Probability for a given site to allow transport.
            real(kind=dp), intent(in), optional :: binsize_base
                !! The edges of the logarithmically distributed bins will be
                !! integer powers of this number. Default: 1.5.
            real(kind=dp), dimension(:), intent(inout), allocatable :: bin_mids
                !! Centres of bins.
            real(kind=dp), dimension(:), intent(inout), allocatable :: results
                !! Cluster number density in \eqref{eq:nsp}.

            integer :: num_bins, i, j, num_labels, sizeindex, spanning_label
            real(kind=dp) :: a, loga

            logical, dimension(:,:), allocatable :: binary_matrix
            integer, dimension(:,:), allocatable :: label_matrix
            integer, dimension(:), allocatable :: clustersizes, histogram
            real(kind=dp), dimension(:), allocatable :: bin_edges, bin_sizes

            if(present(binsize_base)) then
                a = binsize_base
            else
                a = 1.5d0
            end if

            !/cndstart/!
            loga = log(a)

            num_bins  = ceiling(log(1.0d0*L**2)/loga)
            bin_edges = a**[(i, i=0 ,num_bins)]
            bin_mids  = 0.5*(bin_edges(1:num_bins) + bin_edges(2:num_bins+1))
            bin_sizes = bin_edges(2:num_bins+1) - bin_edges(1:num_bins)

            allocate(histogram(1:num_bins))
            histogram = 0

            !$omp  parallel do &
            !$omp& private(binary_matrix, label_matrix, num_labels, &
            !$omp&         clustersizes, spanning_label, sizeindex) &
            !$omp& reduction(+:histogram)
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
            !$omp end parallel do

            results = histogram/(L**2 * num_samples * bin_sizes)
            !/cndend/!
        end subroutine

        function find_spanning_cluster(labelled_matrix, num_labels) result(spanning_label)
            !! Find the label of the percolating cluster, i.e. the one spanning from
            !! one side of the system to the opposing side.

            integer :: spanning_label
                !! Label of the percolating cluster.
            integer, dimension(:,:), intent(in) :: labelled_matrix
                !! Labelled matrix of clusters from [[hoshen_kopelman]]/[[label]].
            integer, intent(in) :: num_labels
                !! Known number of clusters.
            integer :: L

            L = size(labelled_matrix,1)

            spanning_label = find_intersection(labelled_matrix(:,1), labelled_matrix(:,L), num_labels)

            if (spanning_label == -1) then
                spanning_label = find_intersection(labelled_matrix(1,:), labelled_matrix(L,:), num_labels)
            end if
        end function

        function spanning_density_one_sample(p, L) result(spanning_density)
            real(kind=dp) :: spanning_density
            real(kind=dp), intent(in) :: p
            integer, intent(in) :: L

            logical, dimension(:,:), allocatable :: binary_matrix
            integer, dimension(:,:), allocatable :: labelled_matrix

            integer :: num_labels, spanning_label

            binary_matrix = create_binary_matrix(p,L)
            call label(binary_matrix, labelled_matrix, num_labels)

            spanning_label = find_spanning_cluster(labelled_matrix, num_labels)

            if (spanning_label == -1) then
                spanning_density = 0
                return
            end if

            spanning_density = count(labelled_matrix == spanning_label)/real(L**2,kind=dp)
        end function

        function spanning_density(p, L, num_samples)
            !! Density of the spanning/percolating cluster, i.e.
            !! the number of sites on the percolating cluster divided by \\(L^2\\).
            !! Averaged over **num_samples** Monte Carlo samples (with OpenMP).
            real(kind=dp) :: spanning_density
                !! The number of sites on the spanning cluster divided by \\(L^2\\).
            real(kind=dp), intent(in) :: p
                !! The probability for a site to allow transport.
            integer, intent(in) :: L
                !! The size of the system.
            integer, intent(in) :: num_samples
                !! The number of Monte Carlo samples.

            integer :: i
            real(kind=dp), dimension(:), allocatable :: results

            allocate(results(num_samples))

            !$omp parallel do
            do i=1,num_samples
                results(i) = spanning_density_one_sample(p, L)
            end do
            !$omp end parallel do

            spanning_density = sum(results)/num_samples
        end function

        function has_spanning_cluster_one_sample(p, L) result(has_spanning)
            logical :: has_spanning
            real(kind=dp), intent(in) :: p
            integer, intent(in) :: L

            logical, dimension(:,:), allocatable :: binary_matrix
            integer, dimension(:,:), allocatable :: labelled_matrix

            integer :: num_labels, spanning_label

            binary_matrix = create_binary_matrix(p,L)
            call label(binary_matrix, labelled_matrix, num_labels)

            spanning_label = find_spanning_cluster(labelled_matrix, num_labels)

            has_spanning = spanning_label /= -1
        end function

        function spanning_probability(p, L, num_samples)
            !! Calculate the probability of having a spanning/percolating cluster, given
            !! a system size **L** and probability for a site to have transport **p**.
            real(kind=dp) :: spanning_probability
                !! The probability of having a percolating cluster, calculated as the number
                !! of times a percolating cluster is found, divided by the number of attempts
                !! (**num_samples**).
            real(kind=dp), intent(in) :: p
                !! Probability for a site to allow transport.
            integer, intent(in) :: L
                !! Size of the system.
            integer, intent(in) :: num_samples
                !! Number of Monte Carlo samples.

            integer :: i
            logical, dimension(:), allocatable :: results

            allocate(results(num_samples))

            !$omp parallel do
            do i=1,num_samples
                results(i) = has_spanning_cluster_one_sample(p, L)
            end do
            !$omp end parallel do

            spanning_probability = count(results)/real(num_samples,kind=dp)
        end function

        function spanning_probability_inverse(x, L, num_samples, tolerance) result(p_x)
            !! Find the inverse of [[spanning_probability]] by use
            !! of the bisection method.
            real(kind=dp), intent(in) :: x
                !! The value of [[spanning_probability]] for which the inverse
                !! is calculated.
            real(kind=dp), intent(in) :: tolerance
                !! Tolerance of approximation. The return value is within
                !! **tolerance**/2 of the correct (but numerical) value.
            integer, intent(in) :: L
                !! Size of the system.
            integer, intent(in) :: num_samples
                !! Number of Monte Carlo samples to use when evaluating [[spanning_probability]].
            real(kind=dp) :: p_x
                !! The inverse of [[spanning_probability]]. If [[spanning_probability]] is
                !! denoted \\(\Pi(p,L)\\), this function returns \\(p\\) such
                !! that \\(\Pi(p,L)=x\\).


            real(kind=dp) :: lower, upper, lowerPI, upperPI, mid, midPI

            !/invPIstart/!
            lower = 0
            upper = 1
            lowerPI = 0
            upperPI = 1

            do while(upper - lower > tolerance)
                mid = (lower + upper)/2
                midPI = spanning_probability(mid, L, num_samples)
                if(midPI > x) then
                    upper = mid
                    upperPI = midPI
                else
                    lower = mid
                    lowerPI = midPI
                end if
            end do
            p_x = mid
            !/invPIend/!

        end function

        !/labelsubroutinestart/!
        subroutine label_stackoverflow(matrix, labelled_matrix, num_labels)
            logical, dimension(:,:), allocatable, intent(in) :: matrix
            integer, dimension(:,:), allocatable, intent(inout) :: labelled_matrix
            integer, intent(inout) :: num_labels
            integer :: L, i, j

            L = size(matrix,1)
            num_labels = 0

            if(.not. allocated(labelled_matrix)) then
                allocate(labelled_matrix(L,L))
            end if

            labelled_matrix = 0

            do j=1,L
                do i=1,L
                    if(matrix(i,j) .and. labelled_matrix(i,j) == 0) then
                        num_labels = num_labels + 1
                        call growcluster(matrix,labelled_matrix,i,j,num_labels)
                    end if
                end do
            end do
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
                end if
            end if
            if(j<L) then
                if(matrix(i,j+1) .and. labelled_matrix(i,j+1)==0) then
                    call growcluster(matrix, labelled_matrix, i, j+1, label)
                end if
            end if
            if(i>1) then
                if(matrix(i-1,j) .and. labelled_matrix(i-1,j)==0) then
                    call growcluster(matrix, labelled_matrix, i-1, j, label)
                end if
            end if
            if(j>1) then
                if(matrix(i,j-1) .and. labelled_matrix(i,j-1)==0) then
                    call growcluster(matrix, labelled_matrix, i, j-1, label)
                end if
            end if
        end subroutine
        !/growclustersubroutineend/!
end module percolation
