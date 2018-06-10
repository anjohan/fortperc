module randomwalk
    use utilities
    use percolation
    implicit none
    private :: find_random_point
    contains

        function one_random_walker(matrix, num_steps) result(displacement)
            !! Let one random walker do **num_steps** jumps on the `.true.`
            !! values of **matrix**. The initial position is randomly selected.

            logical, dimension(:,:), intent(in) :: matrix
                !! Matrix whose `.true.` values are allowed positions for the
                !! random walker.

            integer :: num_steps
                !! Number of steps for the random walker to take.

            real(kind=dp), dimension(:,:), allocatable :: displacement
                !! Array of displacements, which has dimension
                !! 2 \\(\times\\) (**num_steps** + 1),
                !! and should logically have `dimension(2,0:num_steps)`, such that
                !! the i'th column contains the displacement after i steps.
                !! Averaged over all systems and all walkers.

            integer :: x0, y0, x, y, x_new, y_new, dx, dy, L, i, direction
            real(kind=dp) :: real_direction

            ! North - East - South - West
            integer, dimension(4) :: dxs = [0, 1, 0, -1]
            integer, dimension(4) :: dys = [1, 0, -1, 0]

            L = size(matrix, 1)
            allocate(displacement(2,0:num_steps))

            displacement(:,:) = 0

            call find_random_point(matrix, x0, y0)

            x = x0
            y = y0

            do i = 1, num_steps
                call random_number(real_direction)
                direction = int(4*real_direction + 1)

                dx = dxs(direction)
                dy = dys(direction)

                x_new = x + dx
                y_new = y + dy

                call periodic_wraparound(x_new, L)
                call periodic_wraparound(y_new, L)

                if(matrix(x_new, y_new)) then
                    x = x_new
                    y = y_new
                    displacement(:,i) = displacement(:,i-1) + [dx, dy]
                else
                    displacement(:,i) = displacement(:,i-1)
                end if

            end do
        end function

        function random_walkers(p, L, num_systems, num_walkers, num_steps) result(displacement)
            !! Start **num_walkers** on the percolating cluster of each of **num_systems**
            !! systems, and return the averaged displacement.

            real(kind=dp), intent(in) :: p
                !! The probability for each site to allow transport.
            integer, intent(in) ::  L
                !! The system size.
            integer, intent(in) ::  num_systems
                !! The number of systems over which to average.
            integer, intent(in) ::  num_walkers
                !! The number of random walkers for each system over which to average.
            integer, intent(in) ::  num_steps
                !! The number of steps which the random walkers take.

            real(kind=dp), dimension(:), allocatable :: displacement
                !! Array of squared displacements, which has dimension
                !! (**num_steps** + 1),
                !! and should logically have `dimension(0:num_steps)`, such that
                !! the i'th element contains the displacement after i steps.
                !! Averaged over all systems and all walkers.


            logical, dimension(:,:), allocatable :: binary_matrix, spanning_cluster
            integer, dimension(:,:), allocatable :: label_matrix
            real(kind=dp), dimension(:,:), allocatable :: tmp_displacement

            integer :: i, num_spanning_clusters, num_clusters, spanning_cluster_label

            allocate(displacement(0:num_steps))
            displacement(:) = 0

            num_spanning_clusters = 0
            do while(num_spanning_clusters < num_systems)
                binary_matrix = create_binary_matrix(p, L)
                call label(binary_matrix, label_matrix, num_clusters)
                spanning_cluster_label = find_spanning_cluster(label_matrix, num_clusters)

                if(spanning_cluster_label < 1) then
                    cycle
                end if
                num_spanning_clusters = num_spanning_clusters + 1

                spanning_cluster = mark_percolating_with_periodic( &
                                        label_matrix, spanning_cluster_label, &
                                        num_clusters)


                !$omp parallel do private(tmp_displacement) reduction(+:displacement)
                do i = 1, num_walkers
                    tmp_displacement = one_random_walker(spanning_cluster, num_steps)
                    displacement(:) = displacement(:) + sum(tmp_displacement(:,:)**2,dim=1)
                end do
                !$omp end parallel do
            end do
            displacement(:) = displacement(:)/(num_walkers*num_systems)

        end function
end module randomwalk
