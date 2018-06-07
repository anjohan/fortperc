module randomwalk
    use utilities
    implicit none
    contains
        function one_random_walker(matrix, num_steps) result(displacement)
            !! Let one random walker do **num_steps** jumps on the `.true.`
            !! values of **matrix**. The initial position is randomly selected.

            logical, dimension(:,:), intent(in) :: matrix
                !! Matrix whose `.true.` values are allowed positions for the
                !! random walker.

            integer :: num_steps
                !! Number of steps for the random walker to take.

            real(kind=dp), dimension(:), allocatable :: displacement
            real(kind=dp) :: x0_real, y0_real
            integer :: x0, y0, x, y, dx, dy, L
            integer, dimension(:,:), allocatable :: tmp_displacement

            L = size(matrix, 1)
            allocate(tmp_displacement(2,0:L))

            tmp_displacement(:,0) = 0

            call random_seed()

            do
                call random_number(x0_real)
                call random_number(y0_real)
                x0_real = L*x0_real + 0.5
