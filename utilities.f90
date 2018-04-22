module utilities
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    contains
        function stringfromint(x)
            character(len=:), allocatable :: stringfromint
            integer, intent(in) :: x
            integer :: numdigits

            numdigits = int(log10(real(x))) + 1
            allocate(character(len=numdigits) :: stringfromint)
            write(unit=stringfromint, fmt="(i0)") x
        end function

        function linspace(a,b,N)
            integer, intent(in) :: N
            real(kind=dp), intent(in) :: a,b
            real(kind=dp), dimension(:), allocatable :: linspace
            real(kind=dp) :: dx
            integer :: i

            dx = (b - a)/(N - 1)
            linspace = [ (a + i*dx, i=0, N-1) ]
        end function
end module
