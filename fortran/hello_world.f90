program hello
    implicit none
    real :: t
    integer :: i
    t = 98.6
    do i = 1,10
        print *, "Hello World!"
        print *, t
        t = t + 10.0
    end do
end program hello
