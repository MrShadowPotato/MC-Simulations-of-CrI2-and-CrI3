module variables
    implicit none
    integer, save, public :: seed = -123456789
    real(8), save, public :: kB = 8.617333262D-10
    integer, parameter :: mcs = 100
    real(8) :: temperature = 5
    real(8), parameter :: exchange = 2.7D-3
    integer, parameter :: nx = 3, ny = 3

    
end module variables