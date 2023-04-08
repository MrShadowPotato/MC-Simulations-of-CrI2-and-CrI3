module variables
    implicit none
    character(len=100) :: dummy
    integer, public :: seed, mcs, nx, ny, spins_orientation
    real(8), public :: temperature 
    real(8), save, public :: kB = 8.617333262D-2
    real(8), parameter :: exchange = 2.76

contains

    subroutine read_parameters()
        ! Open the file
        character(len=100) :: filename
        integer :: status
        filename = "params.txt"
        open(unit=10, file=filename, status="old", action="read", iostat=status)
        
        ! Read the parameters
        read(10,*) dummy, dummy, seed
        read(10,*) dummy, dummy, mcs
        read(10,*) dummy, dummy, temperature
        read(10,*) dummy, dummy, nx
        read(10,*) dummy, dummy, ny
        read(10,*) dummy, dummy, spins_orientation
        
        ! Close the file
        close(unit=10)
        
        ! Print the parameters
        write(*,*) "seed = ", seed
        write(*,*) "mcs = ", mcs
        write(*,*) "temperature = ", temperature
        write(*,*) "nx = ", nx
        write(*,*) "ny = ", ny
        write(*,*) "spins_orientation = ", spins_orientation
    end subroutine read_parameters
    
end module variables
