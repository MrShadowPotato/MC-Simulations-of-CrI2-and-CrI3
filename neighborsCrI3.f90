! This program simulates atomic structures by reading atom positions from a file and finding their neighbors.
program simulation
    implicit none

    integer :: n  ! Number of atoms
    
    ! Allocate arrays to store atom positions and neighbor information.
    real(8), allocatable, dimension(:,:) :: atoms
    integer, allocatable, dimension(:,:) :: neighbors
    character(len=2) :: element

    !Primitive vectors
    real(8), parameter :: a1(3) = [6.8911914825, 0.0, 0.0]
    real(8), parameter :: a2(3) = [-3.4488059462999998, 5.9711638719, 0.0]
    
    integer :: i
    real(8) :: x, y, z
    
    ! Open and read the atoms.xyz file containing the number of atoms and their positions.
    open(10, file='structures/SpinPositionsCrI3.xyz', status='old')
    read(10,*) n
    ! Read in atom positions from file
    read(10, '(A)') ! Skip the second line
    allocate(atoms(n,3))
    do i = 1, n
        read(10,*) element, x, y, z
        atoms(i,:) = [x, y, z]
    end do
    close(10)
    
    ! Allocate the neighbors array and find the neighbors for each atom.
    allocate(neighbors(n,3))
    neighbors = find_neighbors(atoms, 4.0d0, 3, a1*3, a2*3)

    ! Print the neighbors for each atom.
    open(1, file='structures/neighborsCrI3.txt', status='replace')
    do i = 1, n
        write(1,*) i, neighbors(i,1), neighbors(i,2), neighbors(i,3)
    end do
    close(1)
contains

    ! This function calculates the Euclidean distance between two points (v1 and v2).
    function distance(v1, v2) result(dist)
        implicit none
        real(8), dimension(:), intent(in) :: v1, v2
        real(8) :: dist
        !dist = sqrt(sum((v1-v2)**2))
        dist = sqrt(dot_product(v1-v2, v1-v2))
    end function distance

    ! This function finds the neighbors of atoms within a given distance (max_distance) and up to a maximum number (max_neighbors).
    function find_neighbors(vectors, max_distance, max_neighbors, vx, vy) result(neighbors)
        implicit none
        real(8), dimension(:,:), intent(in) :: vectors
        real(8), dimension(3), intent(in) :: vx, vy !Lattice vectors for periodic boundary conditions 
        integer, intent(in) :: max_neighbors
        real(8), intent(in) :: max_distance
        integer :: i, j, n, count, n1, n2
        integer, dimension(size(vectors,1), max_neighbors) :: neighbors
        real(8) :: dist ! Distance between two atoms
        n = size(vectors,1)
        
        ! Loop through all atoms and find their neighbors within the specified distance.
        do i = 1, n
            count = 0
            do j = 1, n
                do n1 = -1, 1
                    do n2 = -1, 1
                        ! Calculate the distance between the two atoms.
                        dist = distance(vectors(i,:), vectors(j,:) + n1*vx + n2*vy)
                        ! Check if the atoms are not the same and if their distance is less than the maximum allowed distance.
                        if (i /= j .and. dist <  max_distance) then
                            count = count + 1
                            ! If the atom already has the maximum number of neighbors, exit the loop.
                            if (count > max_neighbors) then
                                exit
                            end if
                            ! Store the index of the neighboring atom.
                            neighbors(i,count) = j
                        end if
                    end do
                end do
            end do
        end do
    end function find_neighbors

end program simulation