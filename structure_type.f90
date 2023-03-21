module structures
    implicit none
    type :: compound
        integer :: natoms, nx, ny
        real(8), dimension(2,3) :: primitive
        real(8), dimension(2,3) :: latticevectors
        real(8), allocatable, dimension(:,:) :: basis 
        character(len=2), allocatable, dimension(:) :: elements
    end type compound
end module structures

