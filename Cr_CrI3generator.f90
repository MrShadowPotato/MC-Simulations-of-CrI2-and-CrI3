program structure
    use variables
    implicit none
    real(8), dimension(:,:), allocatable :: Cr_coordinates, neighbors
    integer :: i
    call read_parameters()
    Cr_coordinates = write_Cr(basis_CrI3, primitiveCrI3, nx, ny, Cr_xyz)
    
    neighbors = find_neighbors(Cr_coordinates, 4.0d0, 3, vlatticeCrI3(1,:), vlatticeCrI3(2,:), Cr_neighbors)

    


contains
function write_Cr(basis, primitive, nx, ny, Cr_file) result(coordinates)
    implicit none
    real(8), dimension(:,:), intent(in) :: basis
    real(8), dimension(2,3), intent(in) :: primitive
    character(len=*), intent(in) :: Cr_file
    integer, intent(in) :: nx, ny
    real(8), dimension(3) :: coordinate
    integer, parameter :: Cr_per_cell = 2
    real(8), dimension(nx*ny*Cr_per_cell,3) :: coordinates
    integer :: i, j, k, count
    open(2, file='Cr_coordinates/'//Cr_file, status='unknown')
    write(2,*) nx*ny*Cr_per_cell
    write(2,*) 'Cr positions for CrI3  with nx and ny respectively', nx, ny
    count = 0
    do i = 1, nx
        do j = 1, ny
            do k = 1, Cr_per_cell
                count = count + 1
                coordinate = basis(k,:) + primitive(1,:)*REAL(i) + primitive(2,:)*REAL(j)
                coordinates(count,:) = coordinate(:)
                write(6, '(A,I8, 3(F16.6))')'Cr#', count, coordinate(1), coordinate(2), coordinate(3)
                write(2, '(A2, 3(F16.6))') 'Cr', coordinate(1), coordinate(2), coordinate(3)
            end do
        end do
    end do
    close(2)
end function write_Cr


function xyz_to_array(file) result(coordinates)
    implicit none
    character(len=*), intent(in) :: file
    integer :: natoms
    real(8), dimension(:,:), allocatable :: coordinates
    integer :: i
    character(len=2) :: element
    open(1, file=file, status='old')
    read(1,*) natoms
    read(1,*) ! Skip the comment line.
    allocate(coordinates(natoms,3))
    do i=1, natoms
        read(1,*) element, coordinates(i,1), coordinates(i,2), coordinates(i,3)
    end do
    close(1)
end function xyz_to_array

function read_neighbors(file, Cr_atoms)  result(neighbors)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(in) :: Cr_atoms
    integer, dimension(:,:), allocatable :: neighbors
    integer :: i
    open(1, file='neighbors/'//file, status='old')
    allocate(neighbors(Cr_atoms, 3))
    do i = 1, Cr_atoms
            read(1, '(3(I8))') neighbors(i,1), neighbors(i,2), neighbors(i,3)
    end do
end function read_neighbors

function find_neighbors(vectors, max_distance, max_neighbors, vx, vy, neighbors_file) result(neighbors)
    implicit none
    real(8), dimension(:,:), intent(in) :: vectors
    character(len=*), intent(in) :: neighbors_file
    real(8), dimension(3), intent(in) :: vx, vy !Lattice vectors for periodic boundary conditions 
    real(8), dimension(3) :: v1, v2
    integer, intent(in) :: max_neighbors
    real(8), intent(in) :: max_distance
    integer :: i, j, n, count, n1, n2
    integer, dimension(size(vectors,1), max_neighbors) :: neighbors
    real(8) :: dist ! Distance between two atoms
    n = size(vectors,1)
    open(15, file='neighbors/'//neighbors_file, status='unknown')
    ! Loop through all atoms and find their neighbors within the specified distance.
    do i = 1, n
        count = 0
        do j = 1, n
            do n1 = -1, 1
                do n2 = -1, 1
                    ! Calculate the distance between the two atoms.
                    v1 = vectors(i,:)
                    v2 = vectors(j,:) + n1*vx + n2*vy
                    dist = sqrt(dot_product(v1 - v2, v1 - v2))
                    ! Check if the atoms are not the same and if their distance is less than the maximum allowed distance.
                    if (i /= j .and. dist <  max_distance) then
                        count = count + 1
                        ! If the atom already has the maximum number of neighbors, exit the loop.
                        if (count > max_neighbors) then !Probably not necessary
                            exit
                        end if
                        ! Store the index of the neighboring atom.
                        neighbors(i, count) = int(j)
                    end if
                end do
            end do
        end do
        ! Print the neighbors of the atom to a file.
        write(15,'(3(I8))') int(neighbors(i,1)), int(neighbors(i,2)), int(neighbors(i,3))
        write(6, *) 'Neighbors of atom ', i, ' are ', int(neighbors(i,1)), int(neighbors(i,2)), int(neighbors(i,3))
    end do
    close(15)
end function find_neighbors


subroutine write_structures(file1, file2, natoms, nx, ny, comment, elements, coordinates)
    implicit none
    character(len=*), intent(in) :: file1, file2
    integer, intent(in) :: natoms
    integer, intent(in) :: nx, ny
    character(len=*), intent(in) :: comment !Optional comment to be written to the file.
    character(len=2), dimension(:), intent(in) :: elements
    real(8), dimension(natoms,3), intent(in) :: coordinates
    integer :: i, j, ncells, count
    ncells = natoms / size(elements)
    open(1, file=file1, status='unknown')
    open(2, file=file2, status='unknown')
    write(1,*) natoms
    write(1,*) 'Numbers of cells for x and y:  ', nx, ny
    write(2,*) nx*ny*2
    write(2,*) 'Cr positions for CrI3  with nx and ny respectively', nx, ny
    count = 0 
    do i=1, ncells
        do j=1, size(elements)
            count = count + 1
            write(6, *) 'Atoms = ', count
            if ( elements(j) =='Cr' ) then
                write(2, '(A2, 3(F16.6))') elements(j), coordinates(count,1), coordinates(count,2), coordinates(count,3)
            end if
            write(1, '(A2, 3(F16.6))') elements(j), coordinates(count,1), coordinates(count,2), coordinates(count,3)
            !write(6, *) count, elements(j), coordinates(count,1), coordinates(count,2), coordinates(count,3)
        end do
    end do
    close(1)
    close(2)
end subroutine write_structures

end program structure