program Cr_CrIx_structures
    use read_input
    implicit none
    real(8), dimension(:,:), allocatable :: Cr_coordinates!, neighbors
    
    call read_parameters()
    
    Cr_coordinates = write_CrIx(basis, primitive, nx, ny)
    neighbors = find_neighbors(Cr_coordinates, neighbor_min_distance, neighbor_max_distance, &
    max_neighbors, vlattice(1,:), vlattice(2,:))

    
contains

function write_CrIx(basis, primitive, nx, ny) result(coordinates)
    implicit none
    real(8), dimension(:,:), intent(in) :: basis
    real(8), dimension(2,3), intent(in) :: primitive
    integer, intent(in) :: nx, ny
    real(8), dimension(3) :: coordinate
    real(8), dimension(nx*ny*Cr_per_cell,3) :: coordinates
    integer :: i, j, k, count
    character(10) :: nx_num
    write(nx_num, '(I0)') nx
    open(1, file='../Cr_coordinates/'//'n'//trim(nx_num)//'Cr_'//compound//'.xyz', status='unknown')
    open(2, file='../CrIx_coordinates/'//'n'//trim(nx_num)//compound//'.xyz', status='unknown')
    write(1,*) nx*ny*Cr_per_cell
    write(1,*) 'Cr positions for CrIx with nx and ny respectively', nx, ny
    write(2,*) nx*ny*size(elements)
    write(2,*) 'Cr positions for CrIx with nx and ny respectively', nx, ny
    count = 0
    do i = 1, nx
        do j = 1, ny
            do k = 1, size(elements)
                coordinate = basis(k,:) + primitive(1,:)*REAL(i) + primitive(2,:)*REAL(j)
                if (elements(k) == 'Cr') then
                    count = count + 1
                    coordinates(count,:) = coordinate(:)
                    write(1, '(A2, 3(F16.6))') elements(k), coordinate(1), coordinate(2), coordinate(3)
                end if

                write(2, '(A2, 3(F16.6))') elements(k), coordinate(1), coordinate(2), coordinate(3)
            end do
        end do
    end do
    close(1)
    close(2)
end function write_CrIx


function find_neighbors(vectors, min_distance, max_distance, max_neighbors, vx, vy) result(neighbors)
    implicit none
    real(8), dimension(:,:), intent(in) :: vectors
    real(8), dimension(3), intent(in) :: vx, vy !Lattice vectors for periodic boundary conditions 
    real(8), dimension(3) :: v1, v2
    integer, intent(in) :: max_neighbors
    real(8), intent(in) :: min_distance, max_distance
    integer :: i, j, n, count, n1, n2
    integer, dimension(size(vectors,1), max_neighbors) :: neighbors
    real(8) :: dist ! Distance between two atoms
    character(10) :: nx_num
    write(nx_num, '(I0)') nx
    n = size(vectors,1)
    open(15, file='../neighbors/n'//trim(nx_num)//'Cr_'//compound//'neighbors.txt', status='unknown')
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
                    if (i /= j .and. dist <  max_distance .and. dist > min_distance) then
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
        write(15,*) neighbors(i,:)
        ! Print the neighbors of the atom to a file.
        !write(15,'(3(I8))') int(neighbors(i,1)), int(neighbors(i,2)), int(neighbors(i,3))
        !write(6, *) 'Neighbors of atom ', i, ' are ', int(neighbors(i,1)), int(neighbors(i,2)), int(neighbors(i,3))
        
    end do
    close(15)
end function find_neighbors
end program Cr_CrIx_structures