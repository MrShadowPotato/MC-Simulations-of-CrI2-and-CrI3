program structuresCrI3
    implicit none


    ! Define variables.
    integer :: natoms, i, j, k, l, m, n, count, n1, n2, nx, ny
    real(8), dimension(2, 3) :: primitiveCrI3, primitiveCrI2, vlatticeCrI3, vlatticeCrI2
    real(8), dimension(8, 3) :: basis_CrI3 ! Basis vectors for the CrI3 structure.
    real(8), dimension(6, 3) :: basis_CrI2 ! Basis vectors for the CrI2 structure.
    character(len=2), dimension(8) :: elementsCrI3 = (/ 'Cr', 'Cr', 'I ', 'I ', 'I ', 'I ', 'I ', 'I ' /)
    character(len=2), dimension(6) :: elementsCrI2 = (/ 'Cr', 'Cr', 'I ', 'I ', 'I ', 'I '/)

    ! Define the primiteve vectors for the CrI3 structure.
    primitiveCrI3(1, :) = [6.8911914825, 0.0, 0.0]
    primitiveCrI3(2, :) = [-3.4488059462999998, 5.9711638719, 0.0]
    ! Define the basis vectors for the CrI3 structure.
    basis_CrI3(1,:) = [3.42245459776438, 1.94691491308500, 9.72065830568899]   !Cr
    basis_CrI3(2,:) = [-0.02339545180016, 3.93735766285690, 9.72080135688876]  !Cr
    basis_CrI3(3,:) = [-1.03538548904423, 5.92724513968807, 11.28527832915810] !I 
    basis_CrI3(4,:) = [2.20273590406674, 3.81437945298318, 11.29027271402490]  !I
    basis_CrI3(5,:) = [5.64639472731378, 2.07416105561567, 11.29166984494761]  !I
    basis_CrI3(6,:) = [0.97989374494890, 5.92913055695139, 8.15211676904270]   !I
    basis_CrI3(7,:) = [1.19604098486642, 2.06898570270974, 8.15546418188776]   !I
    basis_CrI3(8,:) = [4.64469384948989, 3.81029176486010, 8.15107346274473]   !I
    
    ! Define the primitive vectors for the CrI2 structure.
    primitiveCrI2(1, :) = [3.9350889407837069, 0.0, 0.0]
    primitiveCrI2(2, :) = [0.0, 7.793970515734439, -0.0155822353114624]
    ! Degine the basis vectors for the CrI2 structure.
    basis_CrI2(1,:) = [0.00000000000000, 1.18048232888026, 15.66708540526766]  !Cr
    basis_CrI2(2,:) = [1.96754466301406, 5.07735235387682, 15.66009018229279]  !Cr
    basis_CrI2(3,:) = [1.96754466301406, 2.14455628147048, 17.30559462841125]  !I 
    basis_CrI2(4,:) = [1.96754466301406, 0.21647035484051, 14.02945080633759]  !I
    basis_CrI2(5,:) = [0.00000000000000, 6.04092814648801, 17.29708073114964]  !I
    basis_CrI2(6,:) = [0.00000000000000, 4.11332427420235, 14.02184874359694]  !I


    write(*,*) 'This program generates the structures of CrI3.'
    write(*,*) 'Please provide the number of unit cells for the x direction:'
    read(5,*) nx
    write(*,*) 'Please provide the number of unit cells for the y direction:'
    read(5,*) ny




    ! Define the lattice vectors for periodic boundary conditions.
    vlatticeCrI3(1,:) = [primitiveCrI3(1,1)*nx, primitiveCrI3(1,2)*nx, primitiveCrI3(1,3)*nx]
    vlatticeCrI3(2,:) = [primitiveCrI3(2,1)*ny, primitiveCrI3(2,2)*ny, primitiveCrI3(2,3)*ny]
    vlatticeCrI2(1,:) = [primitiveCrI2(1,1)*nx, primitiveCrI2(1,2)*nx, primitiveCrI2(1,3)*nx]
    vlatticeCrI2(2,:) = [primitiveCrI2(2,1)*ny, primitiveCrI2(2,2)*ny, primitiveCrI2(2,3)*ny]

    ! Calculate the number of atoms in the structure.
    natoms = nx*ny*8


    

contains 

!Needs to be modified!!!!!!!!!!!
subroutine write_xyz( filename, natoms, comment, elements, position )
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: natoms
    character(len=*), intent(in) :: comment
    character(len=2), dimension(natoms), intent(in) :: elements
    real(8), dimension(natoms,3), intent(in) :: position
    integer :: i, j, ncells
    ncells = natoms / size(elements)
    open(1, file=filename, status='unknown')
    write(1,*) natoms
    write(1,*) comment
    do i=1, ncells
        do j=1, size(elements)
            write(1, '(A2, 3(F16.10))') elements(j), position(i,1), position(i,2), position(i,3)
        end do
    end do
    close(1)
end subroutine write_xyz



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

function generate_structure(basis, primitive, nx, ny) result(structure)
    implicit none
    real(8), dimension(:,:), intent(in) :: basis
    real(8), dimension(2,3), intent(in) :: primitive
    integer, intent(in) :: nx, ny
    real(8), dimension(nx*ny*size(basis,1),3) :: structure
    integer :: i, j, k, l, m, n, count
    n = size(basis,1)
    count = 0
    do i = 1, nx
        do j = 1, ny
            do k = 1, size(primitive,1)
                    structure(count,:) = basis(l,:) + primitive(1,:)*nx + primitive(2,:)*ny 
            end do
        end do
    end do
end function generate_structure



end program structuresCrI3