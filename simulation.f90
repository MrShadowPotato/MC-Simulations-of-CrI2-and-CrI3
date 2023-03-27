program main
    use variables 
    implicit none
    


    ! Define variables.
    integer :: natoms, i
    real(8), dimension(:,:), allocatable :: coordinates, Cr_coordinates, neighbors, spins
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


    !write(*,*) 'This program generates the structures of CrI3.'
    !write(*,*) 'Please provide the number of unit cells for the x direction:'
    !read(5,*) nx
    !write(*,*) 'Please provide the number of unit cells for the y direction:'
    !read(5,*) ny




    ! Define the lattice vectors for periodic boundary conditions.
    vlatticeCrI3(1,:) = [primitiveCrI3(1,1)*nx, primitiveCrI3(1,2)*nx, primitiveCrI3(1,3)*nx]
    vlatticeCrI3(2,:) = [primitiveCrI3(2,1)*ny, primitiveCrI3(2,2)*ny, primitiveCrI3(2,3)*ny]
    vlatticeCrI2(1,:) = [primitiveCrI2(1,1)*nx, primitiveCrI2(1,2)*nx, primitiveCrI2(1,3)*nx]
    vlatticeCrI2(2,:) = [primitiveCrI2(2,1)*ny, primitiveCrI2(2,2)*ny, primitiveCrI2(2,3)*ny]

    ! Calculate the number of atoms in the structure.
    natoms = nx*ny*8

    coordinates = generate_structure(basis_CrI3, primitiveCrI3, nx, ny)
    !print *, coordinates
    call write_structures('CrI3.xyz', 'CrI3withoutI3.xyz', natoms, nx, ny, 'CrI3 structure', elementsCrI3, coordinates)

    ! Now that we have created the files containing the coordinates,
    ! we can use them to create the array that will hold the neighbors for each Cr atom.
    Cr_coordinates =  xyz_to_array('CrI3withoutI3.xyz')
    neighbors = find_neighbors(Cr_coordinates, 4.0d0, 3, vlatticeCrI3(1,:), vlatticeCrI3(2,:))
    open(1, file='CrI3neighbors.txt', status='unknown')
    do i = 1, nx*ny*2 !cambiar a variable
        write(1,'(3(I7))') int(neighbors(i,1)), int(neighbors(i,2)), int(neighbors(i,3))
    end do
    close(1)
    write(*,*) 'Now we generate the spins'
    spins = generate_spins(nx*ny*2)
    !(mcs, natoms, T, J, H, spins, neighbors) 
    
    call simulation(mcs, nx*ny*2, temperature, exchange, 0.0d0, spins, int(neighbors))
    print *, 'This is a mesagge to check that the program has finished.'


contains 

! This subroutine writes the coordinates to a xyz file.
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
            if ( elements(j) =='Cr' ) then
                write(2, '(A2, 3(F16.10))') elements(j), coordinates(count,1), coordinates(count,2), coordinates(count,3)
            end if
            write(1, '(A2, 3(F16.10))') elements(j), coordinates(count,1), coordinates(count,2), coordinates(count,3)
            !write(6, *) count, elements(j), coordinates(count,1), coordinates(count,2), coordinates(count,3)
        end do
    end do
    close(1)
    close(2)
end subroutine write_structures

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

subroutine write_system_var(energies,  magnetizations)
    implicit none
    real(8), dimension(:), intent(in) :: energies, magnetizations
    integer :: i
    open(1, file='data.txt', status='unknown')
    do i = 1, size(energies)
        write(1, '(2(F16.9))') energies(i), magnetizations(i)
    end do
    close(1)
    end subroutine write_system_var

function calculate_magnetization(spins) result(magnetization)
    !Check logic of this function.
    implicit none
    real(8), dimension(:,:), intent(in) :: spins
    real(8), dimension(3) :: magnetization_vector
    real(8) :: magnetization
    integer :: i
    magnetization_vector = 0.0d0
    do i = 1, size(spins,1)
        magnetization_vector = magnetization_vector + spins(i,:)
    end do
    magnetization_vector =  magnetization_vector / size(spins,1)
    magnetization = sqrt(dot_product(magnetization_vector, magnetization_vector))
    end function calculate_magnetization

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
                        if (count > max_neighbors) then !Probably not necessary
                            exit
                        end if
                        ! Store the index of the neighboring atom.
                        neighbors(i, count) = int(j)
                    end if
                end do
            end do
        end do
    end do
end function find_neighbors

! This function generates the coordinates of the atoms in the structure.
function generate_structure(basis, primitive, nx, ny) result(coordinates)
    implicit none
    real(8), dimension(:,:), intent(in) :: basis
    real(8), dimension(2,3), intent(in) :: primitive
    integer, intent(in) :: nx, ny
    real(8), dimension(nx*ny*size(basis,1),3) :: coordinates
    integer :: i, j, k, n, count
    n = size(basis,1)
    count = 0
    do i = 1, nx
        do j = 1, ny
            do k = 1, size(basis,1)
                count = count + 1
                coordinates(count,:) = basis(k,:) + primitive(1,:)*REAL(i) + primitive(2,:)*REAL(j)
            end do
        end do
    end do
end function generate_structure


function generate_spins(natoms) result(spins)
    integer, intent(in) :: natoms
    real(8), dimension(natoms,3):: spins
    integer :: choice, i
    real(8), dimension(3) :: spin

    write(*,*) "Choose an option: (1) ramdom directions or (2) fixed direction"
    read(*,*) choice
    if (choice == 1) then
        do i = 1, natoms
            !call random_number(spin)
            spin(1) = ran2(seed)
            spin(2) = ran2(seed)
            spin(3) = ran2(seed)
            spins(i,:) = spin / sqrt(dot_product(spin, spin))
        end do
    else if (choice ==2) then
        !call random_number(spin)
        spin(1) = ran2(seed)
        spin(2) = ran2(seed)
        spin(3) = ran2(seed)
        spin = spin / sqrt(dot_product(spin, spin))
        spins(:,1) = spin(1)
        spins(:,2) = spin(2)
        spins(:,3) = spin(3)
    end if
end function generate_spins

function calculate_energy(spins, neighbors, J) result(energy)
    implicit none
    real(8), intent(in) :: spins(:,:) 
    integer, intent(in), dimension(:,:) :: neighbors
    real(8), intent(in):: J
    real(8) :: energy
    integer :: i, k
    energy = 0.0d0
    do i = 1, size(spins,1)
        do k = 1, size(neighbors,2)
            energy = energy - J * dot_product(spins(i,:), spins(neighbors(i,k),:)) !Check whether the sign is correct!!!!!!!!!!!!!!!!
        end do
    end do
    energy = energy/2 ! Each interaction is counted twice
    !print *, energy
end function calculate_energy

function accept_change(old_E, new_E, t) result(accept)
    implicit none
    real(8), intent(in) :: old_E, new_E, t
    real(8) :: delta_E
    logical :: accept
    !real(8), parameter :: kB = 8.617333262E-10 

    delta_E = new_E - old_E
    if (delta_E < 0.0d0) then
        accept = .true.
    else
        accept = exp(-delta_E / (kB * t)) > ran2(seed)
    end if
end function accept_change


! This function flips a random spin in the system.
function flip_spin(spins) result(new_spins)
    real(8), dimension(:,:), intent(inout) :: spins                                                        
    real(8), dimension(size(spins,1),3) :: new_spins
    integer :: index
    real(8), dimension(3) :: spin
    new_spins = spins
    index = random_integer(1, size(spins,1))
    !call random_number(spin)
    spin(1) = ran2(seed)
    spin(2) = ran2(seed)
    spin(3) = ran2(seed)    
    new_spins(index,:) = spin / sqrt(dot_product(spin, spin))
end function flip_spin


! This subroutine simulates a Monte Carlo simulation of a system of atoms with spins.
subroutine simulation(mcs, natoms, T, J, H, spins, neighbors) 
    implicit none

    !Declare input variables
    integer, intent(in) :: mcs, natoms
    real(8), intent(in) :: T, J, H
    real(8), intent(inout) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors

    !Declare local variables
    integer :: i, k
    real(8), dimension(natoms, 3) :: old_system, new_system
    real(8) :: old_E, new_E, total_E, old_magnetization
    real(8), dimension(mcs) :: energies, magnetizations


    open(12, file='data.txt', status='unknown')
    write(12,*) 'Mcs', 'Energy', 'Magnetization'
    write(12,*) '-----------------------------'
    !close(12)
    !stop
    !Loop over the number of Monte Carlo steps
    do i = 1, mcs
        print *, 'Mcs: ', i
        !Loop over the number of atoms
        do k = 1, natoms

            !Store the old system and generate a new one with one spin flipped
            old_system = spins
            new_system = flip_spin(spins)
            
            !Calculate the energy of the old and new system
            old_E = calculate_energy(old_system, neighbors, J)
            new_E = calculate_energy(new_system, neighbors, J)
            old_magnetization = calculate_magnetization(spins)
            
            !Decide wheter to accept the new system or not
            if (accept_change(old_E, new_E, T)) then
                !If the new system is accepted, update the spins and the total energy
                spins = new_system
            end if
            !Record the total energy
            !total_E = total_E + old_E !Maybe it should be new_E
            if (k == natoms - 1) then
                print *, 'k = ', k
            end if
        end do
        !energies(i) = old_E
        !magnetizations(i) = old_magnetization
        write(12, '(I7 ,2(1x, F16.9))') i, old_E, old_magnetization
        !write(6, *) i, old_E, old_magnetization

        !stop
        if ((i == mcs - 1).and.(k == natoms - 1)) then 
            print *, 'i =', mcs - 1
        end if
    end do
    close(12)
    !call write_system_var(energies, magnetizations)

    print *, 'Avg energy: ', total_E/mcs
end subroutine simulation

function random_integer(a, b) result(rand_int)
    integer, intent(in) :: a, b
    integer :: rand_int
    real :: rand_real
    ! Generate a random real number between 0 and 1 using RANDOM_NUMBER
    !call RANDOM_NUMBER(rand_real)
    rand_real = ran2(seed)
    ! Scale the random number to be between a and b
    rand_int = a + int(real(b - a + 1) * rand_real)
end function random_integer

DOUBLE PRECISION FUNCTION ran2(idum)
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  INTEGER(K4B), INTENT(INOUT) :: idum
  !"Minimal" random number generator of Park and Miller combined with a 
  !Marsaglia shift sequence. Returns a uniform random deviate between 0.0 and 
  !1.0 (exclusive of the endpoint values). This fully portable, scalar 
  !generator has the "traditional" (not Fortran 90) calling sequence with a 
  !random deviate as the returned function value: call with idum a negative 
  !integer to initialize; thereafter, do not alter idum except to reinitialize.
  !The period of this generator is about 3.1� 10^18.
  INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  REAL, SAVE :: am
  INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
  if (idum <= 0 .or. iy < 0) then    !Initialize.
     am=nearest(1.0,-1.0)/IM
     iy=ior(ieor(888889999,abs(idum)),1)
     ix=ieor(777755555,abs(idum))
     idum=abs(idum)+1    !Set idum positive.
  end if
  ix=ieor(ix,ishft(ix,13))   !Marsaglia shift sequence with period 232 - 1.
  ix=ieor(ix,ishft(ix,-17))
  ix=ieor(ix,ishft(ix,5))
  k=iy/IQ                         !Park-Miller sequence by Schrage's method,
  iy=IA*(iy-k*IQ)-IR*k            !period 231 - 2.
  if (iy < 0) iy=iy+IM
  ran2=am*ior(iand(IM,ieor(ix,iy)),1)  !Combine the two generators with masking 
END FUNCTION ran2                      !to ensure nonzero value.




end program main