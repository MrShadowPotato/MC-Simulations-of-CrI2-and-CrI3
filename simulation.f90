program main
    use variables 
    implicit none
    


    ! Define variables.
    integer :: natoms, i
    real :: start_time, end_time, elapsed_time

    real(8), dimension(:,:), allocatable :: coordinates, Cr_coordinates, neighbors, spins
  
    call read_parameters()
    ! Calculate the number of atoms in the structure.
    natoms = nx*ny*8
    call cpu_time(start_time)

    coordinates = generate_structure(basis_CrI3, primitiveCrI3, nx, ny)
    !print *, coordinates
    call write_structures('CrI3.xyz', 'CrI3withoutI3.xyz', natoms, nx, ny, 'CrI3 structure', elementsCrI3, coordinates)

    ! Now that we have created the files containing the coordinates,
    ! we can use them to create the array that will hold the neighbors for each Cr atom.
    Cr_coordinates =  xyz_to_array('CrI3withoutI3.xyz')
    neighbors = find_neighbors(Cr_coordinates, 4.0d0, 3, vlatticeCrI3(1,:), vlatticeCrI3(2,:))
    !open(1, file='CrI3neighbors.txt', status='unknown')
    !do i = 1, nx*ny*2 !cambiar a variable
    !    write(1,'(3(I7))') int(neighbors(i,1)), int(neighbors(i,2)), int(neighbors(i,3))
    !end do
    !close(1)
    write(*,*) 'Now we generate the spins'
    spins = generate_spins(nx*ny*2)

    !(mcs, cratoms, nx, ny, T, J, H, spins, neighbors) 
    
    call simulation(mcs, nx*ny*2, nx, ny, temperature, exchange, 0.0d0, spins, int(neighbors))
    print *, 'This is a mesagge to check that the program has finished.'

    call cpu_time(end_time)
    elapsed_time = end_time - start_time

    open(unit=14, file=output_file, status="old", position='append', action="write")
    write(14, *) 'This is the time it took to run the program: ', elapsed_time
    close(14)
    write(6, *) 'This is the time it took to run the program: ', elapsed_time
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
    open(15, file='CrI3neighbors.txt', status='unknown')
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
        ! Print the neighbors of the atom to a file.
        write(15,'(3(I7))') int(neighbors(i,1)), int(neighbors(i,2)), int(neighbors(i,3))
        write(6, *) 'Neighbors of atom ', i, ' are ', int(neighbors(i,1)), int(neighbors(i,2)), int(neighbors(i,3))
    end do
    close(15)
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
    
    !write(*,*) "Choose an option: (1) ramdom directions or (2) fixed direction"
    !read(*,*) choice
    choice = 1
    
    if (choice == 1) then
        !Each spin is a different randomly generated normal vector.
        do i = 1, natoms
            spins(i,:) = random_normal_vector()
        end do
        
    else if (choice ==2) then
        !All the spins are the same randomly generated normal vector.
        spin = random_normal_vector()
        spins(:,1) = spin(1)
        spins(:,2) = spin(2)
        spins(:,3) = spin(3)
    end if
end function generate_spins


function calculate_initial_energy(spins, neighbors, J) result(system_energy)
    implicit none
    real(8), intent(in) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    real(8), intent(in) :: J
    real(8) :: system_energy
    integer :: i, k
    system_energy = 0.0d0
    do i = 1, size(spins,1)
        do k = 1, size(neighbors,2)
            system_energy = system_energy - J * dot_product(spins(i,:), spins(neighbors(i,k),:))
        end do
    end do
    system_energy = system_energy/2 ! Each interaction is counted twice
end function calculate_initial_energy


function calculate_initial_magnetization_vector(spins) result(magnetization_vector)
    implicit none
    real(8), intent(in) :: spins(:,:)
    real(8), dimension(3) :: magnetization_vector
    !real(8) :: magnetization
    integer :: i
    magnetization_vector = 0.0d0
    do i =1, size(spins,1)
        magnetization_vector = magnetization_vector + spins(i,:)
    end do 
    !magnetization_vector = magnetization_vector / size(spins,1)
    end function calculate_initial_magnetization_vector


function calculate_spin_energy(spins, neighbors, index, J) result(energy)
    real(8), intent(in) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    integer, intent(in) :: index
    real(8), intent(in) :: J
    integer :: i
    real(8) :: energy 
    energy = 0.0d0 
    do i = 1, size(neighbors,2)
        energy = energy - 2 * J * dot_product(spins(index,:), spins(neighbors(index,i),:))
    end do
end function calculate_spin_energy


function flip_ith_spin(spins, ith) result(new_spins)
    implicit none
    real(8), dimension(:,:), intent(inout) :: spins
    real(8), dimension(size(spins,1),3) :: new_spins
    integer, intent(in) :: ith
    new_spins = spins
    new_spins(ith,:) = random_normal_vector()
end function flip_ith_spin


function energy_change(old_system, new_system, neighbors, J, index) result(delta_E)
    !Calculates the energy change due to spin flip.
    !Compares the energy of the old system with the energy of the new system.
    !Only the energy of the spin that has been flipped is calculated.
    implicit none 
    real(8), intent(in), dimension(:,:) :: old_system, new_system
    integer, intent(in), dimension(:,:) :: neighbors
    real(8), intent(in) :: J
    integer, intent(in) :: index
    integer :: i
    real(8) :: original_energy, new_energy, delta_E
    real(8) :: new_spin(3), original_spin(3)

    original_energy = 0.0d0  
    new_energy = 0.0d0
    
    original_spin = old_system(index,:)
    new_spin = new_system(index,:)
    
    do i = 1, size(neighbors,2)
        original_energy = original_energy - J * dot_product(original_spin, spins(neighbors(index,i),:))
        new_energy = new_energy - J * dot_product(new_spin, spins(neighbors(index,i),:))
    end do

    delta_E = (new_energy - original_energy)
end function energy_change

 

function accept_change(delta_E, t) result(accept)
    implicit none
    real(8), intent(in) :: delta_E, t    
    logical :: accept
    
    if (delta_E < 0.0d0) then
        accept = .true.
    else
        accept = exp(-delta_E / (kB * t)) > ran2(seed)
    end if
end function accept_change




! This subroutine simulates a Monte Carlo simulation of a system of atoms with spins.
subroutine simulation(mcs, CrAtoms, nx, ny, T, J, H, spins, neighbors) 
    implicit none
    
    !Declare input variables
    integer, intent(in) :: mcs, CrAtoms, nx, ny
    real(8), intent(in) :: T, J, H
    real(8), intent(inout) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    
    !Declare local variables
    integer :: i, k, spins_index
    real(8), dimension(CrAtoms, 3) :: old_system, new_system
    real(8) :: current_energy, current_magnetization, delta_E
    real(8), dimension(3) :: current_magnetization_vector, final_magnetization_vector
    
    
    open(12, file=output_file, status='unknown')
    write(12,*) '#seed - mcs - temperature - CrAtoms - nx - ny - H'
    write(12,*) '#', seed, mcs, T, CrAtoms, nx, ny, H
    write(12,*) '#Mcs', 'Energy', 'Magnetization'
    !close(12)
    !stop
    !Loop over the number of Monte Carlo steps
    current_energy = calculate_initial_energy(spins, neighbors, J)
    current_magnetization_vector = calculate_initial_magnetization_vector(spins)
    current_magnetization = sqrt(dot_product(current_magnetization_vector, current_magnetization_vector))

    do i = 1, mcs
        print *, 'Mcs: ', i
        !Loop over the number of atoms
        do k = 1, CrAtoms
            !Generate the index of a random spin to be flipped
            spins_index = random_integer(1, CrAtoms)
            
            
            !Store the old system and generate a new one with one spin flipped
            old_system = spins
            new_system = flip_ith_spin(spins, spins_index)
            
            !Calculate the energy and magnetization of the old and new system
            
            delta_E = energy_change(old_system, new_system, neighbors, J, spins_index)
            
            !Decide wheter to accept the new system or not
            if (accept_change(delta_E, T)) then
                !If the new system is accepted, update the spins and the total energy
                spins = new_system
                current_energy = current_energy + delta_E
                current_magnetization_vector = current_magnetization_vector &
                - old_system(spins_index,:) + new_system(spins_index,:)
            end if
            current_magnetization = sqrt(dot_product(current_magnetization_vector, current_magnetization_vector))
            current_magnetization = current_magnetization / CrAtoms
            if (k == 1) then 
                write(12, '(I7 ,2(1x, F16.9))') i, current_energy, current_magnetization
            end if
        end do


    end do
    close(12)
    
    print *, 'End of the simulation'
end subroutine simulation




function random_normal_vector() result(vector)
    implicit none
    real(8), dimension(3) :: vector
    vector(1) = 2*ran2(seed)-1
    vector(2) = 2*ran2(seed)-1
    vector(3) = 2*ran2(seed)-1
    vector = vector / sqrt(dot_product(vector, vector))
end function random_normal_vector


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