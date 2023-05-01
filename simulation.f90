program main
    use variables 
    use rng
    implicit none
    

    ! Define variables.
    integer :: i
    real :: start_time, end_time, elapsed_time
    real(8), dimension(:,:), allocatable :: spins
    integer, dimension(:,:), allocatable :: neighbors
  
    call read_parameters()
    call cpu_time(start_time)


    neighbors = read_neighbors(Cr_neighbors, Cr_atoms)
 
    write(*,*) 'Now we generate the spins'
    spins = generate_spins(Cr_atoms, spins_orientation)


    
    call simulation(mcs, Cr_atoms, nx, ny, temperature, exchange, spins, neighbors)
    print *, 'This is a mesagge to check that the program has finished.'

    call cpu_time(end_time)
    elapsed_time = end_time - start_time

    open(unit=14, file='data/stabilizer/'//output_file, status="old", position='append', action="write")
    write(14, *) 'This is the time it took to run the program: ', elapsed_time
    close(14)
    write(6, *) 'This is the time it took to run the program: ', elapsed_time
contains 


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




function generate_spins(natoms, choice) result(spins)
    integer, intent(in) :: natoms, choice
    real(8), dimension(natoms,3):: spins
    integer :: i
    real(8), dimension(3) :: spin
    
    if (choice == 1) then
        !Each spin is a different randomly generated normal vector.
        do i = 1, natoms
            spins(i,:) = random_normal_vector(seed)
        end do
        
    else 
        !All the spins are the same randomly generated normal vector.
        spin = initial_magnetization_vector(:)
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


function flip_ith_spin(spins, ith) result(new_spins)
    implicit none
    real(8), dimension(:,:), intent(inout) :: spins
    real(8), dimension(size(spins,1),3) :: new_spins
    integer, intent(in) :: ith
    new_spins = spins
    new_spins(ith,:) = random_normal_vector(seed)
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
subroutine simulation(mcs, CrAtoms, nx, ny, T, J, spins, neighbors) 
    implicit none
    
    !Declare input variables
    integer, intent(in) :: mcs, CrAtoms, nx, ny
    real(8), intent(in) :: T, J
    real(8), intent(inout) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    
    !Declare local variables
    integer :: i, k, spins_index
    real(8), dimension(CrAtoms, 3) :: old_system, new_system
    real(8) :: current_energy, current_magnetization, delta_E
    real(8), dimension(3) :: current_magnetization_vector
    
    
    open(12, file='data/stabilizer/'//output_file, status='unknown')
    write(12,*) '#seed - mcs - temperature - nx - ny - initiai_magnetization_vector'
    write(12,*) '#', seed, mcs, T, nx, ny, initial_magnetization_vector(:)
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
            spins_index = random_integer(1, CrAtoms, seed)
            
            
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
                write(12, '(I7 , F16.9, F16.9)') i, current_energy/Cr_atoms , current_magnetization
                flush(12)
            end if
        end do


    end do
    close(12)
    
    print *, 'End of the simulation'
end subroutine simulation


end program main