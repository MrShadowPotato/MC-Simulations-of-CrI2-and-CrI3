program temperature_iterator
    use variables 
    use rng
    implicit none
    

    ! Define local variables.
    integer :: i
    real :: start_time, end_time, elapsed_time, current_time
    real(8), dimension(:,:), allocatable :: spins
    integer, dimension(:,:), allocatable :: neighbors
  
    call read_parameters()    
    call cpu_time(start_time)

    neighbors = read_neighbors(Cr_neighbors, Cr_atoms)
 
    write(*,*) 'Now we generate the spins'
    spins = generate_spins(Cr_atoms, spins_orientation)
    

    
    call temperature_sim(mcs, index_avg , Cr_atoms, nx, ny, initial_temperature, final_temperature, &
    temperature_step, exchange, 0.0d0, spins, neighbors)
    print *, 'This is a mesagge to check that the program has finished.'

    call cpu_time(end_time) 
    elapsed_time = end_time - start_time

    
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


function energy_change(old_spin, new_spin, neighbors, J, index) result(dE)
    real(8), intent(in) :: old_spin(3), new_spin(3)
    integer, intent(in), dimension(:,:) :: neighbors
    real(8), intent(in) :: J
    integer, intent(in) :: index
    real(8) :: ennew, enold, dE

    ennew = 0.0d0
    enold = 0.0d0

    do i = 1, size(neighbors,2)
        ennew = ennew - J * dot_product(new_spin, spins(neighbors(index,i),:))
        enold = enold - J * dot_product(old_spin, spins(neighbors(index,i),:))
    end do
    dE = (ennew - enold)

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
subroutine temperature_sim(mcs, index_avg, CrAtoms, nx, ny, initial_temperature, final_temperature, &
    temperature_step, J, H, spins, neighbors) 
    implicit none
    
    !Declare input variables
    integer, intent(in) :: mcs, CrAtoms, nx, ny, index_avg
    real(8), intent(in) :: J, H, initial_temperature, final_temperature, temperature_step
    real(8), intent(inout) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    
    !Declare local variables
    integer :: i, k, spins_index
    real(8) :: current_energy, current_magnetization, delta_E, T, &
    energy_sum, magnetization_sum, energy_avg, magnetization_avg, &
    sqrd_energy_sum, sqrd_magnetization_sum, sqrd_energy_avg, sqrd_magnetization_avg, &
    beta, Cv, Chi
    real(8), dimension(3) :: current_magnetization_vector, old_spin, new_spin
    
    current_energy = calculate_initial_energy(spins, neighbors, J)
    current_magnetization_vector = calculate_initial_magnetization_vector(spins)
    current_magnetization = sqrt(dot_product(current_magnetization_vector, current_magnetization_vector))
    T = initial_temperature

    open(12, file='data/temperature/'//temp_iterator_file, status='unknown')
    write(12,*) '#seed - mcs - indavg - iT - fT - dT - CrAtoms - nx - ny - H'
    write(12,*) '#', seed, mcs, index_avg, initial_temperature, &
    final_temperature, temperature_step, CrAtoms, nx, ny, H
    write(12,97) 
    97 format('#      T',7x,'|mag|',6x,'magx',6x,'magy',6x,'magz',&
        7x,'ener',3x,'Cv',8x,'suscep') 
    flush(12)
    do while (T < final_temperature)
        
        do i = 1, mcs
            !Metropolis
            !Loop over the number of CrAtoms
            do k = 1, CrAtoms
                !Generate the index of a random spin to be flipped
                spins_index = random_integer(1, CrAtoms, seed)
                
                
                !Store the old spin and generate a new one.
                old_spin = spins(spins_index,:)
                new_spin = random_normal_vector(seed)

                !Calculate the energy and magnetization of the old and new system
                
                delta_E = energy_change(old_spin, new_spin, neighbors, J, spins_index)

                !Decide wheter to accept the new system or not
                if (accept_change(delta_E, T)) then
                    !If the new system is accepted, update the spins and the total energy
                    spins(spins_index,:) = new_spin
                    current_energy = current_energy + delta_E
                    current_magnetization_vector = current_magnetization_vector &
                    - old_spin(:) + new_spin(:)
                end if
                current_magnetization = sqrt(dot_product(current_magnetization_vector, current_magnetization_vector))
                current_magnetization = current_magnetization / CrAtoms
            end do
        end do
        energy_sum = 0
        sqrd_energy_sum = 0
        magnetization_sum = 0
        sqrd_magnetization_sum = 0
        do i = 1, index_avg
            !write(6,*) 'Temperature: ', T, 'Index avg: ', i
            !flush(6)
            do k = 1, CrAtoms
                !Generate the index of a random spin to be flipped
                spins_index = random_integer(1, CrAtoms, seed)
                
                
                !Store the old spin and generate a new one.
                old_spin = spins(spins_index,:)
                new_spin = random_normal_vector(seed)
                
                !Calculate the energy and magnetization of the old and new system
                
                delta_E = energy_change(old_spin, new_spin, neighbors, J, spins_index)
                !Decide wheter to accept the new system or not
                if (accept_change(delta_E, T)) then
                    !If the new system is accepted, update the spins and the total energy
                    spins(spins_index,:) = new_spin 
                    current_energy = current_energy + delta_E
                    current_magnetization_vector = current_magnetization_vector &
                    - old_spin(:) + new_spin(:)
                end if
                current_magnetization = sqrt(dot_product(current_magnetization_vector, current_magnetization_vector))
                current_magnetization = current_magnetization / CrAtoms
            end do
            energy_sum = energy_sum + current_energy
            magnetization_sum = magnetization_sum + current_magnetization
            sqrd_energy_sum = sqrd_energy_sum + current_energy**2
            sqrd_magnetization_sum = sqrd_magnetization_sum + current_magnetization**2
        end do
        energy_avg = energy_sum / index_avg
        magnetization_avg = magnetization_sum / index_avg
        sqrd_energy_avg = sqrd_energy_sum / index_avg
        sqrd_magnetization_avg = sqrd_magnetization_sum / index_avg
        beta = 1 / T*kB
       
        Cv = beta**2 * (sqrd_energy_avg - energy_avg**2)/(Cr_atoms**2) * kB
        Chi =  (sqrd_magnetization_avg - magnetization_avg**2) * beta
        write(12,34) T, magnetization_avg, current_magnetization_vector(:)/Cr_atoms, energy_avg/Cr_atoms, Cv, Chi
        34 format(F12.8,4(1x,F9.6),1x,F14.8,1x,E16.9,1x,E16.9)
        flush(12)
        
        T = T + temperature_step
    end do
    close(12)
    print *, 'End of the temperature_sim'
end subroutine temperature_sim




end program temperature_iterator