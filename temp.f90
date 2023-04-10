program temperature_iterator
    use variables 
    !use rng
    implicit none
    

    ! Define local variables.
    integer :: i
    real :: start_time, end_time, elapsed_time
    real(8), dimension(:,:), allocatable :: spins
    integer, dimension(:,:), allocatable :: neighbors
  
    !Call variables from the vairables module.
    call read_parameters()


    ! Calculate the number of atoms in the structure.
    
    call cpu_time(start_time)

    neighbors = read_neighbors(Cr_neighbors, Cr_atoms)
 
    write(*,*) 'Now we generate the spins'
    spins = generate_spins(Cr_atoms, spins_orientation)
    

    
    call temperature_sim(mcs, index_avg , Cr_atoms, nx, ny, initial_temperature, final_temperature, &
    temperature_step, exchange, 0.0d0, spins, neighbors)
    print *, 'This is a mesagge to check that the program has finished.'

    call cpu_time(end_time) 
    elapsed_time = end_time - start_time

    open(unit=14, file='data/temperature/'//temp_iterator_file, status="old", position='append', action="write")
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
            spins(i,:) = random_normal_vector()
        end do
        
    else 
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
    real(8), dimension(CrAtoms, 3) :: old_system, new_system
    real(8) :: current_energy, current_magnetization, delta_E, T, &
    energy_sum, magnetization_sum, energy_avg, magnetization_avg
    real(8), dimension(3) :: current_magnetization_vector
    
    current_energy = calculate_initial_energy(spins, neighbors, J)
    current_magnetization_vector = calculate_initial_magnetization_vector(spins)
    current_magnetization = sqrt(dot_product(current_magnetization_vector, current_magnetization_vector))
    T = initial_temperature


    open(12, file='data/temperature/'//temp_iterator_file, status='unknown')
    write(12,*) '#seed - mcs - indavg - iT - fT - dT - CrAtoms - nx - ny - H'
    write(12,*) '#', seed, mcs, index_avg, initial_temperature, &
    final_temperature, temperature_step, CrAtoms, nx, ny, H
    write(12,*) '#T', 'Energy', 'Magnetization'
    
    do while (T < final_temperature)
        print *, 'Temperature: ', T
        !Loop over the number of Monte Carlo stepsd
        do i = 1, mcs
            print *, 'Temperature:  ', T, 'Mcs: ', i
            !Loop over the number of CrAtoms
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
            end do
        end do
        energy_sum = 0
        magnetization_sum = 0
        do i = 1, index_avg
            print * ,'Temperature: ', T, 'Index avg: ', i
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
            end do
            energy_sum = energy_sum + current_energy
            magnetization_sum = magnetization_sum + current_magnetization
        end do
        energy_avg = energy_sum / index_avg
        magnetization_avg = magnetization_sum / index_avg
        write(12,*) T, energy_avg, magnetization_avg
        T = T + temperature_step
    end do
    close(12)
    print *, 'End of the temperature_sim'
end subroutine temperature_sim




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




end program temperature_iterator