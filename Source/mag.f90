module mag 
    use constants
    use rng
    implicit none

    contains 

    function generate_spins(natoms, choice, initial_direction) result(spins)
        integer, intent(in) :: natoms, choice
        real(8), dimension(natoms,3):: spins
        real(8), dimension(3), intent(in) :: initial_direction
        integer :: i
        real(8), dimension(3) :: spin
        
        if (choice == 1) then
            !Each spin is a different randomly generated normal vector.
            do i = 1, natoms
                spins(i,:) = random_normal_vector(seed)
            end do
            
        else 
            !All the spins are the same randomly generated normal vector.
            spin = initial_direction(:)
            spins(:,1) = spin(1)
            spins(:,2) = spin(2)
            spins(:,3) = spin(3)
        end if
    end function generate_spins


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
        end function calculate_initial_magnetization_vector

end module mag