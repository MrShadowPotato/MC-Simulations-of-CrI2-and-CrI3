
!module metropolis
!    use constants
!    use energy
!    use rng
!   contains

subroutine metropolis_rng
    use rng 
    !use energy
    use constants
    implicit none
    integer :: i, index_spin
    real(8) :: delta_E, old_spin(3), new_spin(3)
    logical :: accept
    
    do i = 1, Cr_atoms
        index_spin = random_integer(1, Cr_atoms, seed)
        old_spin = spins(index_spin, :)
        new_spin = random_normal_vector(seed)
        
        call energy_change(old_spin, new_spin, index_spin, delta_E)


        if (delta_E < 0.0d0) then
            accept = .true.
        else
            accept = exp(-delta_E / (kB * T)) > ran2(seed)
        end if
        if (accept) then
            spins(index_spin, :) = new_spin
            system_energy = system_energy + delta_E
            mag_vec = mag_vec + new_spin(:) - old_spin(:)
            !write(66,*) mag_vec(:), sum(spins(:,1)), sum(spins(:,2)), sum(spins(:,3))

        end if
    end do
    !stop


end subroutine metropolis_rng




!end module metropolis
