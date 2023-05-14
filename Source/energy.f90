
!Calculate the energy of the system consering the interaction between spins, the anisotropy and the Zeeman effect.
subroutine calculate_system_energy()
    use constants    
    implicit none
    real(8) :: ienergy(3)
    integer :: i

    do i = 1, size(spins,1)
        call ienergies(i, ienergy(:))
        system_energy = system_energy  + ienergy(1) * 0.5d0 + ienergy(2) + ienergy(3) 
        !write(6,*) ienergy(:)
    end do

end subroutine calculate_system_energy



subroutine energy_change(old_spin, new_spin, index, delta_E)
    implicit none
    real(8), intent(in), dimension(3) :: old_spin, new_spin 
    integer, intent(in) :: index
    
    
    real(8) :: old_energy(3), new_energy(3)
    real(8), intent(out) :: delta_E
    
    old_energy = 0.0d0
    new_energy = 0.0d0
    call heisenberg_energy(old_spin, index, old_energy(1))
    call anisotropy_energy(old_spin, old_energy(2))
    call zeeman_energy(old_spin, old_energy(3))
    call heisenberg_energy(new_spin, index, new_energy(1))
    call anisotropy_energy(new_spin, new_energy(2))
    call zeeman_energy(new_spin, new_energy(3))
    !delta_E = old_energy(1) + old_energy(2) + old_energy(3) - new_energy(1) - new_energy(2) - new_energy(3)
    delta_E = sum(new_energy) - sum(old_energy)

end subroutine energy_change



subroutine ienergies(index, energies)
    use constants
    implicit none
    integer, intent(in) :: index
    real(8), intent(out) :: energies(3)
    real(8), dimension(3) :: spin

    energies = 0.0d0
    spin = spins(index,:)
    !Energia de interaccion
    call heisenberg_energy(spin, index, energies(1))
    !Energia de anisotropia
    call anisotropy_energy(spin, energies(2))
    !Energia de Zeeman
    call zeeman_energy(spin, energies(3))
end subroutine ienergies

subroutine heisenberg_energy(spin, index, energy)
    use constants
    implicit none
    real(8), intent(in) :: spin(3)
    integer, intent(in) :: index
    real(8), intent(out) :: energy

    integer :: i
    energy = 0.0d0

    !Exchange energy
    do i=1, size(neighbors,2)
        energy = energy - J * dot_product(spin(:), spins(neighbors(index,i),:))
    end do
    !write(6,*) energy, J; stop
end subroutine heisenberg_energy

subroutine anisotropy_energy(spin,energy)
    use constants
    implicit none
    real(8), intent(in) :: spin(3)
    real(8), intent(out) :: energy

    energy = 0.0d0

    !Anisotropy energy
    energy = energy - anisotropy * dot_product(spin(:), anisotropy_vector(:))**2

end subroutine anisotropy_energy

subroutine zeeman_energy(spin, energy)
    use constants
    implicit none
    real(8), intent(in) :: spin(3)
    real(8), intent(out) :: energy

    energy = 0.0d0

    !Zeeman energy
    energy = energy - g * muB * H * dot_product(spin(:), H_vector(:)) !/ h_bar

end subroutine zeeman_energy


