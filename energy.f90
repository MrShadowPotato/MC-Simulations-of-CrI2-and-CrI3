module energy 
    use constants
contains


!Calculate the energy of the system consering the interaction between spins, the anisotropy and the Zeeman effect.
function calculate_system_energy(spins, neighbors, J, anisotropy, anisotropy_vector, H, H_vector, g) result(system_energy)
    implicit none
    real(8), intent(in) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    real(8), intent(in) :: J, anisotropy, H, g
    real(8), intent(in), dimension(3) :: H_vector, anisotropy_vector
    real(8) :: system_energy, ienergy(3)
    integer :: i
    system_energy = 0

    do i = 1, size(spins,1)
        ienergy(:) = ienergies(i, spins, neighbors, J, anisotropy, anisotropy_vector, H, H_vector, g)
        system_energy = system_energy  + ienergy(1) * 0.5d0 + ienergy(2) + ienergy(3) 
    end do

end function calculate_system_energy






function energy_change(old_spin, new_spin, index, spins, neighbors, J, anisotropy, anisotropy_vector, H, H_vector, g) & 
    result(delta_E)
    implicit none
    real(8), intent(in), dimension(3) :: old_spin, new_spin 
    integer, intent(in) :: index
    real(8), intent(in) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    real(8), intent(in) :: J, anisotropy, H, g
    real(8), intent(in), dimension(3) :: H_vector, anisotropy_vector
    real(8) :: old_energy(3), new_energy(3), delta_E

    old_energy = 0.0d0
    new_energy = 0.0d0
    old_energy(1) = heisenberg_energy(old_spin, index, spins, neighbors, J)
    old_energy(2) = anisotropy_energy(old_spin, anisotropy, anisotropy_vector)
    old_energy(3) = zeeman_energy(old_spin, H, H_vector, g)
    new_energy(1) = heisenberg_energy(new_spin, index, spins, neighbors, J)
    new_energy(2) = anisotropy_energy(new_spin, anisotropy, anisotropy_vector)
    new_energy(3) = zeeman_energy(new_spin, H, H_vector, g)
    delta_E = sum(old_energy) - sum(new_energy)

end function energy_change



function ienergies(index, spins, neighbors, J, anisotropy, anisotropy_vector, H, H_vector, g) result(energies)
    implicit none
    integer, intent(in) :: index
    real(8), intent(in) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    real(8), intent(in) :: J, anisotropy, H, g
    real(8), intent(in), dimension(3) :: H_vector, anisotropy_vector
    real(8), dimension(3) :: energies, spin

    energies = 0.0d0
    spin = spins(index,:)
    !Energia de interaccion
    energies(1) = heisenberg_energy(spin, index, spins, neighbors, J)

    !Energia de anisotropia
    energies(2) = anisotropy_energy(spin, anisotropy, anisotropy_vector)
    !Energia de Zeeman
    energies(3) = zeeman_energy(spin, H, H_vector, g)
end function ienergies

function heisenberg_energy(spin, index, spins, neighbors, J) result(energy)
    implicit none
    real(8), intent(in) :: spin(3)
    integer, intent(in) :: index
    real(8), intent(in) :: spins(:,:)
    integer, intent(in), dimension(:,:) :: neighbors
    real(8), intent(in) :: J
    real(8) :: energy

    integer :: k
    energy = 0.0d0

    !Exchange energy
    do k=1, size(neighbors,2)
        energy = energy - J * dot_product(spin(:), spins(neighbors(index,k),:))
    end do
end function heisenberg_energy

function anisotropy_energy(spin, anisotropy, anisotropy_vector) result(energy)
    implicit none
    real(8), intent(in) :: spin(3)
    real(8), intent(in) :: anisotropy
    real(8), intent(in), dimension(3) :: anisotropy_vector
    real(8) :: energy

    energy = 0.0d0

    !Anisotropy energy
    energy = energy - anisotropy * dot_product(spin(:), anisotropy_vector(:))**2

end function anisotropy_energy

function zeeman_energy(spin, H, H_vector, g) result(energy)
    implicit none
    real(8), intent(in) :: spin(3)
    real(8), intent(in) :: H, g
    real(8), intent(in), dimension(3) :: H_vector
    real(8) :: energy

    energy = 0.0d0

    !Zeeman energy
    energy = energy - g * muB * H * dot_product(spin(:), H_vector(:)) !/ h_bar

end function zeeman_energy



end module energy