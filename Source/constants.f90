module constants
    real(8), parameter :: pi = 3.141592653589793238462643383279502884197d0
    real(8), parameter :: h_bar = 6.582119569D-13
    real(8), parameter :: muB = 5.788381756d-2 !meV/T
    real(8), parameter :: kB = 8.617333262D-2 !meV/K  

    integer, save :: Cr_atoms, seed
    real(8), save :: system_energy
    integer, save, allocatable :: neighbors(:,:)
    real(8), save, allocatable :: spins(:,:)
    real(8), save :: J, K, anisotropy, H, g, T
    real(8), save :: anisotropy_vector(3), H_vector(3), mag_vec(3)
    !real(8), save :: avg_energy, avg_mag, avg_energy2, avg_mag2

contains 

function read_neighbors(file, Cr_atoms)  result(neighbors)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(in) :: Cr_atoms
    integer, dimension(:,:), allocatable :: neighbors
    integer :: i
    open(1, file='../neighbors/'//file, status='old')
    allocate(neighbors(Cr_atoms, 3))
    do i = 1, Cr_atoms
            read(1, '(3(I8))') neighbors(i,1), neighbors(i,2), neighbors(i,3)
    end do
end function read_neighbors




end module constants