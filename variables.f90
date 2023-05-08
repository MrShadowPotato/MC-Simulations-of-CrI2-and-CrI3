module variables
    implicit none
    character(len=1000) :: dummy, magnetization_direction
    character(len=100), public :: output_file, Cr_xyz, Cr_neighbors, temp_iterator_file, compound, &
    hLoop_file
    integer, public :: seed, mcs, nx, ny, Cr_atoms, spins_orientation, index_avg, max_neighbors
    real(8), public :: temperature, initial_temperature, final_temperature, temperature_step, neighbor_max_distance, &
    iH, dH
    real(8), dimension(2, 3) :: primitive, vlattice
    real(8), dimension(2, 3) :: primitiveCrI3, primitiveCrI2, vlatticeCrI3, vlatticeCrI2
    real(8), dimension(:,:), allocatable :: basis
    real(8), dimension(8, 3) :: basis_CrI3 ! Basis vectors for the CrI3 structure.
    real(8), dimension(6, 3) :: basis_CrI2 ! Basis vectors for the CrI2 structure.
    real(8), dimension(3) :: easy_vector, easy_vectorCrI3, easy_vectorCrI2, H_vector !Easy axis vector.
    real(8), dimension(3), public :: initial_magnetization_vector
    real(8), public :: kB = 8.617333262D-2
    real(8), public :: exchangeCrI3, exchangeCrI2, exchange
    real(8), parameter :: anisotropy = 0.67 ! x2???
    real(8), parameter :: g = 3.8
    real(8), parameter :: h_bar = 6.582119569D-13
    real(8), parameter :: muB = 5.7883818060D-2
    character(len=2), dimension(8) :: elementsCrI3 = (/ 'Cr', 'Cr', 'I ', 'I ', 'I ', 'I ', 'I ', 'I ' /)
    character(len=2), dimension(6) :: elementsCrI2 = (/ 'Cr', 'Cr', 'I ', 'I ', 'I ', 'I '/)
    character(len=2), dimension(:), allocatable :: elements
    
contains

    subroutine read_parameters()
        ! Open the file
        character(len=100) :: filename
        integer :: status
        filename = "params.txt"
        open(unit=10, file=filename, status="old", action="read", iostat=status)
        
        ! Read the parameters
        read(10,*) dummy, dummy, seed
        read(10,*) dummy, dummy, mcs
        read(10,*) dummy, dummy, temperature
        read(10,*) dummy, dummy, nx
        read(10,*) dummy, dummy, ny
        read(10,*) dummy, dummy, spins_orientation
        read(10,*) dummy, dummy, output_file
        read(10,*) dummy, dummy, Cr_xyz
        read(10,*) dummy, dummy, Cr_neighbors
        read(10,*) dummy, dummy, initial_temperature
        read(10,*) dummy, dummy, final_temperature
        read(10,*) dummy, dummy, temperature_step
        read(10,*) dummy, dummy, index_avg
        read(10,*) dummy, dummy, temp_iterator_file
        read(10,*) dummy, dummy, compound
        read(10,*) dummy, dummy, magnetization_direction
        read(10,*) dummy, dummy, exchangeCrI3
        read(10,*) dummy, dummy, exchangeCrI2
        read(10,*) dummy, dummy, iH
        read(10,*) dummy, dummy, dH
        read(10,*) dummy, dummy, hLoop_file

        
        ! Close the file
        close(unit=10)
        

        Cr_atoms = 2*nx*ny

        ! Define the primitive vectors for the CrI3 and CrI2 structure.
        !--> CrI3
        primitiveCrI3(1, :) = [6.8911914825, 0.0, 0.0]
        primitiveCrI3(2, :) = [-3.4488059462999998, 5.9711638719, 0.0]
        !--> CrI2
        primitiveCrI2(1, :) = [3.9350889407837069, 0.0, 0.0]
        primitiveCrI2(2, :) = [0.0, 7.793970515734439, -0.0155822353114624]


        !Define the lattice vectors for periodic boundary conditions.
        !--> CrI3
        vlatticeCrI3(1,:) = [primitiveCrI3(1,1)*nx, primitiveCrI3(1,2)*nx, primitiveCrI3(1,3)*nx]
        vlatticeCrI3(2,:) = [primitiveCrI3(2,1)*ny, primitiveCrI3(2,2)*ny, primitiveCrI3(2,3)*ny]
        !--> CrI2
        vlatticeCrI2(1,:) = [primitiveCrI2(1,1)*nx, primitiveCrI2(1,2)*nx, primitiveCrI2(1,3)*nx]
        vlatticeCrI2(2,:) = [primitiveCrI2(2,1)*ny, primitiveCrI2(2,2)*ny, primitiveCrI2(2,3)*ny]
       
       
        ! Define the easy axis vector for CrI3 and CrI2.
        !--> CrI3
        easy_vectorCrI3(1) = primitiveCrI3(1,2) * primitiveCrI3(2,3) - primitiveCrI3(1,3) * primitiveCrI3(2,2)
        easy_vectorCrI3(2) = primitiveCrI3(1,3) * primitiveCrI3(2,1) - primitiveCrI3(1,1) * primitiveCrI3(2,3)
        easy_vectorCrI3(3) = primitiveCrI3(1,1) * primitiveCrI3(2,2) - primitiveCrI3(1,2) * primitiveCrI3(2,1)
        easy_vectorCrI3 = easy_vectorCrI3 / sqrt(dot_product(easy_vectorCrI3, easy_vectorCrI3))
        !--> CrI2
        easy_vectorCrI2(1) = primitiveCrI2(1,2) * primitiveCrI2(2,3) - primitiveCrI2(1,3) * primitiveCrI2(2,2)
        easy_vectorCrI2(2) = primitiveCrI2(1,3) * primitiveCrI2(2,1) - primitiveCrI2(1,1) * primitiveCrI2(2,3)
        easy_vectorCrI2(3) = primitiveCrI2(1,1) * primitiveCrI2(2,2) - primitiveCrI2(1,2) * primitiveCrI2(2,1)
        easy_vectorCrI2 = easy_vectorCrI2 / sqrt(dot_product(easy_vectorCrI2, easy_vectorCrI2))
       
       
        ! Define the basis vectors for the CrI3 and CrI2 structure.
        !--> CrI3
        basis_CrI3(1,:) = [3.42245459776438, 1.94691491308500, 9.72065830568899]   !Cr
        basis_CrI3(2,:) = [-0.02339545180016, 3.93735766285690, 9.72080135688876]  !Cr
        basis_CrI3(3,:) = [-1.03538548904423, 5.92724513968807, 11.28527832915810] !I 
        basis_CrI3(4,:) = [2.20273590406674, 3.81437945298318, 11.29027271402490]  !I
        basis_CrI3(5,:) = [5.64639472731378, 2.07416105561567, 11.29166984494761]  !I
        basis_CrI3(6,:) = [0.97989374494890, 5.92913055695139, 8.15211676904270]   !I
        basis_CrI3(7,:) = [1.19604098486642, 2.06898570270974, 8.15546418188776]   !I
        basis_CrI3(8,:) = [4.64469384948989, 3.81029176486010, 8.15107346274473]   !I
        !--> CrI2
        basis_CrI2(1,:) = [0.00000000000000, 1.18048232888026, 15.66708540526766]  !Cr
        basis_CrI2(2,:) = [1.96754466301406, 5.07735235387682, 15.66009018229279]  !Cr
        basis_CrI2(3,:) = [1.96754466301406, 2.14455628147048, 17.30559462841125]  !I 
        basis_CrI2(4,:) = [1.96754466301406, 0.21647035484051, 14.02945080633759]  !I
        basis_CrI2(5,:) = [0.00000000000000, 6.04092814648801, 17.29708073114964]  !I
        basis_CrI2(6,:) = [0.00000000000000, 4.11332427420235, 14.02184874359694]  !I

        
        
        if (compound == 'CrI3') then
            primitive = primitiveCrI3
            vlattice = vlatticeCrI3
            easy_vector = easy_vectorCrI3
            H_vector = easy_vector!!!!
            allocate(basis(8,3))
            basis(:,:) = basis_CrI3(:,:)
            max_neighbors = 3
            allocate(elements(8))
            elements(:) = elementsCrI3(:)
            neighbor_max_distance = 4.0d0
            exchange = exchangeCrI3

        else if (compound == 'CrI2') then
            primitive = primitiveCrI2
            vlattice = vlatticeCrI2
            easy_vector = easy_vectorCrI2
            H_vector = easy_vector!!!!!
            allocate(basis(6,3))
            basis(:,:) = basis_CrI2(:,:)
            max_neighbors = 2
            allocate(elements(6))
            elements(:) = elementsCrI2(:)
            neighbor_max_distance = 4.0d0
            exchange = exchangeCrI2
        else 
            write(6,*) "Error: Compound does not match available options."
            write(6,*) "Closing program..."
            stop
        end if
        
        !Define initial magnetization vector
        if (magnetization_direction == 'easy') then
            initial_magnetization_vector = easy_vector

        else if (magnetization_direction == 'reversed_easy') then
            initial_magnetization_vector = -easy_vector

        else if (magnetization_direction == 'p1') then !takes the direction of the first primitive vector
            initial_magnetization_vector = primitive(1,:)

        else if (magnetization_direction == 'p2') then
            initial_magnetization_vector = primitive(2,:)

        else if (magnetization_direction == 'p1_p2plane') then
            initial_magnetization_vector = primitive(1,:) + primitive(2,:)
        else if (magnetization_direction == '45deg') then 
            initial_magnetization_vector = primitive(1,:) + primitive(2,:)
            initial_magnetization_vector = initial_magnetization_vector /  &
            sqrt(dot_product(initial_magnetization_vector, initial_magnetization_vector))
            initial_magnetization_vector = initial_magnetization_vector + easy_vector
        end if
        

        initial_magnetization_vector = initial_magnetization_vector / &
        sqrt(dot_product(initial_magnetization_vector, initial_magnetization_vector))

    end subroutine read_parameters
    

    
end module variables
