program tLoop2
    use read_input
    use constants
    !use energy
    !use metropolis
    use mag
    implicit none

    real(8) :: start_time, end_time, beta, Cv, Chi
    real(8) :: avg_energy, avg_mag, avg_energy2, avg_mag2
    real(8) :: avg_mag_vec(3)
    integer :: i


    call cpu_time(start_time)
    call read_parameters()
    neighbors = read_neighbors(Cr_neighbors, Cr_atoms)
    spins = generate_spins(Cr_atoms, spins_orientation, initial_magnetization_vector)
    T = initial_temperature
    call calculate_system_energy
    mag_vec = calculate_initial_magnetization_vector(spins)
    write(6, *) "Initial energy: ", system_energy/Cr_atoms
    write(6, *) "Initial magnetization: ", mag_vec(:)/Cr_atoms
    !stop
    open(12, file='output.txt', status='unknown')
    do while (T < final_temperature)
        do i = 1, mcs
            call metropolis_rng
        end do
        avg_energy = 0
        avg_mag = 0
        avg_energy2 = 0
        avg_mag2 = 0
        avg_mag_vec = 0
        do i = 1, index_avg
            call metropolis_rng
            avg_energy = avg_energy + system_energy
            avg_energy2 = avg_energy2 + system_energy**2
            avg_mag = avg_mag + sqrt(dot_product(mag_vec, mag_vec))
            avg_mag2 = avg_mag2 + dot_product(mag_vec, mag_vec)
            avg_mag_vec = avg_mag_vec + mag_vec
        end do
        avg_energy = avg_energy / index_avg
        avg_energy2 = avg_energy2 / index_avg
        avg_mag = avg_mag / index_avg
        avg_mag2 = avg_mag2 / index_avg
        avg_mag_vec = avg_mag_vec / index_avg
        
        beta = 1 / T*kB
        Cv = beta**2 * (avg_energy2 - avg_energy**2)
        Chi = beta * (avg_mag2 - avg_mag**2)
        T = T + temperature_step
        write(12,34) T, avg_mag/Cr_atoms, avg_mag_vec(:)/Cr_atoms, avg_energy/Cr_atoms, Cv, Chi
        34 format(F12.8,4(1x,F9.6),1x,F14.8,1x,E16.9,1x,E16.9)
        flush(12)
    end do
    call cpu_time(end_time)


end program tLoop2
