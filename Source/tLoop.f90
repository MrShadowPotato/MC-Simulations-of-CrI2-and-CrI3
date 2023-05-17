program tLoop
    use read_input
    use constants
    use mag
    implicit none

    real(8) :: start_time, end_time, time_elapsed
    real(8) :: beta, Cv, Chi
    real(8) :: avg_energy, avg_mag, avg_energy2, avg_mag2
    real(8) :: avg_mag_vec(3)
    integer :: i




    call cpu_time(start_time)
    call read_parameters()


    neighbors = read_neighbors(Cr_atoms)
    spins = generate_spins(Cr_atoms, spins_orientation, initial_magnetization_vector)
    T = iT
    H = iH
    call calculate_system_energy
    mag_vec = calculate_initial_magnetization_vector(spins)
    write(6, *) "Initial energy: ", system_energy/Cr_atoms
    write(6, *) "Initial magnetization: ", mag_vec(:)/Cr_atoms
    open(12, file='../data/tLoop/'//output_file, status='unknown')

    
    !Write parameters to screen
    open(13, file='../data/tLoop/z'//output_file, status='unknown') 
    write(13, 50) compound, seed!; write(13,*)
    50 format(/,' compound= ', A4,'    seed= ', I7 /)
    write(13, 51) iT, fT, dT!; write(13,*)
    51 format(/,' iT= ', F7.1,'    fT= ', F7.1, '    dT= ', F7.1  /)  
    write(13, 53) J, K, H!; write(13,*)
    53 format(/,' J= ', F8.3,'    K= ', F8.3, '    H= ', F8.3  /)
    write(13, 52) mcs, neq!; write(13,*)
    52 format(/,' mcs= ', I10,'    neq= ', I10 /)
    write(13, 54) Cr_atoms, nx, ny!; write(13,*)
    54 format(/,' Cr_atoms= ', I10,'    nx= ', I5, '    ny= ', I5  /)
    write(13, 55) spins_orientation, magnetization_direction!; write(13,*)
    55 format(/,' spins_orientation= ', I2,'    magnetization_direction= ', A13  /)
    write(13, 56) mag_vec(:)/Cr_atoms!; write(13,*)
    56 format(/,' initial_magnetization_vector= ', 3(XF8.5)  /)
    write(13, 57) H_vector, easy_vector!; write(13,*)
    57 format(/,' H_vector= ', 3(XF8.5),'    easy_vector= ', 3(XF8.5)  /)
    flush(13)

    do while (T < fT)
        do i = 1, mcs
            call metropolis_rng
        end do
        avg_energy = 0
        avg_mag = 0
        avg_energy2 = 0
        avg_mag2 = 0
        avg_mag_vec = 0
        do i = 1, neq
            call metropolis_rng
            avg_energy = avg_energy + system_energy
            avg_energy2 = avg_energy2 + system_energy**2
            avg_mag = avg_mag + sqrt(dot_product(mag_vec, mag_vec))
            avg_mag2 = avg_mag2 + dot_product(mag_vec, mag_vec)
            avg_mag_vec = avg_mag_vec + mag_vec
        end do
        avg_energy = avg_energy / neq
        avg_energy2 = avg_energy2 / neq
        avg_mag = avg_mag / neq
        avg_mag2 = avg_mag2 / neq
        avg_mag_vec = avg_mag_vec / neq
        
        beta = 1 / T*kB
        Cv = beta**2 * (avg_energy2 - avg_energy**2)
        Chi = beta * (avg_mag2 - avg_mag**2)
        write(12,34) T, avg_mag/Cr_atoms, avg_mag_vec(:)/Cr_atoms, avg_energy/Cr_atoms, Cv, Chi
        34 format(F12.8,4(1x,F9.6),1x,F14.8,1x,E16.9,1x,E16.9)
        flush(12)
        T = T + dT
    end do
    close(12)

    call cpu_time(end_time)
    time_elapsed = end_time - start_time
    
    write(13, *) "The program took:", time_elapsed, " seconds to run."
    close(13)
end program tLoop
