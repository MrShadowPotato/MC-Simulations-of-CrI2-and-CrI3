program sLoop
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
    open(12, file='../data/sLoop/s'//output_file, status='unknown')

    
    !Write parameters to screen
    open(13, file='../data/sLoop/zs'//output_file, status='unknown') 
    write(13, 50) compound, seed!; write(13,*)
    50 format(' compound= ', A4,'    seed= ', I7 )   
    write(13, 51) iT
    51 format(' Temperature= ', F7.1)  
    write(13, 53) J, K, H!; write(13,*)
    53 format(' J= ', F8.3,'    K= ', F8.3, '    H= ', F8.3 )
    write(13, 52) mcs!; write(13,*)
    52 format(' mcs= ', I10)
    write(13, 54) Cr_atoms, nx, ny!; write(13,*)
    54 format(' Cr_atoms= ', I10,'    nx= ', I5, '    ny= ', I5 )
    write(13, 55) spins_orientation, magnetization_direction!; write(13,*)
    55 format(' spins_orientation= ', I2,'    magnetization_direction= ', A13 )
    write(13, 56) mag_vec(:)/Cr_atoms!; write(13,*)
    56 format(' initial_magnetization_vector= ', 3(XF8.5)  )
    write(13, 57) H_vector, easy_vector!; write(13,*)
    57 format(' H_vector= ', 3(XF8.5),'    easy_vector= ', 3(XF8.5)  )
    write(13, 58) dS
    58 format(' dS= ', F3.1  )
    flush(13)

    do i = 1, mcs
        call metropolis_rng
        write(12,34) i, sqrt(dot_product(mag_vec, mag_vec))/Cr_atoms, system_energy/Cr_atoms, &
        mag_vec(:)/Cr_atoms 
        34 format(I10,1x,F14.8,1x,F14.8,3(1x,F9.6))
        flush(12)
    end do

    close(12)

    call cpu_time(end_time)
    time_elapsed = end_time - start_time
    
    write(13, *) "The program took:", time_elapsed, " seconds to run."
    close(13)
end program sLoop
