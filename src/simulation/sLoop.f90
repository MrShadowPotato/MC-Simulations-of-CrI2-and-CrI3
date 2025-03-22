program sLoop
    use read_input
    use constants
    use mag
    implicit none

    real(8) :: start_time, end_time, time_elapsed
    integer :: i, i2, index
    real(8), allocatable, dimension(:,:)  :: spins_ma, spins_mb
    real(8) :: Ma, Mb, Ms, Ma_vec(3), Mb_vec(3), Ms_vec(3)
    


    call cpu_time(start_time)
    call read_parameters()


    neighbors = read_neighbors(Cr_atoms, max_neighbors)
    if (compound == 'CrI3') then 
        spins = generate_spins(Cr_atoms, spins_orientation, initial_magnetization_vector)
    else
        allocate(spins_ma(Cr_atoms/2,3))
        allocate(spins_mb(Cr_atoms/2,3))
        spins_ma = generate_spins(Cr_atoms/2, spins_orientation, initial_magnetization_vector(:))
        spins_mb = generate_spins(Cr_atoms/2, spins_orientation, -1*initial_magnetization_vector(:))
        allocate(spins(Cr_atoms,3))
        do i = 1, Cr_atoms/2
            index = 2*i - 1
            spins(index,:) = spins_ma(i,:)
            Ma_vec = Ma_vec + spins_ma(i,:)
            index = 2*i
            spins(index,:) = spins_mb(i,:)
            Mb_vec = Mb_vec + spins_mb(i,:)
            
        end do
        Ma_vec = Ma_vec/Cr_atoms*2
        Mb_vec = Mb_vec/Cr_atoms*2
        Ms_vec = (Ma_vec/2 - Mb_vec/2)
        write(6, *)'Ma:  ', sqrt(dot_product(Ma_vec, Ma_vec))!/Cr_atoms/2
        write(6, *)'Mb:  ', sqrt(dot_product(Mb_vec, Mb_vec))!/Cr_atoms/2
        write(6, *)'Ms:  ', sqrt(dot_product(Ms_vec, Ms_vec))!/Cr_atoms
        
    end if
    

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
    write(13, 58) dS, g
    58 format(' dS= ', F3.1, '    g= ', F3.1 )

    flush(13)

    if (compound == 'CrI3') then
        write(12, 37) '#MCS', 'M', 'E', 'Mx', 'My', 'Mz'
        37 format(A10, 1x, 5(A14, 1x))

    else 
        write(12, 35) '#MCS', 'M', 'Ms', 'Ma', 'Mb', 'E'
        35 format(A10, 1x, 5(A14, 1x)) 
    end if 


    do i = 1, mcs
        call metropolis_rng
        if (compound == 'CrI3') then 
            write(12,34) i, sqrt(dot_product(mag_vec, mag_vec))/Cr_atoms, system_energy/Cr_atoms, &
            mag_vec(:)/Cr_atoms 
            34 format(I10,1x,F14.8,1x,F14.8,3(1x,F14.6))
        else 
            Ma_vec = 0
            Mb_vec = 0
            do i2 = 1, Cr_atoms/2
                index = 2*i2 - 1
                spins_ma(i2,:) = spins(index,:) 
                Ma_vec = Ma_vec + spins_ma(i2,:)
                index = 2*i2
                spins_mb(i2,:) = spins(index,:) 
                Mb_vec = Mb_vec + spins_mb(i2,:)
            end do
            Ma_vec = Ma_vec/Cr_atoms*2
            Mb_vec = Mb_vec/Cr_atoms*2
            Ms_vec = (Ma_vec/2 - Mb_vec/2)
            Ma = sqrt(dot_product(Ma_vec, Ma_vec))
            Mb = sqrt(dot_product(Mb_vec, Mb_vec))
            Ms = sqrt(dot_product(Ms_vec, Ms_vec))
            
            write(12,36) i, sqrt(dot_product(mag_vec, mag_vec))/Cr_atoms, Ms, Ma, Mb, system_energy/Cr_atoms 
            36 format(I10, 1x, 4(F14.8, 1x), F14.8)
        end if 

        flush(12)
    end do

    close(12)

    call cpu_time(end_time)
    time_elapsed = end_time - start_time
    
    write(13, *) "The program took:", time_elapsed, " seconds to run."
    close(13)
end program sLoop
