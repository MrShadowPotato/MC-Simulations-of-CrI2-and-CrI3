program hLoop
    use read_input
    use constants
    use mag
    implicit none
    real(8) :: start_time, end_time, time_elapsed
    real(8) :: MH, avg_MH, avg_energy, avg_mag, avg_mag_vec(3)
    integer :: i, i2, index
    logical :: first_time = .true.
    real(8), allocatable, dimension(:,:)  :: spins_ma, spins_mb
    real(8) :: Ma, Mb, Ms, Ma_vec(3), Mb_vec(3), Ms_vec(3)

    call cpu_time(start_time)
    call read_parameters()

    neighbors = read_neighbors(Cr_atoms, max_neighbors)
    if (compound == 'CrI3') then 
        spins = generate_spins(Cr_atoms, spins_orientation, initial_magnetization_vector)
    else
        allocate(spins_ma(Cr_atoms/2, 3))
        allocate(spins_mb(Cr_atoms/2, 3))
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

    t = iT
    H = iH
    call calculate_system_energy
    mag_vec = calculate_initial_magnetization_vector(spins)
    write(6, *) "Initial energy: ", system_energy/Cr_atoms
    write(6, *) "Initial magnetization: ", mag_vec(:)/Cr_atoms
    open(12, file='../data/hLoop/h'//output_file, status='unknown')
    if (compound == 'CrI3') then 
        MH = dot_product(mag_vec, H_vector(:))/Cr_atoms
    else
        MH = dot_product(Ms_vec, H_vector(:))
    end if

    
    !Write parameters to screen
    open(13, file='../data/hLoop/zh'//output_file, status='unknown') 
    write(13, 50) compound, seed
    50 format(' compound= ', A4,'    seed= ', I7 )
    write(13, 51) iT
    51 format(' Temperature= ', F7.1)  
    write(13, 58) iH, dH
    58 format(' iH= ', F8.3,'    dH= ', F8.3)
    write(13, 53) J, K
    53 format(' J= ', F8.3,'    K= ', F8.3)
    write(13, 52) mcs, neq
    52 format(' mcs= ', I10,'    neq= ', I10 )
    write(13, 54) Cr_atoms, nx, ny
    54 format(' Cr_atoms= ', I10,'    nx= ', I5, '    ny= ', I5  )
    write(13, 55) spins_orientation, magnetization_direction!
    55 format(' spins_orientation= ', I2,'    magnetization_direction= ', A13  )
    write(13, 56) mag_vec(:)/Cr_atoms
    56 format(' initial_magnetization_vector= ', 3(XF8.5)  )
    write(13, 57) H_vector, easy_vector
    57 format(' H_vector= ', 3(XF8.5),'    easy_vector= ', 3(XF8.5)  )
    write(13, 59) dS, g
    59 format(' dS= ', F3.1, '    g= ', F3.1 )
    flush(13)
    
    
    write(12, 60) '#H', 'MH', 'dH'
    60 format(A8, 2x, A12, 2x, A8)
    
    do while (H <= iH)
        if (((MH < -0.98).and.first_time).or. (H < -30)) then 
            dH = -abs(dH)
            first_time = .false.
        end if
        do i = 1, mcs
            call metropolis_rng
        end do
        avg_energy = 0
        avg_mag = 0
        avg_mag_vec = 0
        avg_MH = 0
        do i = 1, neq
            call metropolis_rng
            avg_energy = avg_energy + system_energy
            if (compound == 'CrI3') then 
                avg_mag = avg_mag + sqrt(dot_product(mag_vec, mag_vec))
                avg_mag_vec = avg_mag_vec + mag_vec
                
                MH = dot_product(mag_vec, H_vector(:))/Cr_atoms
                avg_MH = avg_MH + MH
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
                avg_mag = avg_mag + sqrt(dot_product(Ms_vec, Ms_vec))
                avg_mag_vec = avg_mag_vec + Ms_vec

                MH = dot_product(Ms_vec, H_vector(:))
                avg_MH = avg_MH + MH
            end if
        end do  
        write(12,33) H, avg_MH/neq, dh!avg_energy/neq/Cr_atoms, avg_mag/Cr_atoms/neq
        33 format(F8.3, F12.8, F8.3)
        flush(12)
        H = H - dH
    end do
    close (12)

    call cpu_time(end_time)
    time_elapsed = end_time - start_time
    write(13, *) "The program took:", time_elapsed, " seconds to run."
    close(13)


end program hLoop   


