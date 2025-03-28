program hLoop
    use read_input
    use constants
    use mag
    implicit none
    real(8) :: start_time, end_time, time_elapsed
    real(8) :: Ma_angle, Mb_angle, canting_angle
    real(8) :: canting_angle_avg, canting_angle_avg_f, canting_angle_avg_af
    real(8) :: avg_energy, avg_mag, avg_mag_vec(3)
    real(8) :: MH, avg_MH
    integer :: i, i2, index
    logical :: first_time = .true.
    real(8), allocatable, dimension(:,:)  :: spins_ma, spins_mb
    real(8) :: Ma_vec(3), Mb_vec(3), Ms_vec(3)

    MH = 0.0d0  ! Initialize MH to prevent uninitialized use warning

    call cpu_time(start_time)
    call read_parameters()

    neighbors = read_neighbors(Cr_atoms, max_neighbors)
    if (compound == 'CrI3') then 
        spins = generate_spins(Cr_atoms, spins_orientation, initial_magnetization_vector)
    else if (compound == 'CrI2') then
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

        Ma_angle = dot_product(Ma_vec, easy_vector(:))/sqrt(dot_product(Ma_vec, Ma_vec))
        Ma_angle = max(-1.0, min(1.0, Ma_angle))
        Ma_angle = acos(Ma_angle)

        Mb_angle = dot_product(Mb_vec, - easy_vector(:))/sqrt(dot_product(Mb_vec, Mb_vec))
        Mb_angle = max(-1.0, min(1.0, Mb_angle))
        Mb_angle = acos(Mb_angle)
        !Ma_angle = acos(dot_product(Ma_vec, easy_vector(:))/sqrt(dot_product(Ma_vec, Ma_vec)))
        !cos_angle = max(-1.0, min(1.0, Ma_angle))
        !Mb_angle = acos(dot_product(Mb_vec, - easy_vector(:))/sqrt(dot_product(Mb_vec, Mb_vec)))
        canting_angle = Ma_angle/2 + Mb_angle/2
        write(6,*) 'Easy axis: ', easy_vector(:)
        write(6, *)'|easy|', sqrt(dot_product(easy_vector, easy_vector))
        write(6, *)'Ma:  ', sqrt(dot_product(Ma_vec, Ma_vec))!/Cr_atoms/2
        write(6, *)'Mb:  ', sqrt(dot_product(Mb_vec, Mb_vec))!/Cr_atoms/2
        write(6, *)'Mav:  ', Ma_vec(:)
        write(6, *)'Ma:   ', sqrt(dot_product(Ma_vec, Ma_vec))
        write(6, *)'Mbv:  ', Mb_vec(:)
        write(6, *)'Mb:   ', sqrt(dot_product(Mb_vec, Mb_vec))
        write(6,*) 'Ma_angle: ', Ma_angle !* 180 / pi
        write(6,*) 'Mb_angle: ', Mb_angle !* 180 / pi
        write(6,*) 'canting_angle: ', canting_angle !* 180 / pi

        write(6,*) 'Theta Ma', acos(Ma_vec(3)) * 180 / pi
        write(6,*) 'Phi Ma', atan2(Ma_vec(2), Ma_vec(1)) * 180 / pi
        write(6,*) 'Theta Mb', acos(Mb_vec(3)) * 180 / pi
        write(6,*) 'Phi Mb', atan2(Mb_vec(2), Mb_vec(1)) * 180 / pi
        write(6,*) 'Theta Easy' , acos(easy_vector(3)) * 180 / pi
        write(6,*) 'Phi Easy', atan2(easy_vector(2), easy_vector(1)) * 180 / pi
        !stop 

    end if

    t = iT
    H = iH
    call calculate_system_energy
    mag_vec = calculate_initial_magnetization_vector(spins)
    write(6, *) "Initial energy: ", system_energy/Cr_atoms
    open(12, file='../data/hLoop/h'//output_file, status='unknown')
    if (compound == 'CrI3') then 
        write(6, *) "Initial magnetization: ", mag_vec(:)/Cr_atoms
        MH = dot_product(mag_vec, H_vector(:))/Cr_atoms
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
    
  
  
    
    if (compound == 'CrI3') then 
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
                avg_mag = avg_mag + sqrt(dot_product(mag_vec, mag_vec))
                avg_mag_vec = avg_mag_vec + mag_vec
                MH = dot_product(mag_vec, H_vector(:))/Cr_atoms
                avg_MH = avg_MH + MH
            end do  
            write(12,33) H, avg_MH/neq, dh!avg_energy/neq/Cr_atoms, avg_mag/Cr_atoms/neq
            33 format(F8.3, F12.8, F8.3)
            flush(12)
            H = H - dH
        end do
    else if (compound == 'CrI2') then 
        H = 0
        write(12, 36) '#H', 'CA_avg', 'CA_avg_F', 'CA_avg_AF'
        36 format(A8, 2x, A12, 2x, A12, 2x, A12)
        do while (H < iH)
            do i = 1, mcs
                call metropolis_rng
                
            end do
            canting_angle_avg_f = 0
            canting_angle_avg_af = 0
            canting_angle_avg = 0
            do i = 1, neq
                call metropolis_rng
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

                !write(6,*) 'Theta Ma', acos(Ma_vec(3)) * 180 / pi
                !write(6,*) 'Phi Ma', atan2(Ma_vec(2), Ma_vec(1)) * 180 / pi
                !write(6,*) 'Theta Mb', acos(Mb_vec(3)) * 180 / pi
                !write(6,*) 'Phi Mb', atan2(Mb_vec(2), Mb_vec(1)) * 180 / pi
                !write(6,*) 'Theta Easy' , acos(easy_vector(3)) * 180 / pi
                !write(6,*) 'Phi Easy', atan2(easy_vector(2), easy_vector(1)) * 180 / pi
                !stop 
                Ma_angle = dot_product(Ma_vec, easy_vector(:))/sqrt(dot_product(Ma_vec, Ma_vec))
                Ma_angle = max(-1.0, min(1.0, Ma_angle))
                Ma_angle = acos(Ma_angle)

                Mb_angle = dot_product(Mb_vec, -easy_vector(:))/sqrt(dot_product(Mb_vec, Mb_vec))
                Mb_angle = max(-1.0, min(1.0, Mb_angle))
                Mb_angle = acos(Mb_angle)
                !Ma_angle = acos(dot_product(Ma_vec, easy_vector(:))/sqrt(dot_product(Ma_vec, Ma_vec)))
                !Mb_angle = acos(dot_product(Mb_vec, - easy_vector(:))/sqrt(dot_product(Mb_vec, Mb_vec)))
                
                canting_angle = Ma_angle/2 + Mb_angle/2
                !write(6,*) 'Easy vector ', easy_vector(:)
                !write(6,*) 'Ma: ', Ma_vec(:)
                !write(6,*) 'Mb: ', Mb_vec(:)
                !write(6,*) 'Ma_angle: ', Ma_angle * 180 / pi
                !write(6,*) 'Mb_angle: ', Mb_angle * 180 / pi
                !write(6,*) 'canting_angle: ', canting_angle * 180 / pi
                !stop
                canting_angle_avg_f = canting_angle_avg_f + Ma_angle
                canting_angle_avg_af = canting_angle_avg_af + Mb_angle
                canting_angle_avg = canting_angle_avg + canting_angle

            end do
            write(12,37) H, canting_angle_avg/neq*180/pi, canting_angle_avg_f/neq*180/pi, canting_angle_avg_af/neq*180/pi
            37 format(F8.3, 2x, F12.8, 2x, F12.8, 2x, F12.8)
            flush(12)
        H = H + dh
        end do 
    end if 


    close (12)

    call cpu_time(end_time)
    time_elapsed = end_time - start_time
    write(13, *) "The program took:", time_elapsed, " seconds to run."
    close(13)


end program hLoop   


