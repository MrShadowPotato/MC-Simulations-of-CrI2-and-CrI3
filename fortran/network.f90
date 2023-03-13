program network_generator
    implicit none

    !Declaracion de variables
    integer :: i, j, k, nx, ny, N
    real(8), dimension(8,3) :: basis_vectors
    real(8), parameter :: a1(3) = [6.8911914825, 0.0, 0.0]
    real(8), parameter :: a2(3) = [-3.4488059462999998, 5.9711638719, 0.0]
    character(len=2), dimension(8) :: element = (/ 'Cr', 'Cr', 'I ', 'I ', 'I ', 'I ', 'I ', 'I ' /)
    
    !Input para las dimensiones de la grilla
    write(6,*)'Numero de celdas eje x?'
    read(5,*)nx
    write(6,*)'Numero de celdas eje y?'
    read(5,*)ny

    !Numero total de atomos
    N = nx*ny*8
    


    !Vectores base
    basis_vectors(1,:) = [3.42245459776438, 1.94691491308500, 9.72065830568899]   !Cr
    basis_vectors(2,:) = [-0.02339545180016, 3.93735766285690, 9.72080135688876]  !Cr
    basis_vectors(3,:) = [-1.03538548904423, 5.92724513968807, 11.28527832915810] !I 
    basis_vectors(4,:) = [2.20273590406674, 3.81437945298318, 11.29027271402490]  !I
    basis_vectors(5,:) = [5.64639472731378, 2.07416105561567, 11.29166984494761]  !I
    basis_vectors(6,:) = [0.97989374494890, 5.92913055695139, 8.15211676904270]   !I
    basis_vectors(7,:) = [1.19604098486642, 2.06898570270974, 8.15546418188776]   !I
    basis_vectors(8,:) = [4.64469384948989, 3.81029176486010, 8.15107346274473]   !I

    
    !Abrir archivo
    open(1, file='CrI3structure.xyz', status='unknown')

    !Numero de atomos y comentario
    write(1,*) N
    write(1,*) 'CrI3 monocapa' 
    !Genera la grilla
    do i=1, nx
        do j=1, ny
            !Vector posicion para el origen de las celdas.
            do k = 1, 8
                write(1, '(A2, 3(F12.6))') element(k), &
                    a1(1)*REAL(i) + a2(1)*REAL(j) + basis_vectors(k,1), &
                    a1(2)*REAL(i) + a2(2)*REAL(j) + basis_vectors(k,2), &
                    a1(3)*REAL(i) + a2(3)*REAL(j) + basis_vectors(k,3)
            end do
        end do
    end do


    close(1)
end program network_generator