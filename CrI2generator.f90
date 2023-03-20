program CrI2generator
    implicit none

    !Declaracion de variables
    integer :: i, j, k, nx, ny, N
    real(8), dimension(6,3) :: basis_vectors
    real(8), parameter :: a1(3) = [3.9350889407837069, 0.0, 0.0]
    real(8), parameter :: a2(3) = [0.0, 7.793970515734439, -0.0155822353114624]
    character(len=2), dimension(6) :: element = (/ 'Cr', 'Cr', 'I ', 'I ', 'I ', 'I '/)
    

    !Input para las dimensiones de la grilla
    write(6,*)'Generador de la estructura cristalina de CrI2'
    write(6,*)'Numero de celdas eje x?'
    read(5,*)nx
    write(6,*)'Numero de celdas eje y?'
    read(5,*)ny


    !Numero total de atomos
    N = nx*ny*6
    

    !Vectores base
    basis_vectors(1,:) = [0.00000000000000, 1.18048232888026, 15.66708540526766]  !Cr
    basis_vectors(2,:) = [1.96754466301406, 5.07735235387682, 15.66009018229279]  !Cr
    basis_vectors(3,:) = [1.96754466301406, 2.14455628147048, 17.30559462841125]  !I 
    basis_vectors(4,:) = [1.96754466301406, 0.21647035484051, 14.02945080633759]  !I
    basis_vectors(5,:) = [0.00000000000000, 6.04092814648801, 17.29708073114964]  !I
    basis_vectors(6,:) = [0.00000000000000, 4.11332427420235, 14.02184874359694]  !I

    
    !Abrir archivo
    open(1, file='./structures/CrI2structure.xyz', status='unknown')


    !Numero de atomos y comentario
    write(1,*) N
    write(1,*) 'monocapa CrI2 ' 


    !Genera la grilla
    do i=1, nx
        do j=1, ny
            do k = 1, 6
                write(1, '(A2, 3(F16.9))') element(k), &
                    a1(1)*REAL(i) + a2(1)*REAL(j) + basis_vectors(k,1), &
                    a1(2)*REAL(i) + a2(2)*REAL(j) + basis_vectors(k,2), &
                    a1(3)*REAL(i) + a2(3)*REAL(j) + basis_vectors(k,3)
            end do
        end do
    end do


    close(1)
end program CrI2generator