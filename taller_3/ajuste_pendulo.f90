program pendulo

    implicit none

    ! Constantes en precisión doble
    double precision, parameter :: PI = dacos(-1.0d0) !Calcula pi con arcocoseno
    character(len=*), parameter :: ARCHIVO_ENTRADA = "pendulo_limpio.dat"
    character(len=*), parameter :: ARCHIVO_SALIDA  = "resultados_ajuste.dat"

    ! Variables para lectura
    integer :: medicion_id, numero_oscilaciones
    double precision :: longitud_cm, angulo_deg, tiempo_medido_s
    integer :: io_status

    ! Variables de acumulación y regresión
    integer :: N
    double precision :: L, T, x, y
    double precision :: sx, sy, sxx, sxy, denominador
    double precision :: a, b, g, R2

    ! Variables para el cálculo de R2
    double precision :: y_prom, y_pred, ss_tot, ss_res

    ! verificación y acumulación de sumatorias

    open(unit=10, file=ARCHIVO_ENTRADA, status='old', action='read', iostat=io_status)

    if (io_status /= 0) then
        print *, "Error: No se pudo abrir el archivo ", ARCHIVO_ENTRADA
        stop
    end if

    N   = 0
    sx  = 0.0d0 ! Asegura la doble precision
    sy  = 0.0d0
    sxx = 0.0d0
    sxy = 0.0d0

    do !Lee cada registro del archivo de entrada, linea por linea
        read(10, *, iostat=io_status) medicion_id, longitud_cm, angulo_deg, numero_oscilaciones, tiempo_medido_s

        if (io_status < 0) exit ! EOF
        if (io_status > 0) then
            print *, "Error al leer registro en ", ARCHIVO_ENTRADA
            close(10)
            stop
        end if

        L = longitud_cm / 100.0d0
        T = tiempo_medido_s / dble(numero_oscilaciones)
        x = L
        y = T**2

        N   = N + 1
        sx  = sx + x
        sy  = sy + y
        sxx = sxx + (x**2)
        sxy = sxy + (x * y)
    end do

    close(10)

    if (N < 2) stop

    ! Cálculo pendiente (a), intercepto (b) Y gravedad (g)

    denominador = (dble(N) * sxx) - (sx**2)

    if (dabs(denominador) < 1.0d-12) stop

    a = ((dble(N) * sxy) - (sx * sy)) / denominador
    b = ((sy * sxx) - (sx * sxy)) / denominador

    g = (4.0d0 * (PI**2)) / a

    ! 3. SEGUNDA LECTURA: Cálculo del R^2

    ! Abre el archivo de entrada nuevamente
    open(unit=10, file=ARCHIVO_ENTRADA, status='old', action='read', iostat=io_status)
    if (io_status /= 0) stop

    y_prom = sy / dble(N)
    ss_tot = 0.0d0
    ss_res = 0.0d0

    do 
        read(10, *, iostat=io_status) medicion_id, longitud_cm, angulo_deg, numero_oscilaciones, tiempo_medido_s
        if (io_status < 0) exit

        L = longitud_cm / 100.0d0
        T = tiempo_medido_s / dble(numero_oscilaciones)
        x = L
        y = T**2

        y_pred = (a * x) + b

        ss_tot = ss_tot + ((y - y_prom)**2)
        ss_res = ss_res + ((y - y_pred)**2)
    end do

    close(10)

    R2 = 1.0d0 - (ss_res / ss_tot)

    print *, "                  Resultados                 "

    print '(A, I5)',        " Número de datos (n)     : ", n
    print '(A, F12.6, A)',  " Pendiente (a)           : ", a, " s^2/m"
    print '(A, F12.6, A)',  " Intercepto (b)          : ", b, " s^2"
    print '(A, F12.6)',     " Coef. Determinación (r2): ", r2
    print '(A, F12.6, A)',  " Aceleración grav. (g)   : ", g, " m/s^2"

    
    ! Salida para la veriicación automática
   
    open(unit=20, file=ARCHIVO_SALIDA, status='replace', action='write', iostat=io_status)
    if (io_status == 0) then
        ! Uso de * para guardar toda la precisión sin recortar decimales
        write(20, *) N, a, b, R2, g
        close(20)
        print *, "Resultados guardados exitosamente en: ", ARCHIVO_SALIDA
    end if

end program pendulo
