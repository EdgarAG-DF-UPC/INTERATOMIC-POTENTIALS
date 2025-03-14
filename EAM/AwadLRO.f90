REAL*8 FUNCTION phi_Li_Awad(r)
    IMPLICIT NONE
    REAL*8 r, phiatt, phirep
    REAL*8 w, rc, S, sigma, rs
    COMMON/PARAMS/w, rc, S, sigma, rs
    if (r .ge. 2d0 + rs) then
        phi_Li_Awad = phiatt(r)
    else
        phi_Li_Awad = phirep(r)
    endif
    RETURN
END FUNCTION phi_Li_Awad

REAL*8 FUNCTION phiatt(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    !REAL*8 w, rc, S, sigma, rs
    REAL*8, parameter :: A1=0.290231 , A2=1.383589
    REAL*8, parameter :: B1=1.692128, B2=4.237653
    REAL*8, parameter :: C1=1.128630, C2=-1.802154
    REAL*8, parameter :: D1=0.290181, D2=0.290311
    REAL*8, parameter :: E1=1.307873, E2=-0.000074, E3=1.021786 , E4=0.000066
    REAL*8, parameter :: F1=0.740131, F2=0.753327
    REAL*8, parameter :: r1=0.312355, r2=0.787057, r3=1.212659, r4=1.641811, r5=0.524198, r6=1.204905
    REAL*8, parameter :: w=0.825, S=-0.00046, rc=7.5, sigma=1.21, rs=0.07

    phiatt =  w*A1*(r1/(r-rs))**B1*dcos(E1*((r-rs)/r2+C1)) &
            + E2*dexp(D1+F1*(r-rs)/r3) &
            + w*A2*(r4/(r-rs))**B2*dsin(E3*((r-rs)/r5+C2)) &
            + E4*dexp(D2+F2*(r-rs)/r6) &
            + S*dexp((rc-r)/sigma)


    RETURN
END FUNCTION phiatt

REAL*8 FUNCTION phirep(r)
    IMPLICIT NONE
    REAL*8 r
    REAL*8, parameter :: A=1.209005d0, B=2.218267d0, C=2.115713d0, D=2.280588d0, E=8.313021d0, F=-0.757819d0, G=-7.187612d0

    phirep = A*(r/B)**C + D*(r/E)**F + G

    RETURN
    END FUNCTION phirep

    REAL*8 FUNCTION f_Li_Awad(r)
        IMPLICIT NONE
        REAL*8, intent(in) :: r
        REAL*8, parameter :: p1=3.0511d0, p2=1.2200d0
        
        f_Li_Awad = p1 * dexp(-p2*r)
    RETURN
END FUNCTION f_Li_Awad

REAL*8 FUNCTION Femb_Li_Awad(rho)
    IMPLICIT NONE
    REAL*8, intent(in) :: rho
    REAL*8 param(0:6), a(1:7), b(1:7), c(1:7), m
    INTEGER i
    COMMON/LiLiPARAMS/param,a,b,c,m

    m = 0d0 
    param(0) = 1d0 
    param(1) = 0.900d0 
    param(3) = 0.700d0 
    param(2) = 0.840d0 !Belashchenko
    param(4) = 0.550d0 !Belashchenko
    param(5) = 0.350d0 !Belashchenko
    param(6) = 1.100d0 

    a(1) = -1.168d0 !Awad
    a(2) = -1.166700d0 !Awad
    a(3) = -1.159848d0 !Awad
    a(4) = -1.136608d0 !Awad
    a(5) = -1.065193d0 !Awad
    a(6) = -0.863073d0 !Awad
    a(7) = -1.166700d0 !Belashchenko

    b(1) = 0d0 
    b(2) = -0.026000d0 !Awad
    b(3) = -0.202400d0 !Awad
    b(4) = -0.129600d0 !Awad
    b(5) = -0.82600d0 !Awad
    b(6) = -1.198600d0 !Awad
    b(7) = +0.026000d0 !Awad

    c(1) = 0.13d0 !Awad
    c(2) = 1.47d0 !Awad
    c(3) = -0.26d0 !Awad
    c(4) = 2.31d0 !Awad
    c(5) = 0.94d0!Awad
    c(6) = 2.01d0 !Awad
    c(7) = 0.000d0!Awad

    Femb_Li_Awad = 0d0
    if (rho .gt. param(6)) then ! hauria de ser .ge., però llavors dona NaN quan rho=param(6) i m=0 (cas Awad). Com que la funcio i la derivada són contínues, és totalment lícit.
        Femb_Li_Awad = a(7) + b(7)*(rho-param(6)) + c(7)*(rho-param(6))**m
    else if (rho .ge. param(1)) then
        Femb_Li_Awad = a(1) + c(1)*(rho-param(0))**2d0
    else if (rho .lt. param(5)) then
        Femb_Li_Awad =(a(6)+b(6)*(rho-param(5)) + c(6)*(rho-param(5))**2d0)
    else
        do i = 2,5
        if ((rho.ge.param(i)).and.(rho.lt.param(i-1))) then
            Femb_Li_Awad = a(i) + b(i)*(rho-param(i-1)) + c(i)*(rho-param(i-1))**2d0
        endif
        enddo
    endif
    return
END FUNCTION Femb_Li_Awad

REAL*8 FUNCTION phi_pbli_morse(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8, parameter :: eps=0.115d0, r0=3.06d0, alpha=4.35d0 ! eV, Å, adim
	    ! phi_pbli_morse = eps * (dexp(-2d0*alpha*(r-r0)) - 2d0*dexp(-alpha*(r-r0)))     !!!! Notar error a l'apendix de Awad et al (alpha ha de ser adimensional en comptes de 1/Å!!!!!)
    phi_pbli_morse = eps * (dexp(-2d0*alpha*(r/r0-1d0)) - 2d0*dexp(-alpha*(r/r0-1d0)))
	    ! print*, phi_pbli_morse, &
	    ! eps * (dexp(-2d0*alpha*(r/r0-1)) - 2d0*dexp(-alpha*(r/r0-1)))
    return
END FUNCTION phi_pbli_morse