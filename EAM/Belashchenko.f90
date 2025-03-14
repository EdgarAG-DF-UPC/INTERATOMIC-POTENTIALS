REAL*8 FUNCTION phi_li(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8, dimension(1:6), parameter :: k = (/ -0.161539351212d1, 0.329193195820d2, -0.245830404172d3,&
    0.840217873656d3, -0.136938125679d4, 0.905623694715d3 /)
    REAL*8, parameter :: k1=0.252868d0, k2=0.15252d0, k3=0.38d0, k4=1.96d0
    INTEGER i

    if (r .le. 2.45d0) then
        phi_li = k1 + k2*(2.45d0-r) + k3*(dexp(k4*(2.45d0-r))-1d0)
    elseif (r .gt. 2.45d0) then
        phi_li = 0d0
        do i = 1, 6
            phi_li = phi_li + k(i) * r**(1-i)
        enddo
    endif

    return
END FUNCTION phi_li

REAL*8 FUNCTION phi_li_2009(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8, dimension(1:6), parameter :: k = (/ -0.161539351212d1, 0.329193195820d2, -0.245830404172d3,&
    0.840217873656d3, -0.136938125679d4, 0.905623694715d3 /)
    REAL*8, parameter :: k1=0.252868d0, k2=0.15252d0, k3=0.38d0, k4=1.96d0
    INTEGER i


    if (r .le. 2.25d0) then
        phi_li_2009 = 0.401569d0*dexp(2.100d0*(2.35d0-r))
    elseif (r .gt. 2.25d0) then
        phi_li_2009 = 0d0
        do i = 1, 6
            phi_li_2009 = phi_li_2009 + k(i) * r**(1-i)
        enddo
    endif

    return
END FUNCTION phi_li_2009

REAL*8 FUNCTION phi_li_2011(r)
    ! returns the pairwise potential of lithium-lithium (in eV)
    IMPLICIT NONE
    REAL*8 r, phi
    REAL*8 k(0:5), k1, k2, k3, k4
    INTEGER i
    k(0) = -0.161539351212d1
    k(1) = 0.329193195820d2
    k(2) = -0.245830404172d3
    k(3) = 0.840217873656d3
    k(4) = -0.136938125679d4
    k(5) = 0.905623694715d3
    k1 = 0.252868d0
    k2 = 0.15252d0
    k3 = 0.38d0
    k4 = 1.96d0
    phi = 0d0
    if (r .le. 7.5d0) then
        if (r .gt. 2.45d0) then
            do i=0,5
                phi = phi + k(i)/r**dble(i)
            enddo
        else
            phi = k1+k2*(2.45d0-r)+k3*(dexp(k4*(2.45d0-r))-1d0)
        endif
    endif
    phi_li_2011 = phi
    return
END FUNCTION phi_li_2011
REAL*8 FUNCTION phi_pbli_2019(r)
    IMPLICIT NONE
    REAL*8 r
    phi_pbli_2019 = 4d0 * (r**(-8d0) - r**(-4d0))
    return
END FUNCTION phi_pbli_2019

REAL*8 FUNCTION phi_pb_2012(r)
    IMPLICIT NONE
    INTEGER i, n, m, k
    PARAMETER(n=3,k=8)
    REAL*8 r, phi, a(1:n,0:k), rp(1:4)
    COMMON/BelashchenkoPb/a,rp
    ! print*, rp
    ! do m = 0, k
    !     print*, a(:,m)
    ! enddo
    ! stop
    phi = 0d0
    if ( r .gt. 2.60d0 ) then
       do i = 1, n
          do m = 0, k
          if ((r .gt. rp(i)).and.(r .le. rp(i+1))) then
             phi = phi + a(i,m)*(r-rp(i+1))**m 
          endif
          enddo
       enddo
    else
       phi = 0.438472d0 - 3.99326d0*(2.60d0-r) + 2.8d0 *(dexp(1.96d0*(2.60d0-r)) - 1d0)
    endif
    phi_pb_2012 = phi
    return
END FUNCTION phi_pb_2012

SUBROUTINE electro(r,psi,dpsidr,p1,p2)
! returns the electronic density of the LM
    IMPLICIT NONE
    REAL*8, intent(in) :: r, p1, p2
    REAL*8, intent(out) ::  psi, dpsidr 
    
    psi = p1 * dexp(-p2*r)
    dpsidr = -p2 * psi

    return
END SUBROUTINE electro
    
SUBROUTINE embedding_Li(rho,F,dFdrho)
    IMPLICIT NONE
    REAL*8, intent(in) :: rho
    REAL*8, intent(out) :: F, dFdrho
    REAL*8 param(0:6), a(1:7), b(1:7), c(1:7), m
    INTEGER i
    COMMON/BelashchenkoLi/param,a,b,c,m
    
    F = 0d0
    if (rho .gt. param(6)) then ! hauria de ser .ge., però llavors dona NaN quan rho=param(6) i m=0 (cas Awad). Com que la funcio i la derivada són contínues, és totalment lícit.
        F = a(7) + b(7)*(rho-param(6)) + c(7)*(rho-param(6))**m
        dFdrho = b(7) + m*c(7)*(rho-param(6))**(m-1d0)
    else if (rho .ge. param(1)) then
        F = a(1) + c(1)*(rho-param(0))**2d0
        dFdrho = 2d0*c(1)*(rho-param(0))
    else if (rho .lt. param(5)) then
        F =(a(6)+b(6)*(rho-param(5)) + c(6)*(rho-param(5))**2d0)
        dFdrho = b(6) + 2d0*c(6)*(rho-param(5))
    else
        do i = 2,5
        if ((rho.ge.param(i)).and.(rho.lt.param(i-1))) then
            F = a(i) + b(i)*(rho-param(i-1)) + c(i)*(rho-param(i-1))**2d0
            dFdrho = b(i) + 2d0*c(i)*(rho-param(i-1))
        endif
        enddo
    endif
    return
END SUBROUTINE embedding_Li
    
    
SUBROUTINE embedding_Pb(rho,F,dFdrho)
    IMPLICIT NONE
    REAL*8, intent(in) :: rho
    REAL*8, intent(out) :: F, dFdrho
    REAL*8 param(0:6), a(1:6), b(1:6), c(1:6), m
    REAL*8 dif, frac
    INTEGER i
    
    m = 1.60d0 
    param(0) = 1d0 
    param(1) = 0.90d0 
    param(2) = 0.81d0 
    param(3) = 0.77d0 
    param(4) = 0.71d0 
    param(5) = 0.46d0 
    param(6) = 1.40d0

    a(1) = -1.5186d0 
    a(2) = -1.500978d0 
    a(3) = -1.469082d0 
    a(4) = -1.4485844d0 
    a(5) = -1.425158d0 
    a(6) = -1.297423d0 

    b(1) = 0d0 
    b(2) = -0.35244d0 
    b(3) = -0.35244d0 
    b(4) = -0.67244d0 
    b(5) = -0.10844d0 
    b(6) = -0.91344d0 

    c(1) = 1.7622d0 
    c(2) = 0.0000d0 
    c(3) = 4.00d0 
    c(4) = -4.70d0 
    c(5) = 1.61d0
    c(6) = -5.70d0
    
    F = 0d0
                
    if (rho .le. param(5)) then
        dif = rho - param(5)
        frac = rho / param(5)
        F =(a(6)+b(6)*dif+c(6)*dif**2d0)*(2d0*frac-frac**2d0)
        dFdrho = (b(6) + 2d0*c(6)*dif)*(2d0*frac-frac**2d0) &
        + (a(6)+b(6)*dif+c(6)*dif**2d0)*2d0/param(5)*(1d0-frac)
    else if (rho .ge. param(1)) then
        dif = rho - param(0)
        F = a(1)+c(1)*dif**2d0
        dFdrho = 2d0*c(1)*dif
    else
        do i = 2,5
        if ((rho.le.param(i-1))) then
            dif = rho - param(i-1)
            F = a(i) + b(i)*dif + c(i)*dif**2d0
            dFdrho = b(i) + 2d0*c(i)*dif
        endif
        enddo
    endif
    
    return
END SUBROUTINE embedding_Pb

REAL*8 FUNCTION pot_84(r)
    IMPLICIT NONE
    REAL*8 r
    pot_84 = 4d0 * (r**(-8d0) - r**(-4d0))
    return
END FUNCTION pot_84

SUBROUTINE CORRECCIO_YUKAWA(r, pot, i, j, X1)
    IMPLICIT NONE
    INTEGER, intent(in) :: i, j
    REAL*8, intent(in) :: r, X1
    REAL*8, intent(inout) :: pot
    REAL*8, parameter :: lambda=1.10d0
    REAL*8, dimension(1:2) :: q
    REAL*8 aij

    if (X1 .eq. 0.9d0) then
        q(1) = 0.11d0
        q(2) = -0.99d0
    elseif (X1 .eq. 0.8d0) then
        q(1) = 0.13d0
        q(2) = -0.52d0
    elseif (X1 .eq. 0.7d0) then
        q(1) = 0.141d0
        q(2) = -0.329d0
    elseif (X1 .eq. 0.6d0) then
        q(1) = 0.0926d0
        q(2) = -0.139d0
    elseif (X1 .eq. 0.5d0) then
        q(1) = 0.075d0
        q(2) = -0.075d0
    else
        q(1) = 0d0
        q(2) = 0d0
    endif

    aij = q(i) * q(j) * 14.3996d0 !14.3996 factor conversió (computational units --> eV)

    pot = pot + (aij/r) * dexp(-lambda*r)

    return
END SUBROUTINE CORRECCIO_YUKAWA