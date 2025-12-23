SUBROUTINE PARAMETRES_ZHOU(elem, versio, reO, feO, rhoeO, rhosO, alphaO, betaO, AO, BO, kappaO, lambdaO, FnO, FO, etaO, FFeO)
    IMPLICIT NONE
    CHARACTER*2, intent(in) :: elem
    ! TYPE PARAMETRE
    !     REAL*8 Cu, Ag, Au, Ni, Pd, Pt, Al, Pb, Fe, Mo, Ta, W, Mg, Co, Ti, Zr
    ! END TYPE PARAMETRE
    ! TYPE(PARAMETRE) re, fe, rhoe, alpha, beta, A, B, kappa, lambda, Fn(0:3), F(0:3), eta, FFe
    CHARACTER*2, dimension(1:16) :: elements
    REAL*8, dimension(1:16) :: re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, eta, FFe
    REAL*8, dimension(1:16,0:3) :: Fn, F
    REAL*8, intent(out) :: reO, feO, rhoeO, rhosO, alphaO, betaO, AO, BO, kappaO, lambdaO, FnO(0:3), FO(0:3), etaO, FFeO
    CHARACTER*6 cnt
    CHARACTER*2 versio
    INTEGER ielement
    

    OPEN(14, FILE="PARAMS"//versio//".tab")
    ! READ(14,*) cnt, element
    ! READ(14,*) cnt, re
    ! READ(14,*) cnt, fe
    ! READ(14,*) cnt, rhoe
    ! READ(14,*) cnt, alpha
    ! READ(14,*) cnt, beta
    ! READ(14,*) cnt, A
    ! READ(14,*) cnt, B
    ! READ(14,*) cnt, kappa
    ! READ(14,*) cnt, lambda
    ! READ(14,*) cnt, Fn(0)
    ! READ(14,*) cnt, Fn(1)
    ! READ(14,*) cnt, Fn(2)
    ! READ(14,*) cnt, Fn(3)
    ! READ(14,*) cnt, F(0)
    ! READ(14,*) cnt, F(1)
    ! READ(14,*) cnt, F(2)
    ! READ(14,*) cnt, F(3)
    ! READ(14,*) cnt, eta
    ! READ(14,*) cnt, FFe
    READ(14,*) cnt, elements
    ielement = findloc(array=elements, value=elem, dim=1)
    READ(14,*) cnt, re
    reO = re(ielement)
    READ(14,*) cnt, fe
    feO = fe(ielement)
    READ(14,*) cnt, rhoe
    rhoeO = rhoe(ielement)
    READ(14,*) cnt, rhos
    rhosO = rhos(ielement)
    READ(14,*) cnt, alpha
    alphaO = alpha(ielement)
    READ(14,*) cnt, beta
    betaO = beta(ielement)
    READ(14,*) cnt, A
    AO = A(ielement)
    READ(14,*) cnt, B
    BO = B(ielement)
    READ(14,*) cnt, kappa
    kappaO = kappa(ielement)
    READ(14,*) cnt, lambda
    lambdaO = lambda(ielement)
    READ(14,*) cnt, Fn(:,0)
    READ(14,*) cnt, Fn(:,1)
    READ(14,*) cnt, Fn(:,2)
    READ(14,*) cnt, Fn(:,3)
    FnO(:) = Fn(ielement, :)
    READ(14,*) cnt, F(:,0)
    READ(14,*) cnt, F(:,1)
    READ(14,*) cnt, F(:,2)
    READ(14,*) cnt, F(:,3)
    FO(:) = F(ielement, :)
    READ(14,*) cnt, eta
    etaO = eta(ielement)
    READ(14,*) cnt, FFe
    FFeO = FFe(ielement)
    CLOSE(14)

END SUBROUTINE PARAMETRES_ZHOU


REAL*8 FUNCTION Femb_Zhou(rho_reduced)
      IMPLICIT NONE
      REAL*8, intent(in) :: rho_reduced
      REAL*8 rho, re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn(0:3), F(0:3), eta, FFe
      INTEGER i
      ! REAL*8, parameter :: rhoe=8.906840d0, rhon=0.85d0*rhoe, &
      !    rho0=1.15d0*rhoe
      REAL*8 rhon, rho0
      ! REAL*8, parameter :: eta=1.172361d0, Fe=-1.440494d0
      ! REAL*8, dimension(0:3), parameter :: Fn=(/-1.419644d0, -0.228622d0, 0.630069d0, -0.560952d0 /), &
      !    F=(/-1.44d0, 0d0, 0.921049d0, 0.108847d0 /)
      COMMON /Zhou/ re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn, F, eta, FFe
   
      rho = rho_reduced*rhoe
      rhon=0.85d0*rhoe
      rho0=1.15d0*rhoe

      Femb_Zhou = 0d0
      if (rho .lt. rhon) then
         do i = 0, 3
            Femb_Zhou = Femb_Zhou + Fn(i)*(rho/rhon-1d0)**dble(i)
         enddo
      elseif (rho .lt. rho0) then
         do i = 0, 3
            Femb_Zhou = Femb_Zhou + F(i)*(rho/rhoe-1d0)**dble(i)
         enddo
      elseif (rho .ge. rho0) then
            Femb_Zhou = FFe * (1d0 - eta*dlog(rho/rhos)) * (rho/rhos)**eta
      else
         print*, "ERROR: rho out of bounds"
         stop
      endif

      ! if (rho0 .le. rho) then
      ! Femb = Fe * (1d0 - eta*dlog(rho/rhoe)) * (rho/rhoe)**eta
      ! elseif (rhon .le. rho) then
      ! do i=0,3
      !    Femb = Femb + F(i)*(rho/rhoe-1d0)**dble(i)
      ! enddo
      ! elseif (rho .lt. rhon) then
      ! do i=0,3
      !    Femb = Femb + Fn(i)*(rho/rhon-1d0)**dble(i)
      ! enddo
      ! else
      ! print*, "ERROR"
      ! stop
      ! endif
END FUNCTION Femb_Zhou

REAL*8 FUNCTION phi_Zhou(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8 re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, eta, FFe
    REAL*8, dimension(0:3) :: Fn, F
    REAL*8 funcio
    COMMON /Zhou/ re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn, F, eta, FFe

    ! phi = A*dexp(-alpha*(r/re-1d0))/(1d0+(r/re-kappa)**20d0) - B*dexp(-beta*(r/re-1d0))/(1d0+(r/re-lambda)**20d0)

    phi_Zhou = (funcio(r,re,A,alpha,kappa) - funcio(r,re,B,beta,lambda))

    return
END FUNCTION phi_Zhou

REAL*8 FUNCTION funcio(r, re, C, gamma, delta)
    IMPLICIT NONE
    REAL*8, intent(in) :: r, re, C, gamma, delta

    funcio = C * dexp(-gamma * (r/re - 1d0)) / (1d0 + (r/re - delta)**20d0)

    return
END FUNCTION funcio