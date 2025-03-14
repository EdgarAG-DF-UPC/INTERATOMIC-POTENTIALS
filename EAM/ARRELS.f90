PROGRAM ARRELS
    IMPLICIT NONE
    REAL*8 phi_Zhou, phi_Li_Awad, phi_Pb_2012, phi_Li, r0
    REAL*8 re, fe, rhoe, alpha, beta, A, B, kappa, lambda, Fn(0:3), F0(0:3), eta, FFe
    REAL*8 PSI(0:6), Ae(1:7), Be(1:7), Ce(1:7), me ! EAM Li-Li
    INTEGER m
    REAL*8 aaa(1:3,0:8), rp(1:4) ! EAM Pb-Pb
    CHARACTER*3 aim ! ---------     "    "    -------- 
    EXTERNAL phi_Zhou, phi_Li_Awad, phi_Pb_2012, phi_Li
    COMMON /Zhou/ re, fe, rhoe, alpha, beta, A, B, kappa, lambda, Fn, F0, eta, FFe
    COMMON /BelashchenkoLi/ PSI,Ae,Be,Ce,me
    COMMON /BelashchenkoPb/ aaa,rp

    print*, "### Pb - Zhou ###"
    call PARAMETRES_ZHOU("Pb", "04", re, fe, rhoe, alpha, beta, A, B, kappa, lambda, Fn, F0, eta, FFe)
    call BISECTION(phi_Zhou,3d0,3.5d0,1d-4,r0)
    print*, "SOLUCIO: r0 = ", r0, "Å"

    print*, ""
    print*, ""
    print*, "### Pb - Belashchenko ###"
    rp(1) = 2.60d0 !Angs
    rp(2) = 4.60d0 !Angs
    rp(3) = 7.60d0 !Angs
    rp(4) = 9.01d0 !Angs
    open(97, file = "Aim.Pb.tab")
    do m = 0, 8
        read(97,*) aim, aaa(1,m), aaa(2,m), aaa(3,m)
    enddo
    close(97)
    call BISECTION(phi_Pb_2012,3d0,3.5d0,1d-4,r0)
    print*, "SOLUCIO: r0 = ", r0, "Å"

    print*, ""
    print*, ""
    print*, "### Li - Belashchenko ###"
    me = 1.5d0 
    PSI(0) = 1d0 
    PSI(1) = 0.900d0 
    PSI(2) = 0.840d0 !Belashchenko
    PSI(3) = 0.700d0 
    PSI(4) = 0.550d0 !Belashchenko
    PSI(5) = 0.350d0 !Belashchenko
    PSI(6) = 1.100d0 

    Ae(1) = -0.8948d0 !Belashchenko
    Ae(2) = -0.894474d0 !Belashchenko
    Ae(3) = -0.887963d0 !Belashchenko
    Ae(4) = -0.878482d0 !Belashchenko
    Ae(5) = -0.850369d0 !Belashchenko
    Ae(6) = -0.800385d0 !Belashchenko
    Ae(7) = -0.894474d0 !Belashchenko

    Be(1) = 0d0 
    Be(2) = -0.006520d0 !Belashchenko
    Be(3) = -0.210520d0 !Belashchenko
    Be(4) = 0.07580d0 !Belashchenko
    Be(5) = -0.449920d0 !Belashchenko
    Be(6) = -0.049920d0 !Belashchenko
    Be(7) = 0.006520d0 

    Ce(1) = 0.0326d0 !Belashchenko
    Ce(2) = 1.700d0 !Belashchenko
    Ce(3) =  -1.020d0 !Belashchenko
    Ce(4) =  1.750d0 !Belashchenko
    Ce(5) =  -1.000d0!Belashchenko
    Ce(6) = 11.0d0 !Belashchenko
    Ce(7) =  0.000d0!Belashchenko
    call BISECTION(phi_Li,2.5d0,3d0,1d-4,r0)
    print*, "SOLUCIO: r0 = ", r0, "Å"


    print*, ""
    print*, ""
    print*, "### Li - Awad ###"
    call BISECTION(phi_Li_Awad,2.5d0,3d0,1d-4,r0)
    print*, "SOLUCIO: r0 = ", r0, "Å"

END PROGRAM ARRELS
SUBROUTINE BISECTION(funcio,AA,BB,eps,arrel)
    IMPLICIT NONE
    REAL*8 funcio
    REAL*8 AA, BB, CC, eps, arrel, dif
    INTEGER niter, i, nmax
    REAL*8 fA, fB, fC
    REAL*8 A, B, C
    
    A = AA
    B = BB
    
    fA = funcio(A)
    fB = funcio(B)
    
    nmax = 10000
    
    ! print*, nint(dlog((B-A)/eps)/dlog(2d0)) +1
    
    !nmax = nint(dlog((B-A)/eps)/dlog(2d0)) +1
    
    ! print*, nmax
    
    do i = 1, nmax
       dif = dabs(B-A)
       if (fA*fB .lt. 0d0) then
          C = 0.5d0*(A+B)
          fC = funcio(C)
          if ((fC.eq.0d0).or.(dif.lt.eps)) then
             niter = i
             arrel = C
             exit
          else if (fC*fA .lt. 0d0) then
             B = C
          else if (fC*fB .lt. 0d0) then
             A = C
          endif
       endif
    enddo
    
    print*, "Ha convergit en", niter, "iteracions"
    
    return
END SUBROUTINE BISECTION