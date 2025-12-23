PROGRAM EAMALLOYFILE
    IMPLICIT NONE
    REAL*8 X
    CHARACTER*3 Xc
    INTEGER i, j, k
    CHARACTER*6 modelPb, modelLi
    CHARACTER*7 modelPbLi
    CHARACTER*4 autorLi, autorPb
    CHARACTER*2 verPb, verLi
    CHARACTER*2 dia, mes
    CHARACTER*4 any
    INTEGER, parameter :: nr=9010, nrho=5153
    REAL*8, parameter :: dr=0.001d0, drho=0.001d0, rcutoff=9.01d0
    CHARACTER*10 FORMAT
    REAL*8 r, phi, rho, F!, Femb_Pb, Femb_Li, phi_Pb, phi_Li
    REAL*8, dimension(0:nr) :: rpot_Pb, rpot_Li
    REAL*8 re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn(0:3), F0(0:3), eta, FFe
    REAL*8 PSI(0:6), Ae(1:7), Be(1:7), Ce(1:7), me ! EAM Li-Li
    REAL*8 aaa(1:3,0:8), rp(1:4) ! EAM Pb-Pb
    CHARACTER*3 aim ! ---------     "    "    -------- 
    LOGICAL AWAD, BELASHCHENKO, ZHOU, YUKAWA
    COMMON /Zhou/ re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn, F0, eta, FFe
    COMMON /BelashchenkoLi/ PSI,Ae,Be,Ce,me
    COMMON /BelashchenkoPb/ aaa,rp
    FORMAT = "(E24.16)"

    read(*,*) dia, mes, any
    read(*,*) 
    read(*,*) autorPb, verPb
    modelPb = autorPb//verPb
    read(*,*) 
    read(*,*) autorLi, verLi
    modelLi = autorLi//verLi
    read(*,*)
    read(*,*) modelPbLi
    read(*,*) YUKAWA
    if (YUKAWA) then 
        read(*,*) X
        write(Xc, "(F3.1)") X
    endif

    call PARAMETRES_Pb()
    call PARAMETRES_Li()
    call PARAMETRES_PbLi()

    if (YUKAWA) then
        open(14, file="Pb-"//modelPb//"_Li-"//modelLi//"_"//modelPbLi//"-Yukawa"//Xc//".eam.alloy")
    else
        open(14, file="Pb-"//modelPb//"_Li-"//modelLi//"_"//modelPbLi//".eam.alloy")
    endif
    write(14, *) "# DATE (dd/mm/yyyy): "//dia//"/"//mes//"/"//any
    write(14, *) "# UNITS: metal"
    write(14, *) "# CONTRIBUTORS: E. Alvarez-Galera (mix, files); A.S. Al-Awad (Lithium); X.W. Zhou et al. (Lead);"&
    // "D.K. Belashchenko (Lead, Lithium)"
    write(14, *) "2 Pb Li"
    write(14,*) nrho+1, drho, nr+1, dr, rcutoff

! ATOMIC SECTION

   ! Lead
    write(14,*) "82   207.2   4.902   FCC"
    do k = 0, nrho
       F = Femb_Pb(drho*dble(k))
       write(14,FORMAT) F
    enddo
    do k = 0, nr
       r = dr*dble(k)
       rho = f_Pb(r)
       write(14,FORMAT) rho
       phi = phi_Pb(r)
       if (YUKAWA) call CORRECCIO_YUKAWA(r, phi, 2, 2, X)
       rpot_Pb(k) = r*phi
    enddo
    ! Lithium             
    write(14,*) "3   6.941   3.449   BCC"
    do k = 0, nrho
       F = Femb_Li(drho*dble(k))
       write(14,FORMAT) F
    enddo
    do k = 0, nr
       r = dr*dble(k)
       rho = f_Li(r)
       write(14,FORMAT) rho
       phi = phi_Li(r)
       if (YUKAWA) call CORRECCIO_YUKAWA(r, phi, 1, 1, X)
       rpot_Li(k) = r*phi
    enddo
    print*, "ATOMIC SECTION finished"
! POTENTIAL SECTION
    write(14,FORMAT) 0d0 ! ---> Infinit (segons quina versió de LAMMPS es queixa quan troba un NaN al fitxer) 
    do k=1,nr
        write(14,FORMAT) rpot_Pb(k)
    enddo
    write(14,FORMAT) 0d0 ! ---> Infinit (segons quina versió de LAMMPS es queixa quan troba un NaN al fitxer) 
    do k=1,nr
        r = dr*dble(k)
        phi = phi_PbLi(r)
        if (YUKAWA) call CORRECCIO_YUKAWA(r, phi, 1, 2, X)
        write(14,FORMAT) r*phi
    enddo
    write(14,FORMAT) 0d0 ! ---> Infinit (segons quina versió de LAMMPS es queixa quan troba un NaN al fitxer) 
    do k=1,nr
        write(14,FORMAT) rpot_Li(k)
    enddo
    print*, "POTENTIAL SECTION finished"
    
    close(14)

    print*, "eam.alloy file is ready to use :)"
CONTAINS
REAL*8 FUNCTION phi_Pb(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8 phi_Zhou, phi_pb_2012

    if ((modelPb .eq. "Zhou01") .or. (modelPb .eq. "Zhou04")) then
        phi_Pb = phi_Zhou(r)
    elseif (modelPb .eq. "Bela12") then
        phi_Pb = phi_pb_2012(r)
    else
        print*, "ERROR: no phi_Pb selected"
        stop
    endif
    return
END FUNCTION phi_Pb
REAL*8 FUNCTION Femb_Pb(rho)
    IMPLICIT NONE
    REAL*8, intent(in) :: rho
    REAL*8 Femb_Zhou, Femb_pb_2012, dFemb_pb_2012, dFemb_Pb

    if ((modelPb .eq. "Zhou01") .or. (modelPb .eq. "Zhou04")) then
        Femb_Pb = Femb_Zhou(rho)
    elseif (modelPb .eq. "Bela12") then
        call embedding_Pb(rho, Femb_Pb, dFemb_Pb)
    else
        print*, "ERROR: no Femb_Pb selected"
        stop
    endif
    return
END FUNCTION Femb_Pb
REAL*8 FUNCTION f_Pb(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8 funcio, df_Pb

    if ((modelPb .eq. "Zhou01") .or. (modelPb .eq. "Zhou04")) then
        f_Pb = funcio(r, re, Fe, beta, lambda) / rhoe
    elseif (modelPb .eq. "Bela12") then
        call electro(r, f_Pb, df_Pb, 5.1531d0, 1.2200d0)
    else
        print*, "ERROR: no f_Pb selected"
        stop
    endif
    return
END FUNCTION f_Pb
REAL*8 FUNCTION phi_Li(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8 phi_Li_Awad, phi_li_2009, phi_li_2011

    if (modelLi .eq. "Awad23") then
        phi_Li = phi_Li_Awad(r)
    elseif ((modelLi .eq. "Bela09")&
            .or. &
            (modelLi .eq. "Bela14")) then
        phi_Li = phi_li_2009(r)
    elseif (modelLi .eq. "Bela11") then
        phi_Li = phi_li_2011(r)
    else
        print*, "ERROR: no phi_Li selected"
        stop
    endif
    return
END FUNCTION phi_Li
REAL*8 FUNCTION Femb_Li(rho)
    IMPLICIT NONE
    REAL*8, intent(in) :: rho
    REAL*8 Femb_Li_Awad, dFemb_Li

    if (modelLi .eq. "Awad23") then
        Femb_Li = Femb_Li_Awad(rho)
    elseif ((modelLi .eq. "Bela09") &
            .or. &
            (modelLi .eq. "Bela14") &
            .or. &
            (modelLi .eq. "Bela11")) then
        call embedding_Li(rho, Femb_Li, dFemb_Li)
    else
        print*, "ERROR: no Femb_Li selected"
        stop
    endif
    return
END FUNCTION Femb_Li
REAL*8 FUNCTION f_Li(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8 f_Li_Awad, df_Li
    REAL*8 p1, p2

    if (modelLi .eq. "Awad23") then
        p1 = 3.0511d0
        p2 = 1.2200d0
    elseif ((modelLi .eq. "Bela09") .or. (modelLi .eq. "Bela14")) then
        p1 = 3.0450d0
        p2 = 1.2200d0
    elseif (modelLi .eq. "Bela11") then
        p1 = 3.0511d0
        p2 = 1.2200d0
    else
        print*, "ERROR: no f_Li selected"
        stop
    endif
    call electro(r, f_Li, df_Li, p1, p2)
    return
END FUNCTION f_Li
REAL*8 FUNCTION phi_PbLi(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8 phi_pbli_morse, pot_84
    REAL*8, parameter :: eps = 0.134d0, r0 = 3.102d0, alpha = 4.58d0 !eV, Å, -
    REAL*8, parameter :: eps_84 = 0.0880d0, sig_84 = 2.553d0 !eV, Å

    if (modelPbLi .eq. "Morse23") then
        phi_PbLi = phi_pbli_morse(r)
    elseif (modelPbLi .eq. "Morse19") then
        phi_PbLi = eps * (dexp(-2d0*alpha*(r/r0-1d0)) - 2d0*dexp(-alpha*(r/r0-1d0))) 
    elseif (modelPbLi .eq. "LJ84_19") then
        phi_PbLi = eps_84*pot_84(r/sig_84)
    ! elseif (modelPbLi .eq. "Fraile") then
    !     phi_PbLi = ...
    else
        print*, "ERROR: no phi_PbLi selected"
        stop
    endif
END FUNCTION phi_PbLi
SUBROUTINE PARAMETRES_Li()
    IMPLICIT NONE
    integer i
    if (modelLi .eq. "Bela09") then
        print*, "Li -----> Belashchenko EAM (2009)"
        me = 1.5d0 
        PSI(0) = 1d0 
        PSI(1) = 0.900d0 
        PSI(2) = 0.840d0 
        PSI(3) = 0.700d0 
        PSI(4) = 0.480d0 
        PSI(5) = 0.420d0 
        PSI(6) = 1.100d0 

        Ae(1) = -0.880900d0 
        ! Ae(2) = -0.880470d0 
        ! Ae(3) = -0.872754d0 
        ! Ae(4) = -0.859314d0 
        ! Ae(5) = -0.808846d0 
        ! Ae(6) = -0.776842d0 
        ! Ae(7) = -0.881043d0 

        Be(1) =  0d0 
        ! Be(2) = -0.008600d0 
        ! Be(3) = -0.248600d0 
        ! Be(4) =  0.056600d0 
        ! Be(5) = -0.808846d0 
        ! Be(6) = -0.551400d0 
        ! Be(7) =  0.006520d0
        Be(7) = 0d0 

        Ce(1) =  0.0430d0!0.043000d0 
        Ce(2) =  2.000d0!2.000000d0 
        Ce(3) = -1.090d0!-1.090000d0 
        Ce(4) =  1.300d0!1.300000d0 
        Ce(5) =  0.300d0!-0.515400d0 
        ! Ce(6) =  0.000000d0 
        ! Ce(7) =  0d0!0.018130d0 

        ! Ae(2) = Ae(1) + Be(1)*(PSI(1)-PSI(0)) + Ce(1)*(PSI(1)-PSI(0))**2d0
        ! Be(2) = 2d0*Ce(1)*(PSI(1)-PSI(0))
        
        ! Ae(3) = Ae(2) + Be(2)*(PSI(2)-PSI(0)) + Ce(2)*(PSI(2)-PSI(0))**2d0
        ! Be(3) = 2d0*Ce(2)*(PSI(2)-PSI(0))
        do i = 2,5
            Ae(i) = Ae(i-1) + Be(i-1)*(PSI(i-1)-PSI(i-2)) &
                    + Ce(i-1)*(PSI(i-1)-PSI(i-2))**2d0
            Be(i) = Be(i-1) + Ce(i-1)*(PSI(i-1)-PSI(i-2))
        enddo
        Ce(6) = 0d0!(Ae(6) + Be(6)*(PSI(5)-PSI(1)) + Ce(5)*PSI(5)**2d0) &
                !/ (2d0*PSI(1)*PSI(5) - PSI(1)**2d0)
        Ae(6) = Ae(5) + Be(5)*(PSI(5)-PSI(4)) - Be(6)*(PSI(5)-PSI(1)) &
                + Ce(5)*(PSI(5)-PSI(4))**2d0 - Ce(6)*(PSI(5)-PSI(1))**2d0
        Be(6) = Be(5) + 2d0*Ce(5)*(PSI(5)-PSI(4)) - 2d0*Ce(6)*(PSI(5)-PSI(1))
        

        Ae(7) = Ae(1) - Ce(1)/3d0*(PSI(6)-PSI(0))**2d0
        Ce(7) = 4d0/3d0 * Ce(1) * (PSI(6)-PSI(0))**0.5d0
    elseif (modelLi .eq. "Bela14") then
        print*, "Li -----> Belashchenko EAM (2009) - Modified by Fraile (2014)"
        me = 1.5d0 
        PSI(0) = 1d0 
        PSI(1) = 0.900d0 
        PSI(2) = 0.849d0 
        PSI(3) = 0.700d0 
        PSI(4) = 0.480d0 
        PSI(5) = 0.420d0 
        PSI(6) = 1.100d0 

        Ae(1) = -0.880900d0 
        Ae(2) = -0.880470d0 
        Ae(3) = -0.872754d0 
        Ae(4) = -0.859314d0 
        Ae(5) = -0.808846d0 
        Ae(6) = -0.776842d0 
        Ae(7) = -0.881043d0 

        Be(1) =  0d0 
        Be(2) = -0.008600d0 
        Be(3) = -0.248600d0 
        Be(4) =  0.056600d0 
        Be(5) = -0.808846d0 
        Be(6) = -0.551400d0 
        Be(7) =  0.006520d0 

        Ce(1) =  0.043000d0 
        Ce(2) =  2.000000d0 
        Ce(3) = -1.090000d0 
        Ce(4) =  1.300000d0 
        Ce(5) =  -0.515400d0 
        Ce(6) =  0.000000d0 
        Ce(7) =  0.018130d0 

    elseif (modelLi .eq. "Bela11") then
        print*, "Li -----> Belashchenko EAM (2011-2012)"
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
    elseif (modelLi .eq. "Awad23") then
        print*, "Li -----> Awad EAM"
        me = 0d0 
        PSI(0) = 1d0 
        PSI(1) = 0.900d0 
        PSI(2) = 0.840d0 !Belashchenko
        PSI(3) = 0.700d0 
        PSI(4) = 0.550d0 !Belashchenko
        PSI(5) = 0.350d0 !Belashchenko
        PSI(6) = 1.100d0 

        Ae(1) = -1.168d0 !Awad
        Ae(2) = -1.166700d0 !Awad
        Ae(3) = -1.159848d0 !Awad
        Ae(4) = -1.136608d0 !Awad
        Ae(5) = -1.065193d0 !Awad
        Ae(6) = -0.863073d0 !Awad
        Ae(7) = -1.166700d0 !Belashchenko

        Be(1) = 0d0 
        Be(2) = -0.026000d0 !Awad
        Be(3) = -0.202400d0 !Awad
        Be(4) = -0.129600d0 !Awad
        Be(5) = -0.82600d0 !Awad
        Be(6) = -1.198600d0 !Awad
        Be(7) = +0.026000d0 !Awad

        Ce(1) = 0.13d0 !Awad
        Ce(2) = 1.47d0 !Awad
        Ce(3) = -0.26d0 !Awad
        Ce(4) = 2.31d0 !Awad
        Ce(5) = 0.94d0!Awad
        Ce(6) = 2.01d0 !Awad
        Ce(7) = 0.000d0!Awad
    endif
END SUBROUTINE PARAMETRES_Li
SUBROUTINE PARAMETRES_Pb()
    IMPLICIT NONE
    INTEGER m
    if (modelPb .eq. "Bela12") then ! Belashchenko EAM Pb
        print*, "Pb -----> Belashchenko EAM (2012)"
        rp(1) = 2.60d0 !Angs
        rp(2) = 4.60d0 !Angs
        rp(3) = 7.60d0 !Angs
        rp(4) = 9.01d0 !Angs
        open(97, file = "Aim.Pb.tab")
        do m = 0, 8
            read(97,*) aim, aaa(1,m), aaa(2,m), aaa(3,m)
        enddo
        close(97)
    elseif ((modelPb .eq. "Zhou01") .or. (modelPb .eq. "Zhou04")) then
        print*, "Pb -----> Zhou EAM (20"//verPb//")"
        call PARAMETRES_ZHOU("Pb", verPb, re, fe, rhoe, rhos, alpha, beta, A, B, kappa, lambda, Fn, F0, eta, FFe)
    else
        print*, "ERROR: Pb must be:"
        print*, "---> Bela12"
        print*, "---> Zhou01"
        print*, "---> Zhou04"
        print*, "instead of ", modelPb
        stop
    endif
END SUBROUTINE PARAMETRES_Pb
SUBROUTINE PARAMETRES_PbLi()
    IMPLICIT NONE
    if (modelPbLi .eq. "Morse23") then
        print*, "Pb-Li --> Morse-Awad (2023)"
    ! elseif (modelPbLi .eq. "Morse19") then
    !     phi_PbLi = ...
    elseif (modelPbLi .eq. "LJ84_19") then
        print*, "Pb-Li --> Belashchenko LJ 8-4 (2023)"
    elseif (modelPbLi .eq. "Morse19") then
        print*, "Pb-Li --> Morse-Belashchenko (2019)"
    else
        print*, "ERROR: no phi_PbLi selected. Select:"
        print*, "---> Morse23"
        print*, "---> LJ84_19"
        print*, "instead of ", modelPbLi
        stop
    endif
END SUBROUTINE PARAMETRES_PbLi
END PROGRAM EAMALLOYFILE
