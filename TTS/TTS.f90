PROGRAM TTS
    IMPLICIT NONE
    INTEGER i, j, k, io
	REAL*8 a1(1:3), a2(1:3), a3(1:3)
	REAL*8 alpha(1:3), lambda(1:3), gamma(1:3)
	REAL*8 b(1:3), A(1:3), c(1:3,6:16)
	REAL*8 R_e, D_e, Z_A, Z_B
!	REAL*8 R_e(1:3), D_e(1:3), Z_A(1:3), Z_B(1:3)
	REAL*8 r, Us, Ul,delr,d_Ul,d_Us, r0
	CHARACTER*4 VAR, LABEL(0:20)
	REAL FACTOR
    CHARACTER*2 alk1, alk2
	COMMON/PARAM1/a1, a2, a3, alpha
	COMMON/PARAM2/gamma,lambda,A,b,c
	COMMON/PARAM3/Z_A, Z_B, R_e, D_e

	LABEL(0) = "0.00"
	LABEL(1) = "0.05"
	LABEL(2) = "0.10"
	LABEL(3) = "0.15"
	LABEL(4) = "0.20"
	LABEL(5) = "0.25"
	LABEL(6) = "0.30"
	LABEL(7) = "0.35"
	LABEL(8) = "0.40"
	LABEL(9) = "0.45"
	LABEL(10) = "0.50"
	LABEL(11) = "0.55"
	LABEL(12) = "0.60"
	LABEL(13) = "0.65"
	LABEL(14) = "0.70"
	LABEL(15) = "0.75"
	LABEL(16) = "0.80"
	LABEL(17) = "0.85"
	LABEL(18) = "0.90"
	LABEL(19) = "0.95"
	LABEL(20) =  "1.00"

    delr = 0.02d0

    open(14,file="table2")
        read(14,*)
        read(14,*)
        read(14,*) VAR, a1(1), a1(2), a1(3)
        read(14,*) VAR, a2(1), a2(2), a2(3)
        read(14,*) VAR, a3(1), a3(2), a3(3)
        read(14,*) VAR, alpha(1), alpha(2), alpha(3)
        read(14,*) VAR, lambda(1), lambda(2), lambda(3)
        read(14,*) VAR, gamma(1), gamma(2), gamma(3)
        read(14,*) VAR, b(1), b(2), b(3)
        read(14,*) VAR, A(1), A(2), A(3)
        read(14,*) VAR, C(1,6), C(2,6), C(3,6)
        read(14,*) VAR, C(1,8), C(2,8), C(3,8)
        read(14,*) VAR, C(1,10), C(2,10), C(3,10)
        read(14,*) VAR, C(1,12), C(2,12), C(3,12)
        read(14,*) VAR, C(1,14), C(2,14), C(3,14)
        read(14,*) VAR, C(1,16), C(2,16), C(3,16)
	close(14)

    DO
        READ(*,*,IOSTAT=io) alk1, alk2
        IF (io .gt. 0) THEN
            print*, "... something went wrong ...", io
            continue
        ELSEIF (io .lt. 0) THEN 
            print*, "... end of file ..."
            exit
        ELSE
        DO J = 0,20
            FACTOR = 0.05*real(J)
            call WRITETABLE(alk1, alk2, FACTOR)
            ! call WRITETABLE("Li", FACTOR)
            ! call WRITETABLE("Na", FACTOR)
            ! call WRITETABLE("K#", FACTOR)
            ! call WRITETABLE("Rb", FACTOR)
            ! call WRITETABLE("Cs", FACTOR)
        ENDDO
        ENDIF
    ENDDO

CONTAINS

    SUBROUTINE WRITETABLE(alk1, alk2, fac)
        IMPLICIT NONE
        CHARACTER*2, intent(in) :: alk1, alk2
        REAL, intent(in) :: fac
        CHARACTER*4 LABEL
        REAL*8 UTTS, dUTTS
        REAL*8 UTTS_short, UTTS_long, DUTTS_long, DUTTS_short, rs
        CHARACTER*8 r0s
        INTEGER k
        EXTERNAL UTTS, dUTTS
        COMMON/PARAMS/k

        write(LABEL, "(F4.2)") fac
        ! print*, "Writing the "//alk1//"-"//alk2//" file with factor "//LABEL//"..."

        open(114, file=LABEL//"_"//alk2//alk1//".table")
        IF (alk1 .eq. "He") THEN
            if (alk2 .eq. "Li") then
                k = 2
                call LITHIUM()
            elseif (alk2 .eq. "Na") then
                k = 2
                call SODIUM()
            elseif (alk2 .eq. "K#") then
                k = 2
                call POTASSIUM()
            elseif (alk2 .eq. "Rb") then
                k = 2
                call RUBIDIUM()
            elseif (alk2 .eq. "Cs") then
                k = 2
                call CESIUM()
            elseif (alk2 .eq. "He") then
                k = 3
                call HELIUM()
            endif
        ELSEIF (alk1 .eq. "T#") THEN
            if (alk2.eq. "T#") then
                call TRITIUM("T#", "T#")
                k = 1
            elseif (alk2 .eq. "He") then
                call TRITIUM("T#", "He")
            else
                print*, "... error ..."
                stop
            endif
        ELSE
            print*, "ERROR: variable alk must be either He, T#, Li, Na, K#, Rb or Cs"
            stop
        ENDIF
        do i = 14, nint(15d0/delr)

            r = 0d0 + dble(i)*delr
            r0 = r / 0.5291772d0

            Us = UTTS_short(r0/R_e,a1(k),a2(k),a3(k),alpha(k))
            Ul = UTTS_long(r0/R_e,alpha(k),lambda(k),gamma(k), b(k),A(k),c(k,:))
            d_Us = DUTTS_short(r0/R_e,a1(k),a2(k),a3(k),alpha(k))
            d_Ul = DUTTS_long(r0/R_e,alpha(k),lambda(k),gamma(k),b(k),A(k),c(k,:))

            write(114,*) i, real(r), &
                        FACTOR*real(((Z_A*Z_B/R_e)*Us + D_e*Ul)    *27.2114d0), &
                        FACTOR*real(-((Z_A*Z_B/R_e)*d_Us + D_e*d_Ul)*27.2114d0 &
                        / R_E / 0.5291772d0)
            if (fac.eq.1.) write(42,*) real(r), FACTOR*real(((Z_A*Z_B/R_e)*Us + D_e*Ul)*27.2114d0)

        enddo

        if (label .eq. "1.00") then
            call BISECTION(UTTS, 2d0, 7.5d0, 1d-4, rs)
            write(r0s, "(F8.6)") rs
            print*, "V(r0="//r0s//" Å;  "//alk1//"-"//alk2//") = 0"
            call BISECTION(dUTTS, 5d0, 8d0, 1d-4, rs)
            write(r0s, "(F8.6)") rs
            print*, "V(Rmin="//r0s//" Å;  "//alk1//"-"//alk2//") = V_min"
            print*, "( Re = ", R_e*0.5291772d0, " Å )" 
            print*, ""
        endif
        return
    END SUBROUTINE WRITETABLE

    SUBROUTINE LITHIUM() ! For Li-He (using the H-He parameters, i.e., index 2):
        IMPLICIT NONE
        R_e = 11.47d0
        D_e = 7.36d-6
        Z_A = 2d0
        Z_B = 3d0
        
        ! open(114,file=LABEL(J)//"_LiHe.table")
        write(114,'(a)') "# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
        write(114,'(a)') "# Potencial per a interaccions Li-He (model TTS)"
        write(114,'(a)')
        write(114,'(a)') "SHENGTANG_LiHe"
        write(114,'(a)') "N 736 R 0.28 15.00"
        write(114,'(a)')

        return
    END SUBROUTINE LITHIUM


    SUBROUTINE SODIUM() ! For Na-He (using the H-He parameters, i.e., index 2):
        IMPLICIT NONE
        R_e = 11.85d0
        D_e = 6.96d-6
        Z_A = 2d0
        Z_B = 11d0
        
        ! open(114,file=LABEL(J)//"_NaHe.table")
        write(114,'(a)') "# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
        write(114,'(a)') "# Potencial per a interaccions Na-He (model TTS)"
        write(114,'(a)')
        write(114,'(a)') "SHENGTANG_NaHe"
        write(114,'(a)') "N 736 R 0.28 15.00"
        write(114,'(a)')

        return
    END SUBROUTINE SODIUM


    SUBROUTINE POTASSIUM() ! For K-He (using the H-He parameters, i.e., index 2):
        IMPLICIT NONE
        R_e = 13.50d0
        D_e = 5d-6
        Z_A = 2d0
        Z_B = 19d0
        
        ! open(114,file=LABEL(J)//"_KHe.table")
        write(114,'(a)') "# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
        write(114,'(a)') "# Potencial per a interaccions K-He (model TTS)"
        write(114,'(a)')
        write(114,'(a)') "SHENGTANG_KHe"
        write(114,'(a)') "N 736 R 0.28 15.00"
        write(114,'(a)')

        return
    END SUBROUTINE POTASSIUM


    SUBROUTINE RUBIDIUM() ! For Rb-He (using the H-He parameters, i.e., index 2):
        IMPLICIT NONE
        R_e = 13.86d0
        D_e = 4.76d-6
        Z_A = 2d0
        Z_B = 37d0
        
        ! open(114,file=LABEL(J)//"_RbHe.table")
        write(114,'(a)') "# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
        write(114,'(a)') "# Potencial per a interaccions Rb-He (model TTS)"
        write(114,'(a)')
        write(114,'(a)') "SHENGTANG_RbHe"
        write(114,'(a)') "N 736 R 0.28 15.00"
        write(114,'(a)')

        return
    END SUBROUTINE RUBIDIUM


    SUBROUTINE CESIUM() ! For Cs-He (using the H-He parameters, i.e., index 2):
        IMPLICIT NONE
        R_e = 14.89d0
        D_e = 3.63d-6
        Z_A = 2d0
        Z_B = 55d0
        
        ! open(114,file=LABEL(J)//"_CsHe.table")
        write(114,'(a)') "# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
        write(114,'(a)') "# Potencial per a interaccions Cs-He (model TTS)"
        write(114,'(a)')
        write(114,'(a)') "SHENGTANG_CsHe"
        write(114,'(a)') "N 736 R 0.28 15.00"
        write(114,'(a)')

        return
    END SUBROUTINE CESIUM


    SUBROUTINE HELIUM() ! For He-He (index 3):
        IMPLICIT NONE
        	
        R_e = 5.608d0
        D_e = 3.482d-5
        Z_A = 2d0
        Z_B = 2d0
                
        ! open(145,file=LABEL(J)//"_HeHe.table")
        write(114,'(a)') "# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
        write(114,'(a)') "# Potencial per a interaccions He-He (model TTS)"
        write(114,'(a)')
        write(114,'(a)') "SHENGTANG_HeHe"
        write(114,'(a)') "N 736 R 0.28 15.00"
        write(114,'(a)')

        return
    END SUBROUTINE HELIUM
    SUBROUTINE TRITIUM(e1,e2)
        IMPLICIT NONE
        CHARACTER*2, intent(in) :: e1, e2

        if ((e1 .eq. "T#") .and. (e1 .eq. "T#")) then

            R_e = 7.83d0
            D_e = 2.048d-5
            Z_A = 1d0
            Z_B = 1d0        
            
            ! open(145,file=LABEL(J)//"_HeHe.table")
            write(114,'(a)') "# DATE: 2024-10-01 UNITS: metal CONTRIBUTOR: Edgar"
            write(114,'(a)') "# Potencial per a interaccions T-T (model TTS)"
            write(114,'(a)')
            write(114,'(a)') "SHENGTANG_HH"
            write(114,'(a)') "N 736 R 0.28 15.00"
            write(114,'(a)')
        elseif ((e1 .eq. "T#") .and. (e1 .eq. "He")) then
            R_e = 6.66d0
            D_e = 2.2568d-5
            Z_A = 1d0
            Z_B = 2d0        
            
            ! open(145,file=LABEL(J)//"_HeHe.table")
            write(114,'(a)') "# DATE: 2024-10-01 UNITS: metal CONTRIBUTOR: Edgar"
            write(114,'(a)') "# Potencial per a interaccions T-He (model TTS)"
            write(114,'(a)')
            write(114,'(a)') "SHENGTANG_HHe"
            write(114,'(a)') "N 736 R 0.28 15.00"
            write(114,'(a)')
        else
            print*, "... error ..."
            stop
        endif

        return
    END SUBROUTINE TRITIUM

END PROGRAM TTS

SUBROUTINE BISECTION(funcio,AA,BB,eps,arrel)
    IMPLICIT NONE
    REAL*8 funcio
    REAL*8, intent(in) :: AA, BB, eps
    REAL*8, intent(out) :: arrel
    INTEGER niter, i, nmax
    REAL*8 dif
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

REAL*8 FUNCTION UTTS(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8 UTTS_short, UTTS_long, DUTTS_long, DUTTS_short
    INTEGER k
    REAL*8 r0
    REAL*8 a1(1:3), a2(1:3), a3(1:3)
	REAL*8 alpha(1:3), lambda(1:3), gamma(1:3)
	REAL*8 b(1:3), A(1:3), c(1:3,6:16)
	REAL*8 R_e, D_e, Z_A, Z_B
    COMMON/PARAMS/k
    COMMON/PARAM1/a1, a2, a3, alpha
	COMMON/PARAM2/gamma,lambda,A,b,c
    COMMON/PARAM3/Z_A, Z_B, R_e, D_e
    r0 = r / 0.5291772d0
    UTTS = (Z_A*Z_B/R_e)*UTTS_short(r0/R_e,a1(k),a2(k),a3(k),alpha(k)) &
    + D_e*UTTS_long(r0/R_e,alpha(k),lambda(k),gamma(k), b(k),A(k),c(k,:))
    return
END FUNCTION UTTS
REAL*8 FUNCTION dUTTS(r)
    IMPLICIT NONE
    REAL*8, intent(in) :: r
    REAL*8 UTTS_short, UTTS_long, DUTTS_long, DUTTS_short
    INTEGER k
    REAL*8 r0
    REAL*8 a1(1:3), a2(1:3), a3(1:3)
        REAL*8 alpha(1:3), lambda(1:3), gamma(1:3)
        REAL*8 b(1:3), A(1:3), c(1:3,6:16)
        REAL*8 R_e, D_e, Z_A, Z_B
    COMMON/PARAMS/k
    COMMON/PARAM1/a1, a2, a3, alpha
        COMMON/PARAM2/gamma,lambda,A,b,c
    COMMON/PARAM3/Z_A, Z_B, R_e, D_e
    r0 = r / 0.5291772d0
    dUTTS = (Z_A*Z_B/R_e)*DUTTS_short(r0/R_e,a1(k),a2(k),a3(k),alpha(k)) &
    + D_e*DUTTS_long(r0/R_e,alpha(k),lambda(k),gamma(k), b(k),A(k),c(k,:))
    return
END FUNCTION dUTTS

