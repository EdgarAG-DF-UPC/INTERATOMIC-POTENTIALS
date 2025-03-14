	PROGRAM LiPb_potential
	IMPLICIT NONE
	INTEGER k,m
	REAL*8 r, V, dV, delr
	REAL*8 a,b,c,d,e,f,g,h,i,j
	CHARACTER*4 VAR, LABEL(0:20)
	REAL FACTOR
	REAL*8 rs
	CHARACTER*8 r0s
	EXTERNAL V, dV
	COMMON/PARAM/a,b,c,d,e,f,g,h,i,j
	
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
	
	delr = 1d-2
	
	DO M = 0,20
	
	   FACTOR = 0.05*real(M)
	   
	   
	   open(42,file=LABEL(M)//"_PbHe.table")

	   write(42,'(a)') 
     &	   "# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
	   write(42,'(a)') 
     &	   "# Potencial Pb-He (model Sladek - dU[TQ]Z)"
	   write(42,'(a)')
	   write(42,'(a)') "SLADEK_PbHe"
	   write(42,'(a)') "N 736 R 0.28 15.00"
	   write(42,'(a)')	
!	dU[TQ]Z in Fig.1
!	HF(QZ) + CORR(~r^{-3}, TQ)
	   a = 4.6874d0
	   b = 5.4114d0
	   c = -20.7496d0
	   d = 56.2663d0
	   e = -95.6327d0
	   f = 112.157d0
	   g = -77.8808d0
	   h = 39.4116d0
	   i = 28.2574d0
	   j = -3.7882d0
	   
	   do k=14,750
	      r = float(k)*2d-2
	      write(42,*) k, real(r), 
     &	      FACTOR*real(V(r)*1d-6*27.211386245988), 
     &	      FACTOR*real(-dV(r)*1d-6*27.211386245988/a)
	   enddo
	   
	   close(42)
	
	ENDDO
	call BISECTION(V, 2d0, 7.5d0, 1d-4, rs)
	write(r0s, "(F8.6)") rs
	print*, "V("//r0s//" Å) = 0"
        call BISECTION(dV, 2d0, 7.5d0, 1d-4, rs)
        write(r0s, "(F8.6)") rs
        print*, "V("//r0s//" Å) = Vmin"
	
	delr = 1d-2
	
	
	
!	dU[TQ]Z in Fig.1
!	HF(QZ) + CORR(~r^{-3}, TQ)
	open(14,file="dUTQZ.dat")
	
	open(42,file="PbHe.dUTQZ.table")
	write(42,'(a)') 
     &	"# DATE: 2022-09-07 UNITS: metal CONTRIBUTOR: Edgar"
	write(42,'(a)') 
     &	"# Potencial per a interaccions Pb-He (model Sladek - dU[TQ]Z)"
	write(42,'(a)')
	write(42,'(a)') "SLADEK_PbHe"
	write(42,'(a)') "N 736 R 0.02 14.72"
	write(42,'(a)')
	
	a = 4.6874d0
	b = 5.4114d0
	c = -20.7496d0
	d = 56.2663d0
	e = -95.6327d0
	f = 112.157d0
	g = -77.8808d0
	h = 39.4116d0
	i = 28.2574d0
	j = -3.7882d0
	
	do k=300,1500
	   r = float(k)*delr
	   write(14,*) r, V(r)*1d-6*27.211386245988, 
     &	   dV(r)*1d-6*27.211386245988/a
	enddo
	
	do k=14,750
	   r = float(k)*2d-2
	   write(42,*) k, real(r), real(V(r)*1d-6*27.211386245988), 
     &	   real(-dV(r)*1d-6*27.211386245988/a)
	enddo
	
	close(14)
	close(42)

!	dU[DT]Z in Fig.1
!	HF(QZ) + CORR(~r^{-3}, DT)
	open(15,file="dUDTZ.dat")
	
	open(43,file="PbHe.dUDTZ.table")
	write(43,'(a)') 
     &	"# DATE: 2023-07-03 UNITS: metal CONTRIBUTOR: Edgar"
	write(43,'(a)') 
     &	"# Potencial per a interaccions Pb-He (model Sladek - dU[DT]Z)"
	write(43,'(a)')
	write(43,'(a)') "SLADEK_PbHe"
	write(43,'(a)') "N 736 R 0.02 14.72"
	write(43,'(a)')
	
	a = 4.7861d0
	b = 5.5993d0
	c = -24.411d0
	d = 57.9706d0
	e = 40.4858d0
	f = -281.8289d0
	g = 327.4464d0
	h = 33.4269d0
	i = -140.9091d0
	j = 19.4699d0
	
	do k=300,1600
	   r = float(k)*delr
	   write(15,*) r, V(r)*1d-6*27.211386245988, 
     &	   dV(r)*1d-6*27.211386245988/a
	enddo
	
	do k=1,736
	   r = float(k)*2d-2
	   write(43,*) k, real(r), real(V(r)*1d-6*27.211386245988), 
     &	   real(-dV(r)*1d-6*27.211386245988/a)
	enddo
	
	close(15)
	close(43)


!	dU[DTQ]Z not shown in Fig.1
!	HF(QZ) + CORR(~r^{-3}, DTQ)
	open(16,file="dUDTQZ.dat")
	
	open(44,file="PbHe.dUDTQZ.table")
	write(44,'(a)') 
     &	"# DATE: 2023-07-03 UNITS: metal CONTRIBUTOR: Edgar"
	write(44,'(a)') 
     &	"# Potencial per a interaccions Pb-He (model Sladek - dU[DTQ]Z)"
	write(44,'(a)')
	write(44,'(a)') "SLADEK_PbHe"
	write(44,'(a)') "N 736 R 0.02 14.72"
	write(44,'(a)')
	
	a = 4.7478d0
	b = 5.1242d0
	c = -25.326d0
	d = 63.324d0
	e = -9.3914d0
	f = -178.7298d0
	g = 259.1428d0
	h = 35.8045d0
	i = -136.4474d0
	j = 24.6517d0
	
	do k=300,1600
	   r = float(k)*delr
	   write(16,*) r, V(r)*1d-6*27.211386245988, 
     &	   dV(r)*1d-6*27.211386245988/a
	enddo
	
	do k=1,736
	   r = float(k)*2d-2
	   write(44,*) k, real(r), real(V(r)*1d-6*27.211386245988), 
     &	   real(-dV(r)*1d-6*27.211386245988/a)
	enddo
	
	close(16)
	close(44)
	
	END PROGRAM LiPb_potential
	
	
!	REAL*8 FUNCTION V(r)
!	IMPLICIT NONE
!	REAL*8 x, r
!	REAL*8 a,b,c,d,e,f,g,h,i,j
!	COMMON/PARAM/a,b,c,d,e,f,g,h,i,j
	
!	x = (r-a)/a
	
!	V = -h*exp(-b*x)*(1+b*x+c*(x)**2+
!     &	d*(x)**3+e*(x)**4+f*(x)**5+g*(x)**6+
!     &	i*(x)**7+j*(x)**8)
	
!	return
!	END FUNCTION
	REAL*8 FUNCTION V(r)
	IMPLICIT NONE
	INTEGER k
	REAL*8 x, r
	REAL*8 a(1:8), De, Re
	REAL*8 aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
	COMMON/PARAM/aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
!	COMMON/PARAM/Re,a(1),a(2),a(3),a(4),a(5),a(6),De,a(7),a(8)
	
	Re = aa
	De = hh
	a(1) = bb
	a(2) = cc
	a(3) = dd
	a(4) = ee
	a(5) = ff
	a(6) = gg
	a(7) = ii
	a(8) = jj
	
	x = (r-Re)/Re
	
	V = 1d0
	do k=1,8 
	   V = V + a(k)*x**dble(k)
	enddo
	V = -De*exp(-a(1)*x)*V
	
	return
	END FUNCTION 
 
	
	REAL*8 FUNCTION dV(r)
	IMPLICIT NONE
	INTEGER k
	REAL*8 x, r
	REAL*8 a(1:8), De, Re
	REAL*8 sumk
	REAL*8 aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
	COMMON/PARAM/aa,bb,cc,dd,ee,ff,gg,hh,ii,jj
!	COMMON/PARAM/Re,a(1),a(2),a(3),a(4),a(5),a(6),De,a(7),a(8)
	
	Re = aa
	De = hh
	a(1) = bb
	a(2) = cc
	a(3) = dd
	a(4) = ee
	a(5) = ff
	a(6) = gg
	a(7) = ii
	a(8) = jj
	
	x = (r-Re)/Re
	
	sumk = 0d0
	do k=1,8
	   sumk = sumk + (a(1)*x - k)*a(k)*x**dble(k-1)
	enddo
	dV = De*exp(-a(1)*x)*(a(1) + sumk)
	
	return
	END FUNCTION 

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
