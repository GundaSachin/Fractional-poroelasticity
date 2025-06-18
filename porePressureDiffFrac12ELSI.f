cccRoutine to solve the fractional porepressure diffusion equation to use in a couple thermomechanic problem
ccc    OLGA BARRERA 6/07/2018
C!======================================================================
      module global
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      integer elCount,err
      real*8 deltap(1000,28110,8,3)
      real*8, allocatable :: globalSdv(:,:,:)
      end module global
C!======================================================================
      subroutine UVARM(uvar,direct,t,time,dtime,cmname,orname,
     1 nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord,
     2 jmac,jmatyp,matlayo,laccfla)
      ! This subroutine is used to transfer SDV's from the UEL
      ! onto the dummy mesh for viewing.
      use global
      include 'ABA_PARAM.INC'
      character*80 cmname,orname
      character*3 flgray(15)
      dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
      dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)

       real*8 sigmaH0
!       COMMON hydr(100,1,1,1)  ! hydrostatic stress(nincr, noel, gp, componenti (1))
!
!ccc GET INCREMENT OF sigmaH (at the beginning of the current increment)ccccc
      !Write(*,*)"uvarm called"
      call getvrm('S',array,jarray,flgray,jrcd,
     $     jmac, jmatyp, matlayo, laccfla)
!
       sigmaH0 = globalSdv(noel,npt,1)
      sigmaH = (1.d0/3.d0)*(array(1)+array(2)+array(3))
       !if (noel == 10 .and. npt == 5 ) then
       !  write(*,*)"sigmaH0", sigmaH0
       !  write(*,*)"sigmaH", sigmaH
       !end if
      globalSdv(noel,npt,1)=sigmaH
!!!      hydr(kinc,noel,npt,1)=sigmaH
        dsigmaH = sigmaH-sigmaH0
! !     (dsigmaH= STATEV(1)-hydr(kinc-1,noel,npt,1))
      globalSdv(noel,npt,2)= dsigmaH
!
       ! WRITE(6,*) 'sigmaH', sigmaH
       !WRITE(6,*) 'dsigmaH', dsigmaH
       !WRITE(6,*) 'globalsdv', globalSdv(noel,npt,1)
      return
      end subroutine uvarm
C=======================================================================
      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      use global
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3)
       real*8 COND,DrdSigmaHdt,beta,lambdabeta,kapbalpha
        !COMMON deltap(1000,28110,8,3)  !dimen of this depends on the model c3d8t (4 gauss point)
        real*8 bina,eps  !
      dimension bina(kinc+1),eps(kinc,ntgrd) ! need to see this
C     1 termk(3)
      integer kinc2
      !write(*,*)"started"
C
c      OPEN(14,FILE='C:\temp\porepress.out')
      beta = PROPS(1)  ! K
      lambdabeta=PROPS(2)
      kapbalpha=PROPS(3)
      B = PROPS(4)
      numelem= 28110
      nintp = 8
      if (noel .eq. 1 .and. npt .eq. 1) then
        write(*,*)"kinc",kinc
      end if
      !WRITE(14,*) 'beta ', beta
      ! WRITE(14,*) 'lambdabeta', lambdabeta
      !WRITE(14,*) 'kapbalpha', kapbalpha
C      SKEMP = PROPS(2) ! B
C      PERM= PROPS(3) ! k
C      BIOT= PROPS(4) !alfa
C      VISC= PROPS(5) !fluid viscosity
C      COND=(BULK*SKEMP*PERM)/(BIOT*VISC) ! FLUX=-COND*DP
c      WRITE(14,*) 'BULK', BULK
c     WRITE(14,*) 'SKEMP', SKEMP
c      WRITE(14,*) 'PERM', PERM
c      WRITE(14,*) 'BIOT', BIOT
      !WRITE(14,*) 'kinc', KINC
      ! WRITE(14,*) 'NOEL', NOEL
      !WRITE(14,*) 'GAUSS POINT', NPT
      !WRITE(14,*) 'ntgrd', ntgrd
      !WRITE(14,*) 'DTEMDX', DTEMDX
      !WRITE(14,*) 'DTEMp', DTEMP
      !WRITE(14,*)'dtime',DTIME
CCC ANALOGY with the thermal problem
Ccc   U is  dP/dT= U+(du/dp)dp+ B*dsigmaH/dt
CCCC TEMP=P
CCC  DTEMP=DP
CCC DUDT=DUDP
CC NEED GETVRM TO extrapolate dsigmaH

      if(.not.allocated(globalSdv)) then
         allocate(globalSdv(numElem,nIntp,NSTATV))
      endif
      do i=1,NSTATV
        STATEV(i) = globalSdv(noel,npt,i)
      end do
c      WRITE(14,*) 'kin PRIMA DEL deltap', kinc
c        WRITE(14,*) 'kinc2', kinc2
      deltap(kinc,noel,npt,1:NTGRD)=DTEMDX(1:NTGRD)
      eps(1:kinc,1:ntgrd)=deltap(1:kinc,noel,npt,1:ntgrd)

c       WRITE(14,*) 'deltap', deltap
c       WRITE(14,*) 'eps', eps
C EVALUATE BINOMIAL COEFFICIENTS
C
      bina(1)=1
      do k=1,kinc
          bina(k+1)=(k-1-beta)*bina(k)/(k)
      end do
c      write(14,*) 'bina',bina
      !write(*,*) "state", STATEV(2)
      DUDT = 1.0d0
      DU = DUDT*DTEMP + B*STATEV(2)   !!!-B*sigmaH???
      !write(*,*)"statev",statev(2)
      U = U+DU
      !write(*,*)"DU",DU,"U",U
!      WRITE(14,*) 'STATEV(2)', STATEV(2)

C FLUX TERMS: FLUX=-COND*GRAD(P)- [LAMBDAP*D0(DTEMDX)
CCC GRAD(P)=DTEMDX
c      write(14,*)'FLUXPRIMA',FLUX

      DO I=1, NTGRD
          FLUX(I)=0.D0
          do k1=1,kinc
C    ! FLUX(I)=FLUX(I)-(kapbalpha*lambdabeta)*eps(kinc-k1+1,I)*bina(k1)/
C   !1 (dtime**beta)
      FLUX(I)=FLUX(I)-(kapbalpha*lambdabeta)*eps(kinc-k1+1,I)*bina(k1)/
     1 (dtime**beta)

         DFDG(I,I)= -(kapbalpha*lambdabeta)/(dtime**beta)!-COND
         STATEV(I)=FLUX(I)
         end do
      END DO
!C
!      WRITE(14,*) 'KINC', KINC
!        WRITE(14,*) 'KINC2', KINC2
        !WRITE(14,*)'EPS',eps
        !write(14,*)'flux', FLUX
        !WRITE(14,*)'dfdf', DFDG


      ! WRITE(14,*) 'U', U
      !  WRITE(14,*) 'DU', DU
      !   WRITE(14,*) 'DTEMP', DTEMP
      !   WRITE(14,*) 'TEMP', TEMP
      !  WRITE(14,*) 'DUDT', DUDT
      ! WRITE(14,*) 'NTGRD', NTGRD
      ! WRITE(14,*) 'FLUX', FLUX
      ! WRITE(14,*) 'DTEMDX', DTEMDX
      RETURN
      END

       SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
!C
      use global
      INCLUDE 'ABA_PARAM.INC'
!C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
       real*8 sigmaH0
!       COMMON hydr(100,1,1,1)  ! hydrostatic stress(nincr, noel, gp, componenti (1))
!
!ccc GET INCREMENT OF sigmaH (at the beginning of the current increment)ccccc
      Write(*,*)"usdfld called"
      call getvrm('S',array,jarray,flgray,jrcd,
     $     jmac, jmtyp, matlayo, laccflg)
!
       sigmaH0 = globalSdv(noel,npt,1)
      sigmaH = third*(array(1)+array(2)+array(3))
      globalSdv(noel,npt,1)=sigmaH
!!!      hydr(kinc,noel,npt,1)=sigmaH
        dsigmaH = sigmaH-sigmaH0
! !     (dsigmaH= STATEV(1)-hydr(kinc-1,noel,npt,1))
      globalSdv(noel,npt,2)= dsigmaH
!
       !WRITE(14,*) 'STATEV(1)', STATEV(1)
       ! WRITE(14,*) 'sigmaH', sigmaH
       !WRITE(14,*) 'STATEV(2)', STATEV(2)
       !WRITE(14,*) 'dsigmaH', dsigmaH
      return
      end
