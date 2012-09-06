CC    Chris Godsalve's Hack
C     path:      $Source: /project/rc/rc1/cvsroot/rc/cntnm/src/cntnm_progr.f,v $
C     author:    $Author: jdelamer $
C     revision:  $Revision: 1.2 $
C     created:   $Date: 2011/03/29 20:16:57 $
C 
C
C  --------------------------------------------------------------------------
C |  Copyright ©, Atmospheric and Environmental Research, Inc., 2011         |
C |                                                                          |
C |  All rights reserved. This source code is part of the MT_CKD continuum   |
C |  software and is designed for scientific and research purposes.          |
C |  Atmospheric and Environmental Research, Inc. (AER) grants USER          |
C |  the right to download, install, use and copy this software              |
C |  for scientific and research purposes only. This software may be         |
C |  redistributed as long as this copyright notice is reproduced on any     |
C |  copy made and appropriate acknowledgment is given to AER. This          |
C |  software or any modified version of this software may not be            |
C |  incorporated into proprietary software or commercial software           |
C |  offered for sale.                                                       |
C |                                                                          |
C |  This software is provided as is without any express or implied          |
C |  warranties.                                                             |
C |                       (http://www.rtweb.aer.com/)                        |
C  --------------------------------------------------------------------------
C
C            Used to be     PROGRAM DRCNTNM
         subroutine driver(KntH2O, KntCO2, KntO3, KntO2, KntN2, KntRay
     &, Temp, Press , Xpath, Wave1, Wave2, DWave
     &, WaterMix, CO2Mix, NitroMix, OxyMix, OzoMix)
C
c     The mt_ckd water vapor continuum is a completely new continuum  
c     formulation based on a collision induced  component and a sub-Lorentzian 
c     line wing component.  Both the water vapor continuum coefficients and 
c     those for other molecules are constrained to agree with accurate
c     measurements of continuum absorption in the spectral regions where such
c     measurements exist.
c
c     This is an updated version of the continuum program:
c     this version provides optical depths on file CNTNM.OPTDPT as before:
c     it also provides the continuum coefficients on file  WATER.COEF
c
c     the length of the header records may vary by version:
c         in this version the WATER.COEF header information is 47 records 
c         in this version the CNTNM.OPTDT header information is 34 records 
c
c     presumably the user will want to create an input file to address
c     individual requirements
C
C
      IMPLICIT REAL*8           (V)                                     ! F00030
c
      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(5050)              ? 500060

      COMMON /ABSRBsix/ ABSRB1(5050),ABSRB2(5050)
     &,ABSRB3(5050),ABSRB4(5050),ABSRB5(5050),ABSRB6(5050)
c
      COMMON /FILHDR/ PAVE,TAVE, WK(60), WBROAD, DV ,V1 ,V2 ,NMOL
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2,GRAV,CPDAIR,AIRMWT,SECDY  
      COMMON /LAMCHN/ ONEPL,ONEMI,EXPMIN,ARGMIN                           A03020
      COMMON /XCONT/  V1C,V2C,DVC,NPTC,C(6000) 
c
c********************************************
      COMMON /cnth2o/ V1h,V2h,DVh,NPTh,Ch(5050),csh2o(5050),cfh2o(5050)
c********************************************

c------------------------------------
c for analytic derivative calculation
c note: ipts  = same dimension as ABSRB
c       ipts2 = same dimension as C
      parameter (ipts=5050,ipts2=6000)

C
      double precision Temp,Press,Xpath,WaterMix,Wave1,Wave2,DWave
      double precision CO2Mix,NitroMix,OxyMix,OzoMix
C     KntH2O, KntCO2, KntO3, KntO2, KntN2, KntRay,
      double precision KntH2O(5050),KntCO2(5050),KntO3(5050)
     &, KntO2(5050), KntN2(5050), KntRay(5050)

C                                                                           F00120
      DATA XLOSMT/2.68675E+19/                                              50024
C
      RADCN1=2.*PLANCK*CLIGHT*CLIGHT*1.E-07                                 3070
      RADCN2=PLANCK*CLIGHT/BOLTZ                                            3080
c
      icflg = -999

C
      do 2, i=1,5050
         absrb(i)=0.
         ABSRB1(i)=0.0
         ABSRB2(i)=0.0
         ABSRB3(i)=0.0
         ABSRB4(i)=0.0
         ABSRB5(i)=0.0
         ABSRB6(i)=0.0
 2    continue

      do 3, i=1,60
         wk(i)=0.
 3    continue
c

C   THIS PROGRAM CALCULATES THE CONTINUUM OPTICAL DEPTH
C         FOR AN HOMOGENEOUS LAYER
C
C   THE FOLLOWING QUANTITIES MUST BE SPECIFIED: 
C
C          PRESSURE                   PAVE (MB)
C
C          TEMPERATURE                TAVE ( K)
C
C          COLUMN AMOUNT
C            NITROGEN                 WN2    (MOLEC/CM**2)
C            OXYGEN                   WK(7)  (MOLEC/CM**2)
C            CARBON DIOXIDE           WK(2)  (MOLEC/CM**2)
C            WATER VAPOR              WK(1)  (MOLEC/CM**2)
C
C          NUMBER OF MOLECULES        NMOL
C
C          BEGINNING WAVENUMBER       V1ABS (CM-1)
C
C          ENDING WAVENUMBER          V2ABS (CM-1)
C
C          SAMPLING INTERVAL          DVABS (CM-1)
C
C          NUMBER OF VALUES           NPTABS
C
C
C   THE RESULTS ARE IN ARRAY ABSORB
C
C   NOTE THAT FOR AN ATMOSPHERIC LAYER: 
C
C            WTOT   = XLOSMT * (PAVE/1013) * (273/TAVE) * (PATH LENGTH)
C
C            WBROAD = the column amount for all species not explicitly provided
C
C            WK(M)  = (VOLUME MIXING RATIO) * (COLUMN OF DRY AIR)
C
C
C   THE FOLLOWING IS AN EXAMPLE
c
      PAVE=Press
      TAVE=Temp
      xlength=Xpath
      vmrh2o=WaterMix

c     It may be preferable to specifiy the water column directly!
c
      WTOT = XLOSMT*(PAVE/1013.)*(273./TAVE)* xlength
C
      W_dry = WTOT * (1.-VMRH2O)
c
c     ww is column of dry air;  vol. mix. ratios are based on dry air
c
c argon:
      WA     = 0.009     * W_dry
c nitrogen:
      WN2    = NitroMix      * W_dry
c oxygen:
      WK(7)  = OxyMix      * W_dry

c carbon dioxide:
      WK(2)  = CO2Mix  * W_dry

c Ozone:
      WK(3)  = OzoMix  * W_dry

c water vapor:
      if (abs(vmrh2o-1.) .lt. 1.e-05) then
         wk(1) = wtot
      else
         WK(1) = VMRH2O * W_dry
      endif
C
      WBROAD=WN2+WA
c
      NMOL = 7
c
      V1ABS = Wave1
      V2ABS = Wave2 
 
      DVABS = DWave
c ..........................................................
c      write (*,*) '  v1abs,  v2abs,  dvabs  '
c      read  (*,*)    v1abs,  v2abs,  dvabs
c ..........................................................

      NPTABS =  1. + (v2abs-v1abs)/dvabs


      jrad=0
      jrad=1
c
      v1 = v1abs
      v2 = v2abs
c
      CALL CONTNM(JRAD)

      do 250 i=1,5050
         KntH2O(i)=ABSRB1(i)
         KntCO2(i)=ABSRB2(i)
         KntO3(i)=ABSRB3(i)
         KntO2(i)=ABSRB4(i)
         KntN2(i)=ABSRB5(i)
         KntRay(i)=ABSRB6(i)
250   continue

C

c
      xkt=tave/radcn2
c
      END 

c*******
      Include 'BLOCKS.f'

c*******
      Include 'contnm.f'
