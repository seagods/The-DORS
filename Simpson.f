      Subroutine Simpson(altitude,temper,pressure,T,P,GasPPM
     &,Gass,ENNGAS,AcoeffsGas,Adummy,G,G2,GPPM,tau1Gas,mbar1Gas
     &,mbarphi1Gas,mass1Gas,tau2Gas,mbar2Gas,mbarphi2Gas,mass2Gas
     &,tau3Gas,mbar3Gas,mbarphi3Gas,mass3Gas,tau4Gas,mbar4Gas
     &,mbarphi4Gas,mass4Gas,tau1G,tau2G,tau3G,tau4G,MWGAS,OzoneR2
     &,H2OR3,legpol,legs,xx,heights,Temp,Press,Dens,up,str,tro,bound
     &,ENN,Humids,weight1,weight2,weight3,weight4,linints,T2,P2,upext
     &,strext,str2,troext,tr2,boundex,b2,AeroDryA,AeroDryE,AeroHumA
     &,AeroHumE,tau1R,tau2R,tau3R,tau4R,tau1A,tau2A,tau3A,tau4A
     &,taumix1,taumix2,taumix3,taumix4,omegmix1,omegmix2,omegmix3
     &,omegmix4,omegtau1,omegtau2,omegtau3,omegtau4,tau1,tau2,tau3
     &,tau4,ssalb1,ssalb2,ssalb3,ssalb4,omeg01,omeg02,omeg03,omeg04
     &,Haero,upperprofs550,BHaze,Hum,H2,Hum1,Hum2,Hum3,Hum4,Waves
     &,AeroHumidExt,AeroHumidAbs,AeroDryExt,AeroDryAbs,AA2,AE2
     &,WavePID,Angles,PhaseArray,moments3a,moments4a,moments1
     &,moments2,moments3,moments4,MOMS1,MOMS2,MOMS3,MOMS4,CloudC
     &,aitch,tauheights,taucloud,ssacloud,id_array,icloud,idata
     &,izeds,nprofs,igas,ilambda,ihum,iphasf,iangle,iwave,imodH
     &,imodD,iuprofs,ibprofs,Nh1,Nh2,Nh3,Nh4,nsplit1max,nsplit2max
     &,nsplit3max,nsplit4max,idiv1,idiv2,idiv3,idiv4,iheight,num
     &,istop,Ifour,imax,nmix) 
C
      integer idata,izeds,nprofs,igas,ilambda,ihum,nlayer
     &,iphasf,iangle,iwave,imodH,imodD,iuprofs,ibprofs,Nh1,Nh2,Nh3,Nh4
     &,nsplit1max,nsplit2max,nsplit3max,nsplit4max
     &,idiv1,idiv2,idiv3,idiv4,iheight,num,istop,ifour
     &,imax,nmix,iaddc
C     *********************************************************************
C                   DECLARATIONS
C     *********************************************************************
      real*8 altitude(idata),temper(idata,nprofs),pressure(idata,nprofs)
     &,T(idata),P(idata)
C      Gasseous absorption
       real*8 GasPPM(idata,igas),GPPM(idata),G2(idata),G(iheight) 
     &,Gass(iheight,igas),ENNGAS(iheight,igas) 
       real*8 AcoeffsGas(8,igas),Adummy(8),Tref
     &,mbar1Gas(nsplit1max,igas),mbarphi1Gas(nsplit1max,igas)
     &,mass1Gas(nsplit1max,igas),tau1Gas(nsplit1max,igas) 
     &,mbar2Gas(nsplit2max,igas),mbarphi2Gas(nsplit2max,igas)
     &,mass2Gas(nsplit2max,igas),tau2Gas(nsplit2max,igas) 
     &,mbar3Gas(nsplit3max,igas),mbarphi3Gas(nsplit3max,igas)
     &,mass3Gas(nsplit3max,igas),tau3Gas(nsplit3max,igas) 
     &,mbar4Gas(nsplit4max,igas),mbarphi4Gas(nsplit4max,igas)
     &,mass4Gas(nsplit4max,igas),tau4Gas(nsplit4max,igas) 
       real*8 tau1G(nsplit1max),tau2G(nsplit2max),tau3G(nsplit3max)
     &,tau4G(nsplit4max),Xozo,XH2O
     &,MWGAS(igas),OzoneR2(102),H2OR3(15)
      real*8 legpol(imax),legs(istop,iangle),xx(iangle),out
C                 
C     Form continuous profiles by interpolation                          
      double precision 
     &heights(iheight),Temp(iheight),Press(iheight),Dens(iheight) 
     &,up(iheight),str(iheight),tro(iheight),bound(iheight)
     &,sigma,ENN(iheight),Humids(iheight)
     &,weight1(iheight,nsplit1max),weight2(iheight,nsplit2max)
     &,weight3(iheight,nsplit3max),weight4(iheight,nsplit4max)
     &,aitch(nsplit1max+nsplit2max+nsplit3max+nsplit4max+1)
     &,tauheights(nsplit1max+nsplit2max+nsplit3max+nsplit4max+1)
     &,depth,vis(5),const
C
       integer linints(iheight)
C             
       real*8 TD,T2(idata),P2(idata),upext(izeds)
     &,strext(izeds),str2(izeds),troext(izeds),tr2(izeds)
     &,boundex(izeds),b2(izeds),AeroDryE(ilambda),AeroDryA(ilambda)
     &,AeroHumE(ilambda),AeroHumA(ilambda),exfac,abfac,range1
     &,range2,zstep1,zstep2,zstep3,zstep4,visb,vist
     &,tau1R(nsplit1max),tau2R(nsplit2max),tau3R(nsplit3max)
     &,tau4R(nsplit4max),ssalb1(nsplit1max),ssalb2(nsplit2max)
     &,ssalb3(nsplit3max),ssalb4(nsplit4max),Hum1(nsplit1max)
     &,Hum2(nsplit2max),Hum3(nsplit3max),Hum4(nsplit4max)
     &,tau1A(nsplit1max),tau2A(nsplit2max),tau3A(nsplit3max)
     &,tau4A(nsplit4max),taucloud(nsplit4max+nsplit3max)
     &,ssacloud(nsplit4max+nsplit3max)
     &,taumix1(nsplit1max,nmix),taumix2(nsplit2max,nmix)
     &,taumix3(nsplit3max,nmix),taumix4(nsplit4max,nmix)
     &,omegmix1(nsplit1max,nmix),omegmix2(nsplit2max,nmix)
     &,omegmix3(nsplit3max,nmix),omegmix4(nsplit4max,nmix)
     &,omegtau1(nsplit1max),omegtau2(nsplit2max)
     &,omegtau3(nsplit3max),omegtau4(nsplit4max),tau1(nsplit1max)
     &,tau2(nsplit2max),tau3(nsplit3max),tau4(nsplit4max),weighting 
     &,omeg01(nsplit1max),omeg02(nsplit2max),omeg03(nsplit3max)
     &,omeg04(nsplit4max),tau_Aero,tau_Ray
C     calculate Rayleigh cross section profile and Rayleigh             
C     Optical Depth                                                     
      real*8 R,AV,MWA,lambda,delta     
      real*8 rfact,pi,NSTP,NSP15,recip2,mconst
      real*8 abszero,stanpress
C     Aerosol properties as a function of height                       
      real*8 Haero(izeds),upperprofs550(izeds,iuprofs)
      real*8 BHaze(izeds,ibprofs),Hum(izeds)
      real*8 H2(izeds),thik1,thik2,thik3,thik4
      real*8 Waves(ilambda),AeroHumidExt(ilambda,imodH,ihum)
      real*8 AeroHumidAbs(ilambda,imodH,ihum),AeroDryExt(ilambda,imodD)       
      real*8 AeroDryAbs(ilambda,imodD),AA2(ilambda),AE2(ilambda)
C     Phase function stuff                                             
      real*8 WavePid(iwave),Angles(iangle),PhaseArray(300,iphasf) 
     &,moments3a(imax,ihum),moments4a(imax,ihum)
     &,moments1(imax,2),moments2(imax,2)
     &,moments3(imax,nsplit3max,2),moments4(imax,nsplit4max,2)
     &,MOMS1(imax,nsplit1max),MOMS2(imax,nsplit2max)
     &,MOMS3(imax,nsplit3max),MOMS4(imax,nsplit4max)
     &,CloudC(imax,4)
      integer ifail,n,iw1,iw2
      integer id_array(iwave,26),icloud(nsplit4max+nsplit3max)
     &,iatm,i,k,ix1,ix2,ix3,ix4,itype1
     &,itype2,itypeu,itypes,itypet,itypeb,iupphase,istratphase
     &,iboundphase,itropphase,iocean,itypetro,itypetr
      integer nsplit1,nsplit2,nsplit3,nsplit4,iorder1,iorder2
     &,iorder3,iorder4,j,iray,Inu,IclassW1,IclassW2,IclassW3 

       real*8 groundtemp,groundpress,temppress,heightplane

       integer ihumid,igroundorprof,iplane,inu1,inu2,inu2b,inu3
     &,inu3b,iplanesplit1,iplanesplit2,iplanesplit3,iplanesplit4
     &,istart1,istart2,istart3,istart4,nn,iplanelayer

       character*64 filephase

       logical INRANGE1,INRANGE2,INRANGE3
C
      call Questioner(
     &vist,visb,lambda,nsplit1,nsplit2,nsplit3,nsplit4
     &,iatm,iray,itypeu,itypes,itypet,itypeb,iocean
     &,ihumid,igroundorprof,groundTemp,groundpress
     &,iplane,heightplane,iaddc)

C      Get phase function model identifier
       if(itypeu.eq.13.or.itypeu.eq.16)then
         iupphase=26
       else
         if(itypeu.eq.14)iupphase=22
         if(itypeu.eq.15)iupphase=23
       endif

      vis(1)=50.
      vis(2)=23.
      vis(3)=10.
      vis(4)=5.
      vis(5)=2.
C
      IFAIL=1
      pi=dacos(-1d0)                                                     
C     zero Celcius in Kelvin scale                                      
      abszero=273.15                                                    
C     Universal gas constant (per mole)
      R=8.33143                                                          
C     Avagrado's constant (per mole)                                
      AV=6.022d23                                                       
C     Molecular weight of air                                           
      MWA= 28.964     
C     Depolarisation factor (use "new improved one from 6S") 
      delta=0.029                                                       
C     standard pressure                                                  
      stanpress=1.01325d5   
C     *********************************************************************
C                POINTS FOR QUADRATURE RULE AND INTERPOLATION
C     *********************************************************************
      ix4=idiv4-idiv3+1
      ix3=idiv4-idiv2+1
      ix2=idiv4-idiv1+1
      ix1=idiv4
      do 4 i=1,iheight
        if(i.le.ix4)then
        heights(i)=dble(i-1)/dble(2**Nh4)*2000.
        endif
        if((i.gt.ix4).and.(i.le.ix3))then
        heights(i)=2000.+dble(i-ix4)/dble(2**Nh3)*8000.
        endif
        if((i.gt.ix3).and.(i.le.ix2))then
        heights(i)=10000.+dble(i-ix3)/dble(2**Nh2)*20000.
        endif
        if((i.gt.ix2).and.(i.le.ix1))then
        heights(i)=30000.+dble(i-ix2)/dble(2**Nh1)*90000.
        endif
4     continue

C     This has given us all the Simpson's rule points we need
C     for integrating over each region, But we want to
C     split each region into layers and find the Aerosol
C     optical depth of each.
      zstep1=90000./2.**Nh1
      zstep2=20000./2.**Nh2
      zstep3=8000./2.**Nh3
      zstep4=2000./2.**Nh4
      thik1=90000.0/dble(nsplit1)
      thik2=20000.0/dble(nsplit2)
      thik3=8000.0/dble(nsplit3)
      thik4=2000.0/dble(nsplit4)
C
C     Order of Simpsons Rule within each layer
      iorder1=2**Nh1/nsplit1+1
      iorder2=2**Nh2/nsplit2+1
      iorder3=2**Nh3/nsplit3+1
      iorder4=2**Nh4/nsplit4+1
C     weights for integrations over upper atmosphere
C      **********************************************
C     Note, that since there are an odd number of weights
C     and its the same number within each subregion of a 
C     LOWTRAN aerosol region, 
C     and that they are symmetric wrt exchange of order, thus
C     can take the liberty of writing the weights as if we
C     were taking our sublayers from top to bottom 
C     ,(Not The LOWTRAN REGIONS),
C     and then reverse the order on writing to PHASEDATA.
      do 12 i=1,nsplit1
        do  8 j=1,iorder1
          if(mod(j,2).eq.0)then
            weight1(ix2+j-i+(i-1)*iorder1,i)=4d0/3d0
          else
            weight1(ix2+j-i+(i-1)*iorder1,i)=2d0/3d0
          endif
          if((j.eq.1).or.(j.eq.iorder1))
     &weight1(ix2+j-i+(i-1)*iorder1,i)=1d0/3d0        
8       continue
12     continue
C     weights for integrations over stratosphere
      do 20 i=1,nsplit2
        do 16 j=1,iorder2
          if(mod(j,2).eq.0)then
            weight2(ix3+j-i+(i-1)*iorder2,i)=4d0/3d0
          else
            weight2(ix3+j-i+(i-1)*iorder2,i)=2d0/3d0
          endif
          if((j.eq.1).or.(j.eq.iorder2))
     &weight2(ix3+j-i+(i-1)*iorder2,i)=1d0/3d0 
16       continue
20     continue
C     weights for integrations over upper troposphere
      do 28 i=1,nsplit3
        do  24 j=1,iorder3
          if(mod(j,2).eq.0)then
            weight3(ix4+j-i+(i-1)*iorder3,i)=4d0/3d0
          else
            weight3(ix4+j-i+(i-1)*iorder3,i)=2d0/3d0
          endif
          if((j.eq.1).or.(j.eq.iorder3))
     &weight3(ix4+j-i+(i-1)*iorder3,i)=1d0/3d0  
24       continue
28     continue
C     weights for integrations over boundary layer
      do 36 i=1,nsplit4
        do 32 j=1,iorder4
          if(mod(j,2).eq.0)then
            weight4(1+j-i+(i-1)*iorder4,i)=4d0/3d0
          else
            weight4(1+j-i+(i-1)*iorder4,i)=2d0/3d0
          endif
          if((j.eq.1).or.(j.eq.iorder4))
     &weight4(1+j-i+(i-1)*iorder4,i)=1d0/3d0        
32       continue
36     continue
                                                                       
C     wavelength dependence of real part of refractive index            
C     LIOU (Atmospheric Radiation) p79 and appendix D eqn.17            
                                                                        
C     *********************************************************** 
C            LOWTRAN PROFILES   
C     *********************************************************** 
C     Load LOWTRAN altitudes for temperature, pressure and water vapour    
C     density profiles                                                  
      Call Alts(altitude)                                               
C     Load Lowtran Temperature Profiles                                 
      Call Temps(temper)                                                
C     Load LOWTRAN Pressure Profiles                                    
      Call Presss(pressure)
      Call GetGas(GasPPM,iatm)

C      Load the 34 altitudes which define the aerosol profile (see PRF37
C     (Different from the heights defining temperature and Pressure)
C      These are also used  for user humidity profiles    
       Call HeightHaze(Haero) 
C      convert to metre scale
       do 40 i=1,izeds
         Haero(i)=Haero(i)*1000.
40     continue
        if(ihumid.eq.2)then
         open(1,file="MyHumidity.data",status="old")
         do 41 i=1,izeds
           read(1,*)Hum(i)
41       continue
         close(1,status="keep")
C        SPline to finer resolution
         TD=1d35
         call spline(Haero,Hum,izeds,TD,TD,H2)
         do 42 i=1,iheight
            call splint(Haero,Hum,H2,izeds,heights(i),Humids(i))
C           can get a bit of oscilation at high altitudes
C           where the RH is very small - reset to zero
            if(Humids(i).lt.0)Humids(i)=0.
42       continue
       endif
C      Pick out corresponding atmosphere                                
       do 43 i=1,idata                                                  
         T(i)=Temper(i,iatm)   
C        convert to from mbar to newtons per m**2 
         P(i)=Pressure(i,iatm)*100. 
C        convert to metres
         altitude(i)=altitude(i)*1000.
43     continue 
        if(igroundorprof.eq.1)then
         tempPress=P(1)
         do 44 i=1,idata
           P(i)=P(i)*groundPress/tempPress
44        continue
         T(1)=groundTemp
         T(2)=groundTemp+(T(3)-T(1))/2.
        endif
        if(igroundorprof.eq.2)then
         open(1,file="MyPressures.data",status="old")
         do 45 i=1,izeds
           read(1,*)P(i)
45        continue
         close(1,status="keep")
         open(1,file="MyTemperatures.data",status="old")
         do 46 i=1,idata
           read(1,*)T(i)
46        continue
         close(1,status="keep")
        endif
C      Spline to finer resolution                                       
       TD=1d35                      
       call spline(altitude,T,idata,TD,TD,T2)                         
       call spline(altitude,P,idata,TD,TD,P2)
       do 47 i=1,iheight                                                
          call splint(altitude,T,T2,idata,heights(i),Temp(i))           
          call splint(altitude,P,P2,idata,heights(i),Press(i))
47     continue 
         call splinegas(altitude
     &,GasPPM,Gass,G,G2,GPPM,heights,idata,iheight,igas)
         MWGAS(1)=18.
         MWGAS(2)=48.
         MWGAS(3)=44.
         MWGAS(4)=32.
         MWGAS(5)=44.
         MWGAS(6)=16.
         MWGAS(7)=28.
         do 50 i=1,iheight
         ENN(i)=Press(i)/R/Temp(i)*AV
C        factor of 1000. to get Kg  
         Dens(i)=ENN(i)*MWA/AV/1000.         
C        particles per m**3 * mass per mol = density
C        divide by 1000000 because Gass is in ppm
C        divide by 1000. because AV is in g and want Kg
C        Units of ENNH2O etc Kg/m**2
         if(ihumid.eq.2)call 
     & WaterVap(Press,Temp,Gass(1,1),Humids,AV,MWGAS(1),MWA
     &,iheight,ihumid)
         do 49 ig=1,igas
            ENNGAS(i,ig)=
     & ENN(i)*(Gass(i,ig)/1000000.)*MWGAS(ig)/AV/1000.
49       continue
50       continue
         if(ihumid.eq.1)call 
     & WaterVap(Press,Temp,Gass(1,1),Humids,AV,MWGAS(1),MWA
     &,iheight,ihumid)
        if(iaddc.eq.1)then
C       find where to add cloud
        call cloudy(heights,Humids,taucloud,ssacloud
     &,iheight,icloud,ix4,iorder3,iorder4,nsplit3,nsplit4)
C       recalculate water vapour mass in cloud layers
       ihumid=2
       call WaterVap(Press,Temp,Gass(1,1),Humids,AV
     &,MWGAS(1),MWA,iheight,ihumid)
        endif
        do 51 i=1,nsplit1
           taumix1(i,3)=0.
           omegmix1(i,3)=0.
51      continue
        do 52 i=1,nsplit2
           taumix2(i,3)=0.
           omegmix2(i,3)=0.
52      continue
        do 53 i=1,nsplit3
           taumix3(i,3)=0.
           omegmix3(i,3)=0.
           if(icloud(nsplit4+i).gt.0)then
           taumix3(i,3)=taucloud(i+nsplit4)
           omegmix3(i,3)=ssacloud(i+nsplit4)
           omegtau3(i)=omegtau3(i)
     &+taucloud(i+nsplit4)*ssacloud(i+nsplit4)
           tau3(i)=tau3(i)+taucloud(i+nsplit4)
           endif
53      continue
        do 54 i=1,nsplit4
           taumix4(i,3)=0.
           omegmix4(i,3)=0.
           if(icloud(i).gt.0)then
           taumix4(i,3)=taucloud(i)
           omegmix4(i,3)=ssacloud(i)
           omegtau4(i)=omegtau4(i)+taucloud(i)*ssacloud(i)
           tau4(i)=tau4(i)+taucloud(i)
           endif
54      continue
      open(12,file="CLOUDCDATA1",status="old")
      read(12,*,end=55)(CloudC(i,1),i=1,imax)
55    close(12,status="keep")
      open(12,file="CLOUDCDATA2",status="old")
      read(12,*,end=56)(CloudC(i,2),i=1,imax)
56    close(12,status="keep")
      open(12,file="CLOUDCDATA3",status="old")
      read(12,*,end=57)(CloudC(i,3),i=1,imax)
57    close(12,status="keep")
      open(12,file="CLOUDCDATA4",status="old")
      read(12,*,end=58)(CloudC(i,4),i=1,imax)
58    close(12,status="keep")

      do 67 i=1,nsplit4
      do 66 k=1,imax
         if(icloud(i).gt.0)then
         moments4(k,i,2)=CloudC(k,icloud(i))
         else
         moments4(k,i,2)=0.
         endif
66    continue
67    continue

      do 69 i=1,nsplit3
      do 68 k=1,imax
         if(icloud(i+nsplit4).gt.0)then
         moments3(k,i,2)=CloudC(k,icloud(i+nsplit4))
         else
         moments3(k,i,2)=0.
         endif
68    continue
69    continue



C
         call GetOZO2(OzoneR2)
         call GetH2O3(H2OR3)
C         a factor 10 to small
          Inu=int(1d0/lambda*1000)
C         then multiply by 10 to get steps of 10 per cm
          Inu=Inu*10
          wavynum=1.0/lambda*10000.0
          open(17,file="Wnum",status="unknown")
          write(17,*)WavyNum
          close(17,status="keep")

C     *********************************************************************
C                 Rayleigh Extinction Coefficients
C     *********************************************************************
C      moles per cubic metre                                            
C      factor of 10**24 from lambda being in microns and ENN in per m**3
C      refracive index -1                                              
       rfact=1d24*                                                      
     &8d0*pi**3/3d0/lambda**4*(6d0+3d0*delta)/(6d0-7d0*delta)
C      Number density of molecules at STP                               
       NSTP=stanpress/R/abszero*AV                                      
C      Edlens formula defined at 15 degrees celcius and standard pressur
C      NOT AT STP (see penndorf )                                        
       NSP15=stanpress/R/(abszero+15)*AV                                
       recip2=1d0/lambda**2                                             
C      Edlens Formula                                                   
       mconst=                                                         
     &(6432.8d0+2949810d0/(146d0-recip2)+25540d0/(41d0-recip2))/1d8     
C     Calculate Cross Section as a function of height in order to calculate
C     Rayleigh Scattering.      
         sigma=rfact*(4.*mconst**2+4.*mconst**3+mconst**4)/NSP15**2
         IclassW1=0
         IclassW2=0
         IclassW3=0       
         INRANGE1=.TRUE.                     
         if(Inu.lt.2500)INRANGE1=.FALSE.  
         if(Inu.gt.2500.and.Inu.lt.5060)IclassW1=1
         if(Inu.ge.5060.and.Inu.lt.7620)IclassW1=2
         if(Inu.ge.7620.and.Inu.lt.10180)IclassW1=3
         if(Inu.ge.10180.and.Inu.lt.12740)IclassW1=4
         if(Inu.ge.12740.and.Inu.lt.15300)IclassW1=5
         if(Inu.ge.15300.and.Inu.lt.17860)IclassW1=6
         if(Inu.gt.17860)INRANGE1=.FALSE.  
C        convert from wavenumber to look up table pos'n
         if(IclassW1.eq.1)Inu1=(Inu-2500)/10+1
         if(IclassW1.eq.2)Inu1=(Inu-5060)/10+1
         if(IclassW1.eq.3)Inu1=(Inu-7620)/10+1
         if(IclassW1.eq.4)Inu1=(Inu-10180)/10+1
         if(IclassW1.eq.5)Inu1=(Inu-12740)/10+1
         if(IclassW1.eq.6)Inu1=(Inu-15300)/10+1
C        for ozone bands
         INRANGE2=.TRUE.                     
         if(Inu.lt.13000)INRANGE2=.FALSE.  
         if(Inu.ge.13000.and.Inu.le.24200)iclassW2=1
         if(Inu.gt.24200.and.Inu.lt.27500.)INRANGE2=.FALSE.
         if(Inu.ge.27500.and.Inu.le.50000)iclassW2=2
         if(Inu.gt.50000)INRANGE2=.FALSE.
         if(INRANGE2)then
            if(IclassW2.eq.1)Inu2=(Inu-13000)/200+1
            if(IclassW2.eq.2)Inu2=(Inu-27500)/500+1
            Inu2b=Inu2+1
            if(IclassW2.eq.1)then
               Xozo=(10000./dble(13000+200*(Inu2-1))-lambda)
     &/(10000./dble(13000+200*(Inu2-1))
     &-10000./dble(13000+200*(Inu2b-1)))
            endif
            if(IclassW2.eq.2)then
               Xozo=(10000./dble(27500+500*(Inu2-1))-lambda)
     &/(10000./dble(27500+500*(Inu2-1))
     &-10000./dble(27500+500*(Inu2b-1)))
            Inu2=Inu2+54
            Inu2b=Inu2b+54
            endif
         endif
C        for water vapour at this range
         if(Inu.ge.2350.and.Inu.le.2420)then
            INRANGE3=.TRUE.
            iclassW3=1
         else
            INRANGE3=.FALSE.
         endif
         if(INRANGE3)then
         Inu3=(Inu-2350)/5+1
         Inu3b=Inu3+1
            XH2O=(10000./dble(2350+5*(Inu3-1))-lambda)
     &/(10000./dble(2350+5*(Inu3-1))
     &-10000./dble(2350+5*(Inu3b-1)))
          endif
         Tref=250.0
         do 73 i=1,8
           Adummy(i)=0.
73       continue
         if(INRANGE1)then
C        Get coeffs  water vapour 
         Call GetH2O(Adummy,IclassW1,Inu1)
         do 74 i=1,8
           AcoeffsGas(i,1)=Adummy(i)
           Adummy(i)=0.
74       continue
C        Get coeffs for ozone
         if(IclassW1.eq.1)then
         Call GetOZO(Adummy,IclassW1,Inu1)
         endif
         do 75 i=1,8
           AcoeffsGas(i,2)=Adummy(i)
           Adummy(i)=0.
75       continue
C        Get coeffs for carbon dioxide
         if(IclassW1.le.3)then
         Call GetCO2(Adummy,IclassW1,Inu1)
         endif
         do 76 i=1,8
           AcoeffsGas(i,3)=Adummy(i)
           Adummy(i)=0.
76       continue
C        Get coeffs for oxygen
         if(IclassW1.ge.3)then
         Call GetOXY(Adummy,IclassW1,Inu1)
         endif
         do 77 i=1,8
           AcoeffsGas(i,4)=Adummy(i)
           Adummy(i)=0.
77       continue
C        Get coeffs for nitreous oxide
         Call GetNOX(Adummy,IclassW1,Inu1)
         do 78 i=1,8
           AcoeffsGas(i,5)=Adummy(i)
           Adummy(i)=0.
78       continue
C        Get coeffs for methane
         Call GetMET(Adummy,IclassW1,Inu1)
         do 79 i=1,8
           AcoeffsGas(i,6)=Adummy(i)
           Adummy(i)=0.
79       continue
C        Get coeffs for carbon monoxide
         Call GetCMO(Adummy,IclassW1,Inu1)
         do 80 i=1,8
           AcoeffsGas(i,7)=Adummy(i)
           Adummy(i)=0.
80       continue
         endif

CC        upper atmosphere Optical Depths  
       do 85 ig=1,igas 
       do 84 i=1,nsplit1
         if(ig.eq.1)tau1R(i)=0.
         mbar1Gas(i,ig)=0.
         mbarphi1Gas(i,ig)=0.
         mass1Gas(i,ig)=0.
         do 81 j=ix2+(i-1)*(iorder1-1),ix2+i*(iorder1-1)
           if(ig.eq.1)tau1R(i)=
     &tau1R(i)+zstep1*weight1(j,i)*sigma*ENN(j)
           if(INRANGE1)then
           mbar1Gas(i,ig)=mbar1Gas(i,ig)+zstep1*weight1(j,i)
     &*exp(AcoeffsGas(3,ig)*(Temp(j)-Tref)
     &+AcoeffsGas(4,ig)*(Temp(j)-Tref)**2)*ENNGAS(j,ig)
           mbarphi1Gas(i,ig)=mbarphi1Gas(i,ig)+zstep1*weight1(j,i)
     &*exp(AcoeffsGas(5,ig)*(Temp(j)-Tref)
     &+AcoeffsGas(6,ig)*(Temp(j)-Tref)**2)*ENNGAS(j,ig)*Press(j)
     &/stanpress
         endif
         mass1Gas(i,ig)=mass1Gas(i,ig)+zstep1*weight1(j,i)*ENNGAS(j,ig)
81       continue
C        convert to  g per cm**2
         mbar1Gas(i,ig)=mbar1Gas(i,ig)/10.
         mbarphi1Gas(i,ig)=mbarphi1Gas(i,ig)/10.
         mass1Gas(i,ig)=mass1Gas(i,ig)/10.
         if(INRANGE1)then
C        Optical Depths for Water
         if(ig.eq.1)then
C        Goody Band model for H2O 
         if(AcoeffsGas(2,ig).ne.0.0.and.mbarphi1Gas(i,ig).ne.0.0)
     & tau1Gas(i,ig)=mbar1Gas(i,ig)*AcoeffsGas(1,ig)
     &/sqrt(1.0+AcoeffsGas(1,ig)/AcoeffsGas(2,ig)
     &*mbar1Gas(i,ig)**2/mbarphi1Gas(i,ig)  )
         endif
         if(ig.gt.1)then
         if(AcoeffsGas(2,ig).ne.0.0.and.mbarphi1Gas(i,ig).ne.0.0)
     & tau1Gas(i,ig)=
     & mbarphi1Gas(i,ig)/2./mbar1Gas(i,ig)*AcoeffsGas(2,ig)
     &*(sqrt(  1.0+4.*AcoeffsGas(1,ig)/AcoeffsGas(2,ig)
     &*mbar1Gas(i,ig)**2/mbarphi1Gas(i,ig)  ) -1.)
           if(ig.eq.4)then
           if(tau1Gas(i,4).gt.1d-9)then
C             Reduce total molecular scattering cross section
C             as oxygen absorption cross section grows
              tau1R(i)=(0.8+0.2/(1.+tau1Gas(i,4)/.2))*tau1R(i)
           endif
           endif
         endif
         endif
         if(INRANGE2)then
         if(ig.eq.2)then
C        Dobson Units for this one?
         tau1Gas(i,ig)=mass1Gas(i,ig)*1000./(0.048/0.0224)
     &*(OzoneR2(Inu2)+Xozo*(OzoneR2(Inu2b)-OzoneR2(Inu2)) )
         endif
         endif
         if(INRANGE3)then
         if(ig.eq.1)then
         tau1Gas(i,ig)=mass1Gas(i,ig)
     &*(H2OR3(Inu2)+XH2O*(H2OR3(Inu2b)-H2OR3(Inu2)) )
         endif
         endif
       tau1G(i)=tau1G(i)+tau1Gas(i,ig)
84     continue
85     continue
       do 93 ig=1,igas
       do 92 i=1,nsplit2
         if(ig.eq.1)tau2R(i)=0.
         mbar2Gas(i,ig)=0.
         mbarphi2Gas(i,ig)=0.
         mass2Gas(i,ig)=0.
         do 88 j=ix3+(i-1)*(iorder2-1),ix3+i*(iorder2-1)
           if(ig.eq.1)tau2R(i)=
     &tau2R(i)+zstep2*weight2(j,i)*sigma*ENN(j)
           if(INRANGE1)then
           mbar2Gas(i,ig)=mbar2Gas(i,ig)+zstep2*weight2(j,i)
     &*exp(AcoeffsGas(3,ig)*(Temp(j)-Tref)
     &+AcoeffsGas(4,ig)*(Temp(j)-Tref)**2)*ENNGAS(j,ig)
           mbarphi2Gas(i,ig)=mbarphi2Gas(i,ig)+zstep2*weight2(j,i)
     &*exp(AcoeffsGas(5,ig)*(Temp(j)-Tref)
     &+AcoeffsGas(6,ig)*(Temp(j)-Tref)**2)*ENNGAS(j,ig)*Press(j)
     &/stanpress
           endif
         mass2Gas(i,ig)=mass2Gas(i,ig)+zstep2*weight2(j,i)*ENNGAS(j,ig)
88       continue
C        convert to  g per cm**2
         mbar2Gas(i,ig)=mbar2Gas(i,ig)/10.
         mbarphi2Gas(i,ig)=mbarphi2Gas(i,ig)/10.
         mass2Gas(i,ig)=mass2Gas(i,ig)/10.
       if(INRANGE1)then
         if(ig.eq.1)then
C        Goody Band model for H2O 
         if(AcoeffsGas(2,ig).ne.0.0.and.mbarphi2Gas(i,ig).ne.0.0)
     & tau2Gas(i,ig)=mbar2Gas(i,ig)*AcoeffsGas(1,ig)
     &/sqrt(1.0+AcoeffsGas(1,ig)/AcoeffsGas(2,ig)
     &*mbar2Gas(i,ig)**2/mbarphi2Gas(i,ig)  )
         endif
        if(ig.gt.1)then
        if(AcoeffsGas(2,ig).ne.0.0.and.mbarphi2Gas(i,ig).ne.0.0)
     & tau2Gas(i,ig)=
     & mbarphi2Gas(i,ig)/2./mbar2Gas(i,ig)*AcoeffsGas(2,ig)
     &*(sqrt(  1.0+4.*AcoeffsGas(1,ig)/AcoeffsGas(2,ig)
     &*mbar2Gas(i,ig)**2/mbarphi2Gas(i,ig)  ) -1.)
           if(ig.eq.4)then
           if(tau2Gas(i,4).gt.1d-9)then
C             Reduce total molecular scattering cross section
C             as oxygen absorption cross section grows
              tau2R(i)=(0.8+0.2/(1.+tau2Gas(i,4)/.2))*tau2R(i)
           endif
           endif
         endif
       endif
       if(INRANGE2)then
         if(ig.eq.2)then
C        Dobson Units for this one?
         tau2Gas(i,ig)=mass2Gas(i,ig)*1000./(0.048/0.0224)
     &*(OzoneR2(Inu2)+Xozo*(OzoneR2(Inu2b)-OzoneR2(Inu2)) )
         endif
       endif
       if(INRANGE3)then
         if(ig.eq.1)then
         tau2Gas(i,ig)=mass2Gas(i,ig)
     &*(H2OR3(Inu2)+XH2O*(H2OR3(Inu2b)-H2OR3(Inu2)) )
         endif
       endif
       tau2G(i)=tau2G(i)+tau2Gas(i,ig)
92     continue
93     continue
       do 101 ig=1,igas
       do 100 i=1,nsplit3
         if(ig.eq.1)tau3R(i)=0d0
         mbar3Gas(i,ig)=0.
         mbarphi3Gas(i,ig)=0.
         mass3Gas(i,ig)=0.
         do 96 j=ix4+(i-1)*(iorder3-1),ix4+i*(iorder3-1)
         if(ig.eq.1)tau3R(i)=
     & tau3R(i)+zstep3*weight3(j,i)*sigma*ENN(j)
         if(INRANGE1)then
           mbar3Gas(i,ig)=mbar3Gas(i,ig)+zstep3*weight3(j,i)
     &*exp(AcoeffsGas(3,ig)*(Temp(j)-Tref)
     &+AcoeffsGas(4,ig)*(Temp(j)-Tref)**2)*ENNGAS(j,ig)
           mbarphi3Gas(i,ig)=mbarphi3Gas(i,ig)+zstep3*weight3(j,i)
     &*exp(AcoeffsGas(5,ig)*(Temp(j)-Tref)
     &+AcoeffsGas(6,ig)*(Temp(j)-Tref)**2)*ENNGAS(j,ig)*Press(j)
     &/stanpress
         endif
         mass3Gas(i,ig)=mass3Gas(i,ig)+zstep3*weight3(j,i)*ENNGAS(j,ig)
96       continue
C        convert to  g per cm**2
         mbar3Gas(i,ig)=mbar3Gas(i,ig)/10.
         mbarphi3Gas(i,ig)=mbarphi3Gas(i,ig)/10.
         mass3Gas(i,ig)=mass3Gas(i,ig)/10.
       if(INRANGE1)then
         if(ig.eq.1)then
C        Goody Band model for H2O 
         if(AcoeffsGas(2,ig).ne.0.0.and.mbarphi3Gas(i,ig).ne.0.0)
     & tau3Gas(i,ig)=mbar3Gas(i,ig)*AcoeffsGas(1,ig)
     &/sqrt(1.0+AcoeffsGas(1,ig)/AcoeffsGas(2,ig)
     &*mbar3Gas(i,ig)**2/mbarphi3Gas(i,ig)  )
         endif
         if(ig.gt.1)then
         if(AcoeffsGas(2,ig).ne.0.0.and.mbarphi3Gas(i,ig).ne.0.0)
     & tau3Gas(i,ig)=
     & mbarphi3Gas(i,ig)/2./mbar3Gas(i,ig)*AcoeffsGas(2,ig)
     &*(sqrt(  1.0+4.*AcoeffsGas(1,ig)/AcoeffsGas(2,ig)
     &*mbar3Gas(i,ig)**2/mbarphi3Gas(i,ig)  ) -1.)
           if(ig.eq.4)then
           if(tau3Gas(i,4).gt.1d-9)then
C             Reduce total molecular scattering cross section
C             as oxygen absorption cross section grows
              tau3R(i)=(0.8+0.2/(1.+tau3Gas(i,4)/0.2))*tau3R(i)
           endif
           endif
         endif
       endif
       if(INRANGE2)then
         if(ig.eq.2)then
C        Dobson Units for this one?
         tau3Gas(i,ig)=mass3Gas(i,ig)*1000./(0.048/0.0224)
     &*(OzoneR2(Inu2)+Xozo*(OzoneR2(Inu2b)-OzoneR2(Inu2)) )
         endif
       endif  
       if(INRANGE3)then
         if(ig.eq.1)then
         tau3Gas(i,ig)=mass3Gas(i,ig)
     &*(H2OR3(Inu2)+XH2O*(H2OR3(Inu2b)-H2OR3(Inu2)) )
         endif
       endif
       tau3G(i)=tau3G(i)+tau3Gas(i,ig)
100    continue
101    continue
       do 109 ig=1,igas
       do 108 i=1,nsplit4
         if(ig.eq.1)tau4R(i)=0.
         mbar4Gas(i,ig)=0.
         mbarphi4Gas(i,ig)=0.
         mass4Gas(i,ig)=0.
         do 104 j=1+(i-1)*(iorder4-1),1+i*(iorder4-1)
           if(ig.eq.1)tau4R(i)=
     & tau4R(i)+zstep4*weight4(j,i)*sigma*ENN(j)
           if(INRANGE1)then
           mbar4Gas(i,ig)=mbar4Gas(i,ig)+zstep4*weight4(j,i)
     &*exp(AcoeffsGas(3,ig)*(Temp(j)-Tref)
     &+AcoeffsGas(4,ig)*(Temp(j)-Tref)**2)*ENNGAS(j,ig)
           mbarphi4Gas(i,ig)=mbarphi4Gas(i,ig)+zstep4*weight4(j,i)
     &*exp(AcoeffsGas(5,ig)*(Temp(j)-Tref)
     &+AcoeffsGas(6,ig)*(Temp(j)-Tref)**2)*ENNGAS(j,ig)*Press(j)
     &/stanpress
         endif
         mass4Gas(i,ig)=mass4Gas(i,ig)+zstep4*weight4(j,i)*ENNGAS(j,ig)
104       continue 
C        convert to  g per cm**2
         mbar4Gas(i,ig)=mbar4Gas(i,ig)/10.
         mbarphi4Gas(i,ig)=mbarphi4Gas(i,ig)/10.
         mass4Gas(i,ig)=mass4Gas(i,ig)/10.
         if(INRANGE1)then
         if(ig.eq.1)then
C        Goody Band model for H2O 
         if(AcoeffsGas(2,ig).ne.0.0.and.mbarphi4Gas(i,ig).ne.0.0)
     & tau4Gas(i,ig)=mbar4Gas(i,ig)*AcoeffsGas(1,ig)
     &/sqrt(1.0+AcoeffsGas(1,ig)/AcoeffsGas(2,ig)
     &*mbar4Gas(i,ig)**2/mbarphi4Gas(i,ig)  )
         endif
         if(ig.gt.1)then
        if(AcoeffsGas(2,ig).ne.0.0.and.mbarphi4Gas(i,ig).ne.0.0)
     & tau4Gas(i,ig)=
     & mbarphi4Gas(i,ig)/2./mbar4Gas(i,ig)*AcoeffsGas(2,ig)
     &*(sqrt(  1.0+4.*AcoeffsGas(1,ig)/AcoeffsGas(2,ig)
     &*mbar4Gas(i,ig)**2/mbarphi4Gas(i,ig)  ) -1.)   
           if(ig.eq.4)then
           if(tau4Gas(i,4).gt.1d-9)then
C             Reduce total molecular scattering cross section
C             as oxygen absorption cross section grows
              tau4R(i)=(0.8+.2/(1.+tau4Gas(i,4)/0.2))*tau4R(i)
           endif
           endif
       endif
         endif
         if(INRANGE2)then
         if(ig.eq.2)then
C        Dobson Units for this one?
         tau4Gas(i,ig)=mass4Gas(i,ig)*1000./(0.048/0.0224)
     &*(OzoneR2(Inu2)+Xozo*(OzoneR2(Inu2b)-OzoneR2(Inu2)) )
         endif
         endif  
         if(INRANGE3)then
         if(ig.eq.1)then
         tau4Gas(i,ig)=mass4Gas(i,ig)
     &*(H2OR3(Inu2)+XH2O*(H2OR3(Inu2b)-H2OR3(Inu2)) )
         endif
         endif
         tau4G(i)=tau4G(i)+tau4Gas(i,ig)
108     continue
109    continue
       tau_Ray=0.
       tau_Aero=0.

       do 112 i=1,nsplit1
        if(iray.eq.15)tau1R(i)=0.
        taumix1(i,1)=tau1R(i)+tau1G(i)
        omegmix1(i,1)=tau1R(i)/taumix1(i,1)
        write(6,*)'Effective Rayliegh SSA',omegmix1(i,1)
        tau_Ray=tau_Ray+tau1R(i)
        tau1(i)=tau1(i)+taumix1(i,1)
        omegtau1(i)=omegtau1(i)+taumix1(i,1)*omegmix1(i,1)
112     continue
       do 116 i=1,nsplit2
        if(iray.eq.15)tau2R(i)=0.
        taumix2(i,1)=tau2R(i)+tau2G(i)
        omegmix2(i,1)=tau2R(i)/taumix2(i,1)
        write(6,*)'Effective Rayliegh SSA',omegmix2(i,1)
        tau_Ray=tau_Ray+tau2R(i)
         tau2(i)=tau2(i)+taumix2(i,1)
        omegtau2(i)=omegtau2(i)+taumix2(i,1)*omegmix2(i,1)
116     continue
       do 120 i=1,nsplit3
        if(iray.eq.15)tau3R(i)=0.
        taumix3(i,1)=tau3R(i)+tau3G(i)
        omegmix3(i,1)=tau3R(i)/taumix3(i,1)
        write(6,*)'Effective Rayliegh SSA',omegmix3(i,1)
        tau_Ray=tau_Ray+tau3R(i)
        tau3(i)=tau3(i)+taumix3(i,1)
        omegtau3(i)=omegtau3(i)+taumix3(i,1)*omegmix3(i,1)
120     continue
       do 124 i=1,nsplit4
        if(iray.eq.15)tau4R(i)=0.
        taumix4(i,1)=tau4R(i)+tau4G(i)
        omegmix4(i,1)=tau4R(i)/taumix4(i,1)
        write(6,*)'Effective Rayliegh SSA',omegmix4(i,1)
        tau_Ray=tau_Ray+tau4R(i)
        tau4(i)=tau4(i)+taumix4(i,1)
        omegtau4(i)=omegtau4(i)+taumix4(i,1)*omegmix4(i,1)
124     continue
C       **************************************************************
C      Get necessary integers for linear interpolation in height
       do 132 i=1,iheight
         do 128 j=1,izeds-1
            if(heights(i).le.Haero(j+1)
     &.and.heights(i).ge.Haero(j))then
            linints(i)=j
            go to 132
            endif    
128       continue
132     continue
C      **************************************************************                                       C      Now get the LOWTRAN profiles of the extinction coefficient       
C      at 550 nm as a function of Height, First get all the stuff       
C      above the  boundary layer. The j in upperprofs 550 is an         
C      aerosol identifier and fall into 3 groups                        
C                                                                       
C **** Upper Tropospheric Background (2-10km) zero elsewhere            
C      1 for Autumn/Winter 50km visibility                              
C      2 for Autumn/Winter 23km visibility                              
C      3 for Spring/Summer 50km  visibility                             
C      4 for Spring/Summer 23km  visibility                             
C **** Stratospheric Aerosols (10-30km) zero elsewhere                  
C      5 Background Stratospheric  Autumn/Winter                        
C      6 Background Stratospheric  Spring/Summer                        
C      7 Moderate Volcanic Autumn/Winter                                
C      8 Moderate Volcanic Spring/Summer                                
C      9 High Volcanic  Autumn/Winter                                   
C      10 High Volcanic Spring/Summer                                   
C      11 Extreme Volcanic  Autumn/Winter                               
C      12 Extreme Volcanic Spring/Summer                                
C **** Upper Stratosphere/Mesosphere/Thermosphere (30-70km) zero below 
C      13 Normal Upper Atmospheric                                      
C      14 Transition from volcanic to normal                             
C      15 Transition from volcanic to extreme volcanic                   
C      16 Transition from Extreme volcanic to volcanic                   
       Call upperprofiles(upperprofs550)                                
C                                                                       
C      Now get the Boundary Layer profile for the extinction            
C      coefficient at 50km,23km,10km,5km,2km                            
C      ( Although defined for all 34 heights, = 0 above 2km)            
       Call BoundaryHaze(BHaze)                                         
                                                                        
C      Now We have the profile of the extinction coefficient at 550nm   
C      We want the extinction AND absorption coefficents at other       
C      Wavelengths PLUS the single scattering albedo. LOWTRAN stores    
C      the wavelength dependent scaling factors to get the first two    
C      First we get the wavelengths (in microns)                        
       Call WaveL(Waves)                                                
                                                                        
C      Next We get the conversion factors that depend on relative Humidity
C      These are in the arrays AeroHumidext(i,j,k) and AeroHumidAbs(i,j,k)
C      The i denotes the wavelength, j is an aerosol identifier and     
C      k =1, 2, 3, 4 are for 0%, 70%, 80% and 99% relative Humidity     
C      The Aerosol identifiers are                                      
C      1 for rural boundary layer (0-2km) zero elsewhere                
C      2 for urban boundary layer                                       
C      3 for maritime or oceanic boundary layer                                     
C      4 for for upper troposphere                                                         C                                                                       
       Call AeroHumid(AeroHumidext,AeroHumidabs)                        
C                                                                       
C      The extinction and absorption coefficients for the upper         
C      regions (humidity variations not effective)                      
C      AeroDryext/abs(ilambda,k) k= identifier                               
C      1  Background stratospheric  (10-30km) zero elsewhere            
C      2  Aged Volcanic                                                 
C      3  Fresh Volcanic                                                
C      4  space dust   (30-100km)                                       
       Call AeroDry(AeroDryExt,AeroDryAbs)
                             
                                                                        
C      Now for the Phase functions.                                     
C      There are iphasf phase functions in LOWTRAN, these are identified    
C      by an integer array which we shall call phaseid(iwave,26)           
C      The column number identifies a model, the values 1-26 are        
C      1=RURAL     0%RH   2=RURAL    70%RH   3=RURAL    80%RH           
C      4=RURAL    99%RH   5=MARITIME  0%RH   6=MARITIME 70%RH           
C      7=MARITIME 80%RH   8=MARITIME 99%RH   9=URBAN     0%RH           
C     10=URBAN    70%RH  11=URBAN    80%RH  12=URBAN    99%RH           
C     13=OCEANIC   0%RH  14=OCEANIC  70%RH  15=OCEANIC  80%RH           
C     16=OCEANIC  99%RH  17=TROPOSPH  0%RH  18=TROPOSPH 70%RH           
C     19=TROPOSPH 80%RH  20=TROPOSPH 99%RH  21=STRATOSPHERIC            
C     22=AGED VOLCANIC   23=FRESH VOLCANIC  24=RADIATION FOG            
C     25=ADVECTIVE FOG   26=METEORIC DUST                               
C     and the row number identifies the wavelength from the different   
C     wavelength array which we now get.                                
                                                                        
       call WavePH(WavePid)
C      get iw1 and 2 for use in linear interpolation
       do 136 i=1,iwave-1
	  if((WavePid(i).lt.lambda).and.(WavePid(i+1).ge.lambda))then
	     iw1=i
	     iw2=i+1
	  endif
136     continue
C                                                                       
C     and then the array of identifiers                                 
      Call phase_id(id_array) 
C     and then the iphasf LOWTRAN Phase functions
      Call GetFFase(PhaseArray,istop)
C
C     interpolate upper atmosphere extinction profiles at 550nm
C **** Upper Stratosphere/Mesosphere/Thermosphere (30-120km) zero below 
C     *********************************************************************
C                Upper Atmosphere Aerosol Profiles 
C     *********************************************************************
       do 140 i=1,izeds
         upext(i)=upperprofs550(i,itypeu) 
140     continue   
                   
C      spline to finer resolution 
C      These data unsuitable for cubic spline - use linear
C      interpolation                                      
       range1=30*1000
       range2=120*1000                     
       do 144 i=1,iheight
        if((heights(i).lt.range1).or.(heights(i).gt.range2))then
            up(i)=0d0
        else
          up(i)=upext(linints(i))+
     &(heights(i)-Haero(linints(i)))
     &/(Haero(linints(i)+1) - Haero(linints(i)) )
     &*(upext(linints(i)+1)-upext(linints(i)) ) 
        endif 
144    continue   
                                   
C      spline to finer resolution                                       
       TD=1d35   
       range1=30*1000
       range2=120*1000 
C
C      Now we need to interpolate to the correct wavelength
       if(itypeu.eq.13)itype1=1
       if(itypeu.eq.14)itype1=2
       if(itypeu.eq.15)itype1=3   
       if(itypeu.eq.16)itype1=1               
       TD=1d35                      
       do 152 i=1,ilambda
          AeroDryE(i)=AeroDryExt(i,itype1)
          AeroDryA(i)=AeroDryAbs(i,itype1)
152     continue
       call spline(Waves,AeroDryE,ilambda,TD,TD,AE2)  
       call spline(Waves,AeroDryA,ilambda,TD,TD,AA2) 
       call splint(Waves,AeroDryE,AE2,ilambda,lambda,exfac) 
       call splint(Waves,AeroDryA,AA2,ilambda,lambda,abfac)
        
       do 160 i=1,nsplit1
         tau1A(i)=0
         do 156 j=ix2+(i-1)*(iorder1-1),ix2+i*(iorder1-1)
           tau1A(i)=
     &tau1A(i)+zstep1*weight1(j,i)*up(j)*exfac/1000.0
           Hum1(i)=Hum1(i)+zstep1*weight1(j,i)*Humids(j)
156      continue
C        tau1A(i)=0
         Hum1(i)=Hum1(i)/thik1
         if(iray.eq.14)tau1A(i)=0
         taumix1(i,2)=tau1A(i)
         tau_Aero=tau_Aero+tau1A(i)
         tau1(i)=tau1(i)+taumix1(i,2)
         ssalb1(i)=(exfac-abfac)/exfac
         omegmix1(i,2)=ssalb1(i)
         omegtau1(i)=omegtau1(i)+taumix1(i,2)*omegmix1(i,2)
160     continue
C     *********************************************************************
C                 Stratospheric Profiles 
C     *********************************************************************
       do 164 i=1,izeds
         strext(i)=upperprofs550(i,itypes) 
164    continue 
C       Get phase function model identifier
       if(itypes.eq.5.or.itypes.eq.6)istratphase=21
       if(itypes.eq.7.or.itypes.eq.8)istratphase=22
       if(itypes.eq.9.or.itypes.eq.10)istratphase=23
       if(itypes.eq.11.or.itypes.eq.12)istratphase=23
C
                                                                         
C      spline to finer resolution                                       
       TD=1d35 
       range1=10*1000
       range2=30*1000                     
       call spline(Haero,strext,izeds,TD,TD,str2)  
       do 170 i=1,iheight                                                
          if((heights(i).lt.range1).or.(heights(i).gt.range2))then
            str(i)=0d0
          else
          call splint(Haero,strext,str2,izeds,heights(i),str(i)) 
          endif  
170     continue   
C
C      Now we need to interpolate to the correct wavelength
       if((itypes.eq.5).or.(itypes.eq.6))itype2=1
       if((itypes.eq.7).or.(itypes.eq.8))itype2=2
       if((itypes.eq.9).or.(itypes.eq.10))itype2=3
       if((itypes.eq.11).or.(itypes.eq.12))itype2=3                  
       TD=1d35                      
       do 174 i=1,ilambda
          AerodryE(i)=AeroDryExt(i,itype2)
          AerodryA(i)=AeroDryAbs(i,itype2)
174     continue
       call spline(Waves,AeroDryE,ilambda,TD,TD,AE2)  
       call spline(Waves,AeroDryA,ilambda,TD,TD,AA2) 
       call splint(Waves,AeroDryE,AE2,ilambda,lambda,exfac) 
       call splint(Waves,AeroDryA,AA2,ilambda,lambda,abfac) 
C
       do 180 i=1,nsplit2
         tau2A(i)=0
         do 176 j=ix3+(i-1)*(iorder2-1),ix3+i*(iorder2-1)
           tau2A(i)=tau2A(i)+zstep2*weight2(j,i)*str(j)*exfac/1000.
           Hum2(i)=Hum2(i)+zstep2*weight2(j,i)*Humids(j)
176      continue
C        tau2A(i)=0

        Hum2(i)=Hum2(i)/thik2
        if(iray.eq.14)tau2A(i)=0
        taumix2(i,2)=tau2A(i)
        tau_Aero=tau_Aero+tau2A(i)
        tau2(i)=tau2(i)+taumix2(i,2)
        ssalb2(i)=(exfac-abfac)/exfac
        omegmix2(i,2)=ssalb2(i)
        omegtau2(i)=omegtau2(i)+taumix2(i,2)*omegmix2(i,2)
180     continue
C     *********************************************************************
C                 Upper Tropospheric Profiles 
C     ********************************************************************* 
C **** Upper Tropospheric Background (2-10km) zero elsewhere            
C **** Upper Tropospheric Background (2-10km) zero elsewhere  
C      tropospheric phase function id
       itypetro=16
       itropphase=itypetro-1
C     *********************************************************************
C      Do Lowtran interpolation for visibility
       if(vist.lt.23.)then
       write(6,*)'visibility too low'
       stop
       endif
       if(vist.gt.50.)then
       write(6,*)'visibility too high'
       stop
       endif

187    const=1./(1./23.-1/50.) 
       do 188 i=1,izeds
         troext(i)=const*(
     & (upperprofs550(i,itypet+1)-upperprofs550(i,itypet))/vist
     & +upperprofs550(i,itypet)/23. -upperprofs550(i,itypet+1)/50.)
188     continue                                                         
C     *********************************************************************                                                C      spline to finer resolution                                       
       TD=1d35        
       range1=2*1000
       range2=10*1000             
       call spline(Haero,troext,izeds,TD,TD,tr2)  
       do 196 i=1
     &,iheight                                                
          if((heights(i).lt.range1).or.(heights(i).gt.range2))then
            tro(i)=0d0
          else
          call splint(Haero,troext,tr2,izeds,heights(i),tro(i)) 
          endif     
196    continue 
C
C      Now find extinction coefficients as function of wavelength
C      and relative humidity.
C      select tropospheric
       itypetr=4
C     *********************************************************************            
C      LINEAR INTERPOLATION FOR HUMIDITY                
       do 198 i=1,nsplit3
       do 197 k=ix4+(i-1)*(iorder3-1),ix4+i*(iorder3-1)
         Hum3(i)=Hum3(i)+zstep3*weight3(k,i)*Humids(k)
197    continue
       Hum3(i)=Hum3(i)/thik3
198    continue
       do 220 i=1,nsplit3
       tau3A(i)=0.
C      was taking integral properly, but mid point humidity
C      is being used to represent the sub layer
       if(Hum3(i).le.70)then
       do 204 j=1,ilambda
          AeroHumE(j)=AeroHumidExt(j,itypetr,1)
     &+(AeroHumidExt(j,itypetr,2)-AeroHumidExt(j,itypetr,1))
     &/70d0*(Hum3(i)-0d0)
          AeroHumA(j)=AeroHumidAbs(j,itypetr,1)
     &+(AeroHumidAbs(j,itypetr,2)-AeroHumidAbs(j,itypetr,1))
     &/70d0*(Hum3(i)-0d0)
204     continue
       endif
       if((Hum3(i).gt.70).and.(Hum3(i).le.80))then
        do 208 j=1,ilambda
          AeroHumE(j)=AeroHumidExt(j,itypetr,2)
     &+(AeroHumidExt(i,itypetr,3)-AeroHumidExt(j,itypetr,2))
     &/10d0*(Hum3(i)-70d0)
          AeroHumA(j)=AeroHumidAbs(j,itypetr,2)
     &+(AeroHumidAbs(j,itypetr,3)-AeroHumidAbs(j,itypetr,2))
     &/10d0*(Hum3(i)-70d0)
208     continue
       endif
       if((Hum3(i).gt.80).and.(Hum3(i).le.100))then
       do 212 j=1,ilambda
          AeroHumE(j)=AeroHumidExt(j,itypetr,3)
     &+(AeroHumidExt(j,itypetr,4)-AeroHumidExt(j,itypetr,3))
     &/20d0*(Hum3(i)-80d0)
          AeroHumA(j)=AeroHumidAbs(j,itypetr,3)
     &+(AeroHumidAbs(j,itypetr,4)-AeroHumidAbs(j,itypetr,3))
     &/20d0*(Hum3(i)-80d0)
212     continue
       endif
       call spline(Waves,AeroHumE,ilambda,TD,TD,AE2)  
       call spline(Waves,AeroHumA,ilambda,TD,TD,AA2) 
       call splint(Waves,AeroHumE,AE2,ilambda,lambda,exfac) 
       call splint(Waves,AeroHumA,AA2,ilambda,lambda,abfac)
       do 216 k=ix4+(i-1)*(iorder3-1),ix4+i*(iorder3-1)
       tau3A(i)=tau3A(i)+zstep3*weight3(k,i)*tro(k)*exfac/1000.0
       ssalb3(i)=ssalb3(i)+zstep3*weight3(k,i)*(exfac-abfac)/exfac
216    continue  
C      tau3A(i)=0.
220    continue  
       do 224 i=1,nsplit3
          ssalb3(i)=ssalb3(i)/thik3
          omegmix3(i,2)=ssalb3(i)
          if(iray.eq.14)tau3A(i)=0
          taumix3(i,2)=tau3A(i)
          tau_Aero=tau_Aero+tau3A(i)
          tau3(i)=tau3(i)+taumix3(i,2)
          omegtau3(i)=omegtau3(i)+taumix3(i,2)*omegmix3(i,2)
224     continue
C     *********************************************************************
C                 Boundary Layer Profiles 
C     ********************************************************************* 



C      Do Lowtran Interpolation between visibilities 
       if(vist.gt.50.)then
       write(6,*)'visibility too high'
       stop
       endif

       do 230 j=1,5  
       if(vist.gt.vis(j))go to 235
230    continue
       write(6,*)'boundary visibility too low'
       stop
235    const=1./(1./vis(j)-1/vis(j-1))
       do  240 i=1,izeds
         boundex(i)=const*(
     & (BHaze(i,j)-BHaze(i,j-1))/visb
     & +BHaze(i,j-1)/vis(j) -BHaze(i,j)/vis(j-1))
240    continue      
                                                                                                                       C      spline to finer resolution                                       
       TD=1d35       
       range1=0
       range2=2*1000              
       call spline(Haero,boundex,izeds,TD,TD,b2)  
       do 248 i=1,iheight                                                
          if((heights(i).lt.range1).or.(heights(i).gt.range2))then
            bound(i)=0d0
          else
          call splint(Haero,boundex,tr2,izeds,heights(i),bound(i))
          endif     
248    continue 
C
       
C      for phase function information 
       if(itypeb.eq.1)iboundphase=0
       if(itypeb.eq.2)iboundphase=8
       if(itypeb.eq.3)then
          if(iocean.eq.1)iboundphase=12
          if(iocean.eq.2)iboundphase=4 
       endif      
C      LINEAR INTERPOLATION FOR HUMIDITY
       do 250 i=1,nsplit4
         Hum4(i)=0.
         do 249 k=1+(i-1)*(iorder4-1),1+i*(iorder4-1)
         Hum4(i)=Hum4(i)+zstep4*weight4(k,i)*Humids(k)
249      continue
         Hum4(i)=Hum4(i)/thik4
250    continue
       do 268 i=1,nsplit4
       tau4A(i)=0.
       if(Hum4(i).le.70)then
       do 252 j=1,ilambda
          AeroHumE(j)=AeroHumidExt(j,itypeb,1)
     &+(AeroHumidExt(j,itypeb,2)-AeroHumidExt(j,itypeb,1))
     &/70d0*(Hum4(i)-0d0)
          AeroHumA(j)=AeroHumidAbs(j,itypeb,1)
     &+(AeroHumidAbs(j,itypeb,2)-AeroHumidAbs(j,itypeb,1))
     &/70d0*(Hum4(i)-0d0)
252     continue
       endif
       if((Hum4(i).gt.70).and.(Hum4(i).le.80))then
        do 256 j=1,ilambda
          AeroHumE(j)=AeroHumidExt(j,itypeb,2)
     &+(AeroHumidExt(j,itypeb,3)-AeroHumidExt(j,itypeb,2))
     &/10d0*(Hum4(i)-70d0)
          AeroHumA(j)=AeroHumidAbs(j,itypeb,2)
     &+(AeroHumidAbs(j,itypeb,3)-AeroHumidAbs(j,itypeb,2))
     &/10d0*(Hum4(i)-70d0)
256     continue
       endif
       if((Hum4(i).gt.80).and.(Hum4(i).le.100))then
       do 260 j=1,ilambda
          AeroHumE(j)=AeroHumidExt(j,itypeb,3)
     &+(AeroHumidExt(j,itypeb,4)-AeroHumidExt(j,itypeb,3))
     &/20d0*(Hum4(i)-80d0)
          AeroHumA(j)=AeroHumidAbs(j,itypeb,3)
     &+(AeroHumidAbs(j,itypeb,4)-AeroHumidAbs(j,itypeb,3))
     &/20d0*(Hum4(i)-80d0)
260     continue
       endif
       call spline(Waves,AeroHumE,ilambda,TD,TD,AE2)  
       call spline(Waves,AeroHumA,ilambda,TD,TD,AA2) 
       call splint(Waves,AeroHumE,AE2,ilambda,lambda,exfac) 
       call splint(Waves,AeroHumA,AA2,ilambda,lambda,abfac)   
       do 264 k=1+(i-1)*(iorder4-1),1+i*(iorder4-1)
          tau4A(i)=
     &tau4A(i)+zstep4*weight4(k,i)*bound(k)*exfac/1000.
         ssalb4(i)=ssalb4(i)+zstep4*weight4(k,i)*(exfac-abfac)/exfac
264     continue
268    continue 
       do 271 i=1,nsplit4
         if(iray.eq.14)tau4A(i)=0
         taumix4(i,2)=tau4A(i)
         tau_Aero=tau_Aero+tau4A(i)
         tau4(i)=tau4(i)+taumix4(i,2)
         ssalb4(i)=ssalb4(i)/thik4
         omegmix4(i,2)=ssalb4(i)
         omegtau4(i)=omegtau4(i)+taumix4(i,2)*omegmix4(i,2)
271    continue
        write(95,*)real(tau_Ray),"	Rayleigh"
        write(95,*)real(tau_Aero),"	   Aerosol"
        close(95,status="keep")
       write(6,*)"wavelength =",lambda
       write(6,*)"Upper Atmosphere Rayleigh Optical Depths"
       write(6,*)(tau1R(i),i=1,nsplit1)
       write(6,*)"Stratosphere Rayleigh Optical Depths"
       write(6,*)(tau2R(i),i=1,nsplit2)
       write(6,*)"Upper Troposphere Rayleigh Optical Depths"
       write(6,*)(tau3R(i),i=1,nsplit3)
       write(6,*)"Boundary Layer Rayleigh Optical Depths"
       write(6,*)(tau4R(i),i=1,nsplit4)
       write(6,*)"total Rayleigh optical depth = ", tau_Ray             
C     &,"   IQBAL =", .008735*lambda**(-4.08)
       write(6,*)"upper atmosphere Aerosol Optical Depths"
       write(6,*)(tau1A(i),i=1,nsplit1)
       write(6,*)"Stratospheric Aerosol Optical Depths"
       write(6,*)(tau2A(I),i=1,nsplit1)
       write(6,*)"Tropospheric Aerosol Optical Depths"
       write(6,*)(tau3A(I),i=1,nsplit3)
       write(6,*)"Boundary layer Optical Depths" 
       write(6,*)(tau4A(i),i=1,nsplit4)
       write(6,*)"Boundary layer Gasseous absorption"
       write(6,*)(tau4G(i),i=1,nsplit4)
       write(6,*)"of which the Water vapour gives"
       write(6,*)(tau4Gas(i,1),i=1,nsplit4)
C    ****************************************************
C        INTERPOLATE FOR PHASE FUNCTIONS
C    **************************************************** 
C     Humidity independent phase functions

       do 300 i=1,imax
	   moments1(i,1)=
     &   (WavePid(iw2)-lambda)/(WavePid(iw2)-WavePid(iw1))
     &*PhaseArray(i,id_array(iw1,iupphase))
     &    +(lambda-WavePid(iw1))/(WavePid(iw2)-WavePid(iw1))
     &*PhaseArray(i,id_array(iw2,iupphase))
300    continue
       do 320 i=1,imax
	   moments2(i,1)=
     &   (WavePid(iw2)-lambda)/(WavePid(iw2)-WavePid(iw1))
     &*PhaseArray(i,id_array(iw1,istratphase))
     &    +(lambda-WavePid(iw1))/(WavePid(iw2)-WavePid(iw1))
     &*PhaseArray(i,id_array(iw2,istratphase))
320    continue

C      Now the humidity dependent ones
       do 350 k=1,ihum
       do 340 i=1,imax
	   moments3a(i,k)=
     &   (WavePid(iw2)-lambda)/(WavePid(iw2)-WavePid(iw1))
     &*PhaseArray(i,id_array(iw1,itropphase+k))
     &    +(lambda-WavePid(iw1))/(WavePid(iw2)-WavePid(iw1))
     &*PhaseArray(i,id_array(iw2,itropphase+k))
340    continue
350    continue
 
       do 380 k=1,ihum
       do 370 i=1,imax
	   moments4a(i,k)=
     &   (WavePid(iw2)-lambda)/(WavePid(iw2)-WavePid(iw1))
     &*PhaseArray(i,id_array(iw1,iboundphase+k))
     &    +(lambda-WavePid(iw1))/(WavePid(iw2)-WavePid(iw1))
     &*PhaseArray(i,id_array(iw2,iboundphase+k))
370    continue
380    continue

       do 440 k=1,nsplit3
       if(Hum3(k).le.70)then
       do 410 i=1,imax
	  moments3(i,k,1)=moments3a(i,1)
     &+Hum3(k)/70.*(moments3a(i,2)-moments3a(i,1))
410    continue
       endif
       if((Hum3(k).gt.70).and.(Hum3(k).le.80))then
        do 420 i=1,imax
	  moments3(i,k,1)=moments3a(i,2)
     &+(Hum3(k)-70.)/10.*(moments3a(i,3)-moments3a(i,2))
420     continue
       endif
       if((Hum3(k).gt.80).and.(Hum3(k).le.100))then
       do 430 i=1,imax
	  moments3(i,k,1)=moments3a(i,3)
     &+(Hum3(k)-80.)/20.*(moments3a(i,4)-moments3a(i,3))
430     continue
       endif
440    continue

       do 490 k=1,nsplit4
       if(Hum4(k).le.70)then
       do 450 i=1,imax
	  moments4(i,k,1)=moments4a(i,1)
     &+Hum4(k)/70.*(moments4a(i,2)-moments4a(i,1))
450    continue
       endif
       if((Hum4(k).gt.70).and.(Hum4(k).le.80))then
        do 460 i=1,imax
	  moments4(i,k,1)=moments4a(i,2)
     &+(Hum4(k)-70.)/10.*(moments4a(i,3)-moments4a(i,2))
460     continue
       endif
       if((Hum4(k).gt.80).and.(Hum4(k).le.100))then
       do 470 i=1,imax
	  moments4(i,k,1)=moments4a(i,3)
     &+(Hum4(k)-80.)/20.*(moments4a(i,4)-moments4a(i,3))
470     continue
       endif
490    continue
       do 498 i=1,iangle
         Angles(i)=pi*(i-1)/(iangle-1)
         xx(i)=cos(Angles(i))
         call lpolcalc(1,istop,legpol,xx(i))
         do 494 j=1,istop
          Legs(j,i)=legpol(j)
494      continue
498     continue

 

C      write(6,*)
C     &"total Optical Depths including gasseous absorption"
      do 817 i=1,nsplit1
        omeg01(i)=0d0
        do 816 j=1,nmix
          omeg01(i)=omeg01(i)+taumix1(i,j)*omegmix1(i,j)/tau1(i)
816       continue
817     continue
       do 819 i=1,nsplit2
        omeg02(i)=0d0
        do 818 j=1,nmix
          omeg02(i)=omeg02(i)+taumix2(i,j)*omegmix2(i,j)/tau2(i)
818       continue
819     continue
       do 821 i=1,nsplit3
        omeg03(i)=0d0
        do 820 j=1,nmix
          omeg03(i)=omeg03(i)+taumix3(i,j)*omegmix3(i,j)/tau3(i)
820     continue
821     continue
        do 823 i=1,nsplit4
        omeg04(i)=0d0
        do 822 j=1,nmix
          omeg04(i)=omeg04(i)+taumix4(i,j)*omegmix4(i,j)/tau4(i)
822     continue
823     continue

      open(21,file="PHASEDATA",status="unknown")
      do 1000 n=1,nsplit1
        nn=nsplit1+1-n
        do 880 i=1,imax
         MOMS1(i,nn)=0d0
880     continue
        do 900 j=1,2
           weighting=omegmix1(nn,j)*taumix1(nn,j)/omegtau1(nn)
              if(j.eq.1)then
              MOMS1(1,nn)=MOMS1(1,nn)+weighting/4d0/pi
              MOMS1(3,nn)=MOMS1(3,nn)+weighting/2d0/4d0/pi
     &*(1.-delta)/(1.+delta/2.)
              endif
              if(j.eq.2)then
              do 890 k=1,imax
                MOMS1(k,nn)=MOMS1(k,nn)+weighting*moments1(k,j-1)
890           continue
              endif
900        continue

           write(21,*)omeg01(nn),tau1(nn),Hum1(nn)
           write(21,*)(MOMS1(k,nn)/MOMS1(1,nn),k=1,imax)
C
       if(nn.lt.10)write(filephase,8881),nn
       if(nn.ge.10)write(filephase,8882),nn
8881   format('upphasedata',i1,'.dat')
8882   format('upphasedata',i2,'.dat')
      open(55,file=filephase,status="unknown")
      open(55,file=filephase,status="unknown")
      eff=MOMS1(imax,nn)/DBLE(2*(istop+1)+1)/MOMS1(1,nn)
      nn=1
      do 506 i=1,iangle
         out=0
         do 505 j=1,istop
           out=out+
     &(MOMS1(j,nn)/MOMS1(1,nn)-DBLE((2*J-1))*eff)*Legs(j,i)/(1D0-eff)
505      continue
         write(55,*)real(Angles(i)*180./pi),real(log10(out))
506   continue
      close(55,status='keep')

1000       continue
      do 2000 n=1,nsplit2
        nn=nsplit2+1-n
        do 1880 i=1,imax
         MOMS2(i,nn)=0d0
1880     continue
        do 1900 j=1,2
           weighting=omegmix2(nn,j)*taumix2(nn,j)/omegtau2(nn)
              if(j.eq.1)then
              MOMS2(1,nn)=MOMS2(1,nn)+weighting/4d0/pi
              MOMS2(3,nn)=MOMS2(3,nn)+weighting/2d0/4d0/pi
     &*(1.-delta)/(1.+delta/2.)
              endif
              if(j.eq.2)then
              do 1800 k=1,imax
                MOMS2(k,nn)=MOMS2(k,nn)+weighting*moments2(k,j-1)
1800           continue
              endif
1900        continue
           write(21,*)omeg02(nn),tau2(nn),Hum2(nn)
           write(21,*)(MOMS2(k,nn)/MOMS2(1,nn),k=1,imax)

       if(nn.lt.10)write(filephase,8883),nn
       if(nn.ge.10)write(filephase,8884),nn
8883   format('stratphasedata',i1,'.dat')
8884   format('stratphasedata',i2,'.dat')
      open(55,file=filephase,status="unknown")
      eff=MOMS2(imax,nn)/DBLE(2*(istop+1)+1)/MOMS2(1,nn)
      do 606 i=1,iangle
         out=0
         do 605 j=1,istop
           out=out+
     &(MOMS2(j,nn)/MOMS2(1,nn)-DBLE((2*J-1))*eff)*Legs(j,i)/(1D0-eff)
605      continue
         write(55,*)real(Angles(i)*180./pi),real(log10(out))
606   continue
      close(55,status='keep')


2000       continue
      do 3000 n=1,nsplit3
         nn=nsplit3+1-n
        do 2880 i=1,imax
         MOMS3(i,nn)=0d0
2880     continue
        do 2900 j=1,nmix
           weighting=omegmix3(nn,j)*taumix3(nn,j)/omegtau3(nn)
              if(j.eq.1)then
              MOMS3(1,nn)=MOMS3(1,nn)+weighting/4d0/pi
              MOMS3(3,nn)=MOMS3(3,nn)+weighting/2d0/4d0/pi
     &*(1.-delta)/(1.+delta/2.)
              endif
              if(j.ge.2)then
              do 2800 k=1,imax
                MOMS3(k,nn)=MOMS3(k,nn)+weighting*moments3(k,nn,j-1)
C                write(6,*)moments3(k,nn,j-1),k-1,nmix,nn
2800           continue
              endif
2900        continue
           write(21,*)omeg03(nn),tau3(nn),Hum3(nn)
           write(21,*)(MOMS3(k,nn)/MOMS3(1,nn),k=1,imax)
       if(nn.lt.10)write(filephase,8885),nn
       if(nn.ge.10)write(filephase,8886),nn
8885   format('tropphasedata',i1,'.dat')
8886   format('tropphasedata',i2,'.dat')
      open(55,file=filephase,status="unknown")
      eff=MOMS3(IMAX,nn)/DBLE(2*(istop+1)+1)/MOMS3(1,nn)
      do 707 k=1,nsplit3 
      do 706 i=1,iangle
         out=0
         do 705 j=1,istop
           out=out+
     &(MOMS3(j,nn)/MOMS3(j,nn)-DBLE((2*J-1))*eff)*Legs(j,i)/(1D0-eff)
705      continue
         write(55,*)real(Angles(i)*180./pi),real(log10(out))
706   continue
707   continue
      close(55,status='keep')

3000       continue
C
      do 4000 n=1,nsplit4
        nn=nsplit4+1-n
        do 3880 i=1,imax
         MOMS4(i,nn)=0d0
3880     continue
        do 3900 j=1,nmix
           weighting=omegmix4(nn,j)*taumix4(nn,j)/omegtau4(nn)
              if(j.eq.1)then
              MOMS4(1,nn)=MOMS4(1,nn)+weighting/4d0/pi
              MOMS4(3,nn)=MOMS4(3,nn)+weighting/2d0/4./pi
     &*(1.-delta)/(1.+delta/2.)
              endif
              if(j.ge.2)then
              do 3800 k=1,imax
                MOMS4(k,nn)=MOMS4(k,nn)+weighting*moments4(k,nn,j-1)
3800           continue
              endif
3900        continue
           write(21,*)omeg04(nn),tau4(nn),Hum4(nn)
           write(21,*)(MOMS4(k,nn)/MOMS4(1,nn),k=1,imax)


       if(nn.lt.10)write(filephase,8887),nn
       if(nn.ge.10)write(filephase,8888),nn
8887   format('boundphasedata',i1,'.dat')
8888   format('boundphasedata',i2,'.dat')
      open(55,file=filephase,status="unknown")
      eff=MOMS4(IMAX,nn)/DBLE(2*(istop+1)+1)/MOMS4(1,nn)
      do 807 k=1,nsplit4 
      do 806 i=1,iangle
         out=0
         do 805 j=1,istop
           out=out+
     &(MOMS4(j,nn)/MOMS4(1,nn)-DBLE((2*J-1))*eff)*Legs(j,i)/(1D0-eff)
805      continue
         write(55,*)real(Angles(i)/pi*180.),real(log10(out))
806   continue
807   continue
      close(55,status='keep')


4000       continue





       write(6,*)"total aerosol optical depth=",tau_Aero
       write(6,*)"total Rayleigh optical depth=",tau_Ray
       close(21,status="keep")  
       open(21,file='tmptau')
       write(21,*)tau_Aero
       close(21,status="keep")  
       open(21,file="tmptransmit")
       do 4001 nn=1,nsplit4
         total4=total4+tau4(nn)-tau4R(nn)-tau4A(nn)
4001   continue
       do 4002 nn=1,nsplit4
         total3=total3+tau3(nn)-tau3R(nn)-tau3A(nn)
4002   continue
       do 4003 nn=1,nsplit4
         total2=total2+tau2(nn)-tau2R(nn)-tau2A(nn)
4003   continue
       do 4004 nn=1,nsplit4
         total1=total1+tau1(nn)-tau1R(nn)-tau1A(nn)
4004   continue
       write(6,*)
     &"GAS TRANSMITTANCE =",exp(-total4-total3-total2-total1),lambda
       write(21,*)exp(-total4-total3-total2-total1)
       close(21,status="keep")

       if(iplane.eq.1)then
       write(6,*)"Aeroplane Data"
       iplanesplit1=nsplit1
       iplanesplit2=nsplit2
       iplanesplit3=nsplit3
       iplanesplit4=nsplit4

       do 5000 i=1,nsplit1
         if(heightplane.ge.heights(ix2+(i-1)*(iorder1-1) ).and.
     &   heightplane.lt.heights(ix2+i*(iorder1-1)) )then
         iplanelayer=1
         iplanesplit1=i
       endif
5000   continue
       do 5001 i=1,nsplit2
         if(heightplane.ge.heights(ix3+(i-1)*(iorder2-1)).and.
     &   heightplane.lt.heights(ix3+i*(iorder2-1)) )then
         iplanelayer=2
         iplanesplit2=i
         endif
5001   continue
       do 5002 i=1,nsplit3
         if(heightplane.ge.heights(ix4+(i-1)*(iorder3-1)).and.
     &   heightplane.lt.heights(ix4+i*(iorder3-1)) )then
         iplanelayer=3
         iplanesplit3=i
       endif
5002   continue
       do 5003 i=1,nsplit4
         if(heightplane.ge.heights(1+(i-1)*(iorder4-1)).and.
     &   heightplane.lt.heights(1+iorder4+i*(iorder4-1)) )then
         iplanelayer=4
         iplanesplit4=i
       endif
5003   continue
       write(6,*)"iplanelayer and iplanesplits"
     &,iplanelayer,iplanesplit1,iplanesplit2,iplanesplit3,iplanesplit4


      open(27,file="PHASEPLANE",status="unknown")
      if(iplanelayer.lt.2)then
      istart1=1
      if(iplanelayer.eq.1)then
      istart1=iplanesplit1+1
        n=iplanesplit1
        nn=nsplit2+1-n
           write(27,*)omeg01(nn),
     &(heightplane-heights(ix2+(n-1)*(iorder1-1)))/thik1*tau1(nn),1,n 
           write(27,*)(MOMS1(k,nn)/MOMS1(1,nn),k=1,imax)
      endif
      do 6000 n=istart1,nsplit1
         nn=nsplit1+1-n
           write(27,*)omeg01(nn),tau1(nn),1,n
           write(27,*)(MOMS1(k,nn)/MOMS1(1,nn),k=1,imax)
6000       continue
      endif
      if(iplanelayer.le.2)then
      istart2=1
      if(iplanelayer.eq.2)then
      istart2=iplanesplit2+1
        n=iplanesplit2
        nn=nsplit2+1-n
           write(27,*)omeg02(nn),
     &(heightplane-heights(ix3+(n-1)*(iorder2-1)))/thik2*tau2(nn),2,n 
           write(27,*)(MOMS2(k,nn)/MOMS2(1,nn),k=1,imax)
      endif
      do 6100 n=istart2,nsplit2
        nn=nsplit2+1-n
           write(27,*)omeg02(nn),tau2(nn),2,n
           write(27,*)(MOMS2(k,nn)/MOMS2(1,nn),k=1,imax)
6100       continue
       endif
      if(iplanelayer.le.3)then
      istart3=1
      if(iplanelayer.eq.3)then
      istart3=iplanesplit3+1
        n=iplanesplit3
         nn=nsplit3+1-n
           write(27,*)omeg03(nn),
     &(heightplane-heights(ix4+(n-1)*(iorder3-1)))/thik3*tau3(nn),3,n 
           write(27,*)(MOMS3(k,nn)/MOMS3(1,nn),k=1,imax)
      endif
      do 6200 n=istart3,nsplit3
         nn=nsplit3+1-n
           write(27,*)omeg03(nn),tau3(nn),3,n
           write(27,*)(MOMS3(k,nn)/MOMS3(1,nn),k=1,imax)
6200       continue
      endif
      if(iplanelayer.le.4)then
      istart4=1
      if(iplanelayer.eq.4)then
      istart4=iplanesplit4+1
        n=iplanesplit4
        nn=nsplit4+1-n
           write(27,*)omeg04(nn),
     &(heightplane-heights(1+(n-1)*(iorder4-1)))/thik4*tau4(nn),4,n 
           write(27,*)(MOMS4(k,nn)/MOMS4(1,nn),k=1,imax)
      endif
      do 6300 n=istart4,nsplit4
        nn=nsplit4+1-n
           write(27,*)omeg04(nn),tau4(nn),4,n
           write(27,*)(MOMS4(k,nn)/MOMS4(1,nn),k=1,imax)
6300   continue
       endif
       endif
       open(31,file='ATMDATA',status='old')
       read(31,*)inewat
       read(31,*)ibuild
       read(31,*)isimpsonG
       read(31,*)isimpsonNG
       read(31,*)num
       read(31,*)igauss
       close(31,status='delete')
       ItrunkA=Ifour
       if(iray.eq.14)ItrunkA=3
       nlayer=nsplit1+nsplit2+nsplit3+nsplit4
       open(31,file='ATMDATA',status='new')
       write(31,*)inewat,"       inewat"
       write(31,*)ibuild,"       ibuild"
       write(31,*)isimpsonG,"       isimpsonG"
       write(31,*)isimpsonNG,"       isimpsonNG"
       write(31,*)num,"       num"
       write(31,*)igauss,"       igauss"
       write(31,*)nlayer,"       nlayer"
       write(31,*)ItrunkA,"       ItrunkA"
       close(31,status='keep')
       depth=0.
       do 7000 i=1,nsplit1
         aitch(i)=
     & heights(ix2+(nsplit1+1-i)*(iorder1-1))
         tauheights(i)=depth
         depth=depth+tau1(nsplit1+1-i)
7000   continue
       do 7001 i=1,nsplit2
         aitch(nsplit1+i)=
     & heights(ix3+(nsplit2+1-i)*(iorder2-1))
         tauheights (nsplit1+i)=depth
         depth=depth+tau2(nsplit2+1-i)
7001   continue
       do 7002 i=1,nsplit3
         aitch(nsplit1+nsplit2+i)=
     & heights(ix4+(nsplit3+1-i)*(iorder3-1))
         tauheights(nsplit1+nsplit2+i)=depth
         depth=depth+tau3(nsplit3+1-i)
7002   continue
       do 7003 i=1,nsplit4
         aitch(nsplit1+nsplit2+nsplit3+i)=
     & heights(1+(nsplit4+1-i)*(iorder4-1))
         tauheights(nsplit1+nsplit2+nsplit3+i)=depth
         depth=depth+tau4(nsplit4+1-i)
7003   continue
        aitch(nlayer+1)=0.
        tauheights(nlayer+1)=depth
      open(1,file='tau_heights',status='unknown')

        do 7004 i=1,nlayer+1
          write(1,*)
     &aitch(i)
     &,tauheights(i)
7004    continue
        close(1,status='keep')
       return      
       end  

                                    
