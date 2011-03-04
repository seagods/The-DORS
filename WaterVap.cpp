#include <math.h>
#include <stdio.h>
void WaterVap(double* Press, double* Temp, double* Gas,
              double* RelHum, double &AV, double &nucleon, double &MWH2O, 
              int &iheight, int &ihumid,int &iwaterice){

         double C1W,C2W,C3W,C4W,C5W;
         double C1I,C2I,C3I,C4I,C5I;
         C1W=1.0007;C2W=3.46E-6;C3W=6.1121; C4W=17.502;C5W=240.97;
         C1I=1.0003;C2I=4.18E-6;C3I=6.1115; C4I=22.452;C5I=272.55;
         double PZero=1013.25;
         double TZero=273.15;
         double Losch=2.6868E25;   //per cubic metre not per cc
         double ENNAIR,ENNH2O;
         double DensH, DensA, mixrat,satmixrat;
         double MWAIR=28.946;

         double ratP,ratT,arg,part,partsat,TC,enhance,exparg;

/* *************************************************************  

         PV=n R T, if m= molec weight and M is mass and rho is density
         PV/M=P/rho=mRT     
         moist air

         P/rho_air =m_air RT,  P/rho_vap=m_vap RT

         given mixing ratio in Kg per Kg

         mix=rho_vap/rho_air=(P_vap m_vap)/(P_air m_air)    ---RTs cancel
         so P_vap=mix m_air/m P_air

         Data and formula below are from Arden L. Buck
         Journal of Applied Met Vol20 pp 1527-1532 Dec 1981


     *************************************************************  */
         if(ihumid==1){
             for(int i=0; i<iheight; i++){
               ratP=Press[i]/PZero;
               ratT=TZero/Temp[i];
               TC=Temp[i]-TZero;
               ENNAIR=Losch*ratP*ratT;
               ENNH2O=ENNAIR*Gas[i]/1E6;  //Since ppm

               DensA=ENNAIR*MWAIR*nucleon/1000.0;
               DensH=ENNH2O*MWH2O*nucleon/1000.0;

               //   partial pressure in mbar
               mixrat=DensH/DensA;
               part=mixrat*MWAIR/MWH2O*Press[i];
               if(iwaterice==1){
                 enhance=C1W+C2W*Press[i];
                 arg= C4W*TC/(C5W+TC);
                 exparg=exp(arg);
                 partsat=enhance*C3W*exparg;}
               else{
                 enhance=C1I+C2I*Press[i];
                 arg= C4I*TC/(C5I+TC);
                 exparg=exp(arg);
                 partsat=enhance*C3I*exparg;}

               satmixrat=partsat/Press[i]*MWH2O/MWAIR;
               RelHum[i]=part/partsat*100.0;

               }
         }
         if(ihumid==2){
             for(int i=0; i<iheight; i++){
               ratP=Press[i]/PZero;
               ratT=TZero/Temp[i];
               TC=Temp[i]-TZero;
               ENNAIR=Losch*ratP*ratT;
               DensA=ENNAIR*MWAIR*nucleon/1000.0;

               if(iwaterice==1){
                 enhance=C1W+C2W*Press[i];
                 arg= C4W*TC/(C5W+TC);
                 exparg=exp(arg);
                 partsat=enhance*C3W*exparg;}
               else{
                 enhance=C1I+C2I*Press[i];
                 arg= C4I*TC/(C5I+TC);
                 exparg=exp(arg);
                 partsat=enhance*C3I*exparg;}

               satmixrat=partsat/Press[i]*MWH2O/MWAIR;
               part=partsat*RelHum[i]/100.;

               mixrat=part*MWH2O/MWAIR/Press[i];
               DensH=DensA*mixrat;
               ENNH2O=DensH/MWH2O*AV*1000.0;
               Gas[i]=ENNH2O/ENNAIR*1E6;
               }
         }
 }    
        
