#include <stdio.h>
#include <iostream>
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision
#include <sstream>  //string stream
#include <string>   // nice string operations
#include <vector>   // a bit like a souped up stack
#include <sys/stat.h>  //POSIX stuff including file queries
#include <new>   // use nothrow
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

using namespace std;

extern "C"{

void driver_( double*, double*, double*, double*, double*, double*,
             double&, double&, double&, double&, double&, double&,
             double&, double&, double&, double&, double&);
}

int main(){

      double Press, Temp, Xpath, WaterMix;
      double Wave1, Wave2, DWave;
      double CO2Mix, NitroMix, OxyMix, OzoMix;

      double KNTH2O[5050], KNTCO2[5050], KNTO3[5050];
      double KNTO2[5050], KNTN2[5050], KNTRay[5050];

      double R,AV;

      int NumPts;

//     Pressure in mbar
      Press=1035.25;
//     Temperature in Kelvin
      Temp=303.0;
//     Path Lengh in cm
      Xpath=1.0;

      R=8.314472;
      AV=6.0221415e23;             //Avogadro's number;

//     1 bar = 10^5 Pascals

       double moles_per_vol;
       moles_per_vol=Press*1.0e2/R/Temp;  //*100 -> Pascals

       double Mols;
       Mols=moles_per_vol*AV*1.0e-6;  //1e-6 convert from per m**3 to per cm**3

       //Note Mols 

       cout << "Number of Molecules per cubic metre=" <<Mols << endl;


//     Lowest wave number
      Wave1=1.0;
      Wave2=2500.0;
      Wave1=1460.0;
      Wave2=1481.0;
      Wave1=8000.0;
      Wave2=16000.0;
//     step in wave number
      DWave=2.0;
//     Number of Data points
//     NPTABS =  1 + (V2abs-V1abs)/dvabs
     
      NumPts=1+(int)((Wave2-Wave1)/DWave);

      cout << NumPts << "  NumPts\n";
      if(NumPts<=0){"Number of Points wrong in Continuum \n"; exit(1);}
      if(NumPts>5040){"Number of Points too large in Continuum \n"; exit(1);}

//     WaterMix=Volume mixing ratio of H2O
       WaterMix=0.0217;  //obviously not ppmv
       double WaterMols=WaterMix*Mols; 
       cout << "Number of Water Molecules per CC=" <<WaterMols << endl;
       double WaterMolSquare;
       // number of molecules per cm^2 in path length Xpath
       WaterMolSquare=WaterMols*Xpath;

       //Output of Driver (KNTH2O) proportional to Xpath
       //Output to file divided by WaterMols (so is per molecule)
       //and Xpath -- so Xpath always 1cm!

       //Daft, but did it just to check KNTH2O really was \propto XPath


//     in ppmv
      CO2Mix=370.0;
      OzoMix=10.0;
//     convert to fraction
      CO2Mix=CO2Mix*1E-6; 
      OzoMix=OzoMix*1E-6;
//     plain fraction to start with
      NitroMix= 0.78;
      OxyMix= 0.21;

//    Now for rest
      double CO2Mols, OzoMols, OxyMols, NitroMols;
      CO2Mols=CO2Mix*Mols; OzoMols=OzoMix*Mols; 
      OxyMols=OxyMix*Mols; NitroMols=NitroMix*Mols;
      

      driver_(KNTH2O, KNTCO2, KNTO3, KNTO2, KNTN2, KNTRay
             , Temp, Press, Xpath, Wave1, Wave2, DWave
             , WaterMix, CO2Mix, NitroMix, OxyMix, OzoMix);

      ofstream OutPut1;
      ofstream OutPut2;
      ofstream OutPut3;
      ofstream OutPut4;
      ofstream OutPut5;


        OutPut1.open("H2Odata.dat",ios::out);
        OutPut2.open("CO2data.dat",ios::out);
        OutPut3.open("OZdata.dat",ios::out);
        OutPut4.open("O2data.dat",ios::out);
        OutPut5.open("N2data.dat",ios::out);

        OutPut1 << "1\n";
        OutPut1 << NumPts << "  0 1 0 1\n";
        OutPut2 << "1\n";
        OutPut2 << NumPts << "  0 1 0 1\n";
        OutPut3 << "1\n";
        OutPut3 << NumPts << "  0 1 0 1\n";
        OutPut4 << "1\n";
        OutPut4 << NumPts << "  0 1 0 1\n";
        OutPut5 << "1\n";
        OutPut5 << NumPts << "  0 1 0 1\n";

        double WV;
        double logK;

        cout << "Number of points=" << NumPts << endl;

        int kountpos[5];
        kountpos[0]=0;kountpos[1]=0;kountpos[2]=0;kountpos[3]=0;
        kountpos[4]=0;

        for(int i=0; i<NumPts; i++){
           WV=Wave1+i*DWave;
           if(KNTH2O[i]>0){
             logK=log10(KNTH2O[i]/WaterMolSquare);
             OutPut1 << WV << "  " << logK << endl;
             kountpos[0]++;
           }
           if(KNTCO2[i]>0){
             logK=log10(KNTCO2[i]/CO2Mols);
             OutPut2 << WV << "  " << logK << endl;
             kountpos[1]++;
           }
           if(KNTO3[i]>0){
             logK=log10(KNTO3[i]/OzoMols);
             OutPut3 << WV << "  " << logK << endl;
             kountpos[2]++;
           }
           if(KNTO2[i]>0){
             logK=log10(KNTO2[i]/OxyMols);
             OutPut4 << WV << "  " << logK << endl;
             kountpos[3]++;
           }
           if(KNTN2[i]>0){
             logK=log10(KNTN2[i]/NitroMols);
             OutPut5 << WV << "  " << logK << endl;
             kountpos[4]++;
           }
       }
        cout << "Number of points H2O=" << kountpos[0] << endl;
        cout << "Number of points CO2=" << kountpos[1] << endl;
        cout << "Number of points O3=" << kountpos[2] << endl;
        cout << "Number of points O2=" << kountpos[3] << endl;
        cout << "Number of points N2=" << kountpos[4] << endl;

  
       OutPut1 << 1 << "  " <<  1 << endl;
       OutPut1 << "Wave@Number   @per_cm" << endl;

       OutPut1 << 1 << "  " <<  0 << endl;  //for vertical axis
       OutPut1 << "log10@tau" << endl;

       OutPut2 << 1 << "  " <<  1 << endl;
       OutPut2 << "Wave@Number   @per_cm" << endl;

       OutPut2 << 1 << "  " <<  0 << endl;  //for vertical axis
       OutPut2 << "log10@tau" << endl;

       OutPut3 << 1 << "  " <<  1 << endl;
       OutPut3 << "Wave@Number   @per_cm" << endl;

       OutPut3 << 1 << "  " <<  0 << endl;  //for vertical axis
       OutPut3 << "log10@tau" << endl;

       OutPut4 << 1 << "  " <<  1 << endl;
       OutPut4 << "Wave@Number   @per_cm" << endl;

       OutPut4 << 1 << "  " <<  0 << endl;  //for vertical axis
       OutPut4 << "log10@tau" << endl;

       OutPut5 << 1 << "  " <<  1 << endl;
       OutPut5 << "Wave@Number   @per_cm" << endl;

       OutPut5 << 1 << "  " <<  0 << endl;  //for vertical axis
       OutPut5 << "log10@tau" << endl;

}


