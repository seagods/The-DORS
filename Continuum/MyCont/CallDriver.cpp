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
             double&, double&, double&, double&, double& );

}



int main(){

      double Press, Temp, Xpath, WaterMix;
      double Wave1, Wave2, DWave;
      double CO2Mix, NitroMix, OxyMix, OzoMix;

      double KNTH2O[5050], KNTCO2[5050], KNTO3[5050];
      double KNTO2[5050], KNTN2[5050], KNTRay[5050];

      int NumPts;

//     Pressure in mbar
      Press=1012.0;
//     Temperature in Kelvin
      Temp=290.0;
//     Path Lengh in cm
      Xpath=1.0;



//     Lowest wave number
      Wave1=1.0;
      Wave2=10000.0;
//     step in wave number
      DWave=2.0;
//     Number of Data points
//     NPTABS =  1 + (V2abs-V1abs)/dvabs
      NumPts=1+(int)((Wave2-Wave1)/DWave);
      if(NumPts<=0){"Number of Points wrong in Continuum \n"; exit(1);}
      if(NumPts>5040){"Number of Points too large in Continuum \n"; exit(1);}

//     WaterMix=Volume mixing ratio of H2O
      WaterMix=0.01;
//     in ppmv
      CO2Mix=370.0;
      OzoMix=10.0;
//     convert to fraction
      CO2Mix=CO2Mix*1E-6; 
      OzoMix=OzoMix*1E-6;
//     plain fraction to start with
      NitroMix= 0.78;
      OxyMix= 0.21;

      driver_(KNTH2O, KNTCO2, KNTO3, KNTO2, KNTN2, KNTRay
             , Temp,Press,Xpath, Wave1, Wave2,DWave
             , WaterMix, CO2Mix, NitroMix, OxyMix, OzoMix);

      ofstream OutPut;

        OutPut.open("data.dat",ios::out);

        OutPut << "1\n";
        OutPut << NumPts << "  0 1 0 1\n";

        double WV;
        double logK;

        cout << "Number of points=" << NumPts << endl;

        int kountpos=0;

        for(int i=0; i<NumPts; i++){
           WV=Wave1+i*DWave;
           if(KNTH2O[i]>0){
             logK=log(KNTH2O[i]);
             OutPut << WV << "  " << logK << endl;
             kountpos++;
           }
       }
        cout << "Number of points=" << kountpos << endl;

  
       OutPut << 1 << "  " <<  1 << endl;
       OutPut << "Wave@Number   @per_cm" << endl;

       OutPut << 1 << "  " <<  0 << endl;  //for vertical axis
       OutPut << "log10@tau" << endl;

}


