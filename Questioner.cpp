#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

void Questioner(bool& verbose,
               int &nsplit1, int &nsplit2, int &nsplit3, int &nsplit4, bool& vis, double &visb,
               double &vist, double &lambda1, double &lambda2, bool &rfun, int &iatm, bool &switchR, bool  &switchA, int& itypeu,
               int& itypes, int& itypet, int& itypeb, bool& ocean, bool& groundP, bool& groundT,
               int& ihumid, double& groundtemp, double& groundpress,
               bool& default_pause, double& HG, double& HB, double& HT, double& HS, double& HU, int& ngaussL,
               bool& calcspec, bool& outspec, bool& PlotIt, bool& logplot, bool& calcRef, bool& outRef,
               int& ngaussLs, int& icutL, int& icutD, int& istep,
               const char* ReadInput){

vist=-1;vist=-1;lambda1=-1; lambda2=-1; nsplit1=-1;nsplit2=-1;nsplit3=-1;nsplit4=-1;
iatm=-1;switchR=0;switchA=0;itypeu=-1;itypes=-1;itypet=-1;itypeb=-1;ocean=0;
ihumid=-1; groundtemp=-1e6;groundpress=-1; default_pause=true;

  cout << "verbose=" << verbose << " Input file =" << ReadInput << endl;

  bool summer=false;
  bool know=false;

  if(!verbose)
  {
   ifstream fp_in;
   fp_in.open(ReadInput, ios::in);
   if(!fp_in.is_open()){ cerr << "Failed to open file (" << ReadInput << ")" << endl;  exit(1);}
      fp_in >> lambda1 >> lambda2 >> rfun;
      fp_in >> calcspec >> outspec >> PlotIt, logplot, calcRef, outRef;
      fp_in >> ngaussLs;
      fp_in >> icutL >> icutD >> istep;
      fp_in >> vis >>  vist  >> visb;
      fp_in >>  nsplit1 >> nsplit2 >> nsplit3 >> nsplit4;
      fp_in >> iatm;
      fp_in >> switchR >> switchA;
      fp_in >> itypeu >> itypes >> itypet >> itypeb >> ocean;
      fp_in >> ihumid >> groundP >> groundT >> groundtemp >> groundpress;
      fp_in >> default_pause >> HG >> HB >> HT >> HS >> HU;
      fp_in >> ngaussL;
      fp_in.close();
  }
  else  // then verbose
 {


  FILE *fp0;
  fp0=fopen("Dors1in.dat", "w");

  cout << "enter 1 if we are using the instrument response function in file=myresp.dat \n";
  cout << "Enter 0 for user defined range with top hat response\n";
  cin >> rfun;
  if(!rfun){
  cout << "Enter lowest wavelength value\n";
  cin >> lambda1;
  cout << "Enter highest wavelength value\n";
  cin >> lambda2;}
  cout << "Enter 1 if you want to calculate the molecular absorption, 0 if you don't\n";
  cin >> calcspec;
  cout << "Enter 1 if you want to output a fine detailed spectrum, enter 0 otherwise\n";
  cin >> outspec;
  if(outspec){
  cout <<" Enter 1 if you are using the Chris Godsalve's PlotIt program to look at the output - 0 otherwise\n";
  cin >>PlotIt;
  cout << "Enter 1 if you want to plot the log of the spectrum- 0 otherwise\n";
  cin >> logplot;
  }
  else{
    PlotIt=false; logplot=false;
  }
  cout << "Enter 1 if you want to calculate the refractive index via the Kramers-Kronig relations\n";
  cout << "and 0 otherwise (calculating it is computationally expensive!)\n";
  cin >> calcRef;
  cout <<"Enter 1 if you want a to output the fine detailed refractive index spectrum\n";
  cin >> outRef;
  cout << "Enter order of Gauss-Legendre quadrature rule for integration over distributions\n";
  cin >> ngaussLs;
  cout << "Enter (integer) number of linewidths for cut-off for a Lorentz profile\n";
  cin >> icutL;
  cout << "Enter (integer) number of linewidths for cut-off for a Doppler profile\n";
  cin >> icutD;
  cout << "Enter istep (integer that sets the resolution for each line --- make it 5 or more)\n";
  cin >> istep;


  cout << "Enter 1 if you want to use visibility, 0 for optical depths instead\n";
  cin >> vis;

  if(vis){
  cout << "Enter a boundary layer visibility\n";
  cout << "Range must be between 2 and 50km.\n";
  cin >>  visb;
  cout << "Enter a troposphere  visibility\n";
  cout << "Range must be between 23 and 50km.\n";
  cin >>  vist;}
  else{
  cout << "Enter the total boundary layer optical depth at 550 micro-metres\n";
  cin >>  visb;
  cout << "Enter the troposphere  optical depth (exluding bdry layer) \n";
  cin >>  vist;}

  cout << "The upper atmosphere shall be split into of nsplit1 layers\n";
  cout << "enter nsplit1 ...\n";
  cin >> nsplit1;

  cout << "The stratosphere shall be split into nsplit2 layers\n";
  cout << "enter nsplit2\n";
  cin >> nsplit2;
  
  cout << "The upper troposphere shall be split into nsplit3 layers\n";
  cin >> nsplit3;

  cout << "The boundary layer shall be split into nsplit4\n";
  cout << "enter nsplit4\n";
  cin >> nsplit4;

  cout <<  "Atmosphere type?\n";
  int isummer=0;
  cout << "Enter 1 for Spring/Summer\n";
  cout << "Enter 0 for Autumn/Winter\n";
  cin >> isummer;
  if(isummer==1)summer=true;

  cout << "Enter 1 for Tropical \n";
  if(summer){
  cout << "Enter 2 for Mid Latitude Spring/Summer\n ";
  cout << "Enter 4 for Sub-Arctic Spring/Summer \n";
  }
  else
  {
  cout << "Enter 3 for Mid Latitude Autumn/Winter \n";
  cout << "Enter 5 for Sub-Arctic Autumn/Winter \n";
  }
  cout << "Enter 6 for Standard \n";
  cout << "Enter 7 for User-Defined \n";
  cin >> iatm;

  cout <<  "enter 1 to switch off Rayliegh Scattering, 0 otherwise \n";
  cin >> switchR;
  cout << "enter 1 to switch off Aerosol Scattering, 0 otherwise \n";
  cin >> switchA;


  cout << "Upper Atmosphere Type\n";

  cout << "Enter 13 for normal upper atmosphere\n";
  cout << "Enter 14 for transition from volcanic to normal\n";
  cout << "Enter 15 for transition from volcanic to extreme volcanic\n";
  cout << "Enter 16 for extreme volcanic\n";
  cin >> itypeu;

   

  if(summer)
      {
      cout << "Enter 6 for Background Stratospheric Spring Summer\n";
      cout << "Enter 8 for Moderate Volcanic  Spring Summer\n";
      cout << "Enter 10 for High Volcanic Spring Summer\n";
      cout << "Enter 12 for Extreme Volcanic   Spring Summer\n";
      }
     else
      {
      cout << "Enter 5 for Background Stratospheric Autumn/Winter\n";
      cout << "Enter 7 for Moderate Volcanic Autumn/Winter\n";
      cout << "Enter 9 for High Volcanic Autumn/Winter\n";
      cout << "Enter 11 for Extreme Volcanic Autumn/Winter\n";
      }
      cin >> itypes;

      if(summer){
          itypet=3;}
      else{
          itypet=1;}

      cout << "Boundary layer type\n";
      cout << "Enter 1 for rural\n";
      cout << "Enter 2 for urban\n";
      cout << "Enter 3 for maritime or ocean\n";
      cin >> itypeb;

     if(itypeb==3){
          cout << "Enter 1 for open  ocean, 0 for maritime\n";
          cin >> ocean;
         }

     cout << "Enter 1 to use humidity profile as per Atmospheres/Mols\n";
     cout << "Enter 2 to use your user defined humidity profile\n";
     cout << "contained in  Atmospheres/Hum1_N.dat and Atmospheres/Hum2_N.dat\n";
     cin >> ihumid;

     cout << "Enter 1 if you have ground pressure and/or temperature data, enter 0 otherwise\n";
     cin >> know;
    if(know){
  //      cin >>  know;
            if(know){
               int iground, ipress;
               cout << "Enter 1 if you know the ground temperature, 0 otherwise\n";
               cin >> iground;
               if(iground==1){
               groundT=true;
               cout << "Enter ground temperature\n";
               cin >> groundtemp;}
               cout << "Enter 1 if you know the ground pressure, 0 otherwise\n";
               cin >> ipress;
               if(ipress==1){
               groundP=true;
               cout << "Enter ground pressure\n";
               cin >> groundpress;}
               }
        }

      cout << "Enter 1 for  default ground,  bdry layer, tropopause, stratopause, TOA heights\n";
      cout << "Enter 2 to enter user defined values for these\n";
      cin >> default_pause;
      if(default_pause){
        HB=2.0; HT=10.0; HS=30.0; HU=120.0;}else
        {
        cout << "All heights refer to sea-level\n";
        cout << "Enter height of ground (in km)\n";
        cout << "(greater than zero unless user defined atmosphere)\n"; 
        cout << "(Default=0km)\n";
        cin >> HG;
        cout << "Enter height of Boundary Layer (in km)\n";
        cout << "(Default=2km)\n";
        cin >> HB;
        cout << "Enter height of Tropopause (in km)\n";
        cout << "(Default=10km)\n";
        cin >> HT;
        cout << "Enter height of Stratopause (in km)\n";
        cout << "(Default=30km)\n";
        cin >> HS;
        cout << "Enter height of Top of Atmosphere (km)\n";
        cout << "(Default=120km)\n";
        cin >> HU;}
        
        cout << "Enter order of Gauss-Legendre quadrature rule for integration over height\n";
        cin >> ngaussL;

        fprintf(fp0, "%lf %lf %d \n", lambda1, lambda2, rfun);
        fprintf(fp0, "%d %d\n", calcspec, outspec);
        fprintf(fp0, "%d\n", ngaussLs);
        fprintf(fp0, "%d %d\n", icutL, icutD, istep);

        fprintf(fp0, "%d %lf %lf\n", vis, vist, visb);
        fprintf(fp0, "%d %d %d %d\n", nsplit1, nsplit2, nsplit3, nsplit4);
        fprintf(fp0, "%d\n", iatm);
        fprintf(fp0, "%d %d\n", switchR, switchA);
        fprintf(fp0, "%d %d %d %d %d\n", itypeu, itypes, itypet, itypeb, ocean);
        fprintf(fp0, "%d %d %d %lf %lf\n", ihumid, groundP, groundT, groundtemp, groundpress);
        fprintf(fp0, "%d %lf %lf %lf %lf %lf\n", default_pause, HG, HB, HT, HS, HU);
        fprintf(fp0, "%d\n", ngaussL);
        fclose(fp0);

     }  //end for else verbose
        

}




