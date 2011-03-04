#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

void Questioner(int& imode, bool& verbose,
     int &nsplit1, int &nsplit2, int &nsplit3, int &nsplit4, bool& vis, double &visb,
     double &vist, double &lambda1, double &lambda2, int &iatm, bool &switchR, bool  &switchA, int& itypeu,
     int& itypes, int& itypet, int& itypeb, bool& ocean, int& ihumid, bool&  know,
     double& groundtemp, double& groundpress, bool& aeroplane, double& heightplane,
     bool& default_pause, double& HG, double& HB, double& HT, double& HS, double& HU,
     int& ngaussL, int& ngaussH, const char* ReadInput){
imode=1;
vist=-1;vist=-1;lambda1=-1; lambda2=-1; nsplit1=-1;nsplit2=-1;nsplit3=-1;nsplit4=-1;
iatm=-1;switchR=0;switchA=0;itypeu=-1;itypes=-1;itypet=-1;itypeb=-1;ocean=0;
ihumid=-1;know=false;groundtemp=-1e6;groundpress=-1;
aeroplane=0;heightplane=-1e6, default_pause=true;

  cout << "verbose=" << verbose << " Input file =" << ReadInput << endl;

  bool summer=false;

  if(!verbose)
  {
   ifstream fp_in;
   fp_in.open(ReadInput, ios::in);
   if(!fp_in.is_open()){ cerr << "Failed to open file (" << ReadInput << ")" << endl;  exit(1);}
      fp_in >> imode;
      if(imode==1){
      fp_in >> vis >>  vist  >> visb;
      fp_in >> lambda1 >> lambda2;
      fp_in >>  nsplit1 >> nsplit2 >> nsplit3 >> nsplit4;
      fp_in >> iatm;
      fp_in >> switchR >> switchA;
      fp_in >> itypeu >> itypes >> itypet >> itypeb >> ocean;
      fp_in >> ihumid >> know >> groundtemp >> groundpress;
      fp_in >> aeroplane >> heightplane;
      fp_in >> default_pause >> HG >> HB >> HT >> HS >> HU;
      fp_in >> ngaussL >> ngaussH;
     fp_in.close();
     }
     else{  
        cout << "imode !=1 not coded yet\n"; exit(1);
     }

  }
  else  // then verbose
 {


  FILE *fp0;
  fp0=fopen("Dors1in.dat", "w");
  cout << "Enter 1 for Earth Atmosphere\n";
  cout << "Enter 2 for Arbitrary Atmosphere\n";
  cin >> imode;
     if(imode==1){

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
  cout << "Enter the total boundary layer optical depth\n";
  cin >>  visb;
  cout << "Enter the troposphere  optical depth (exluding bdry layer) \n";
  cin >>  vist;}

  cout << "Enter lowest wavelength value\n";
  cin >> lambda1;
  cout << "Enter highest wavelength value\n";
  cin >> lambda2;

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
               cout << "Enter ground temperature\n";
               cin >> groundtemp;}
               cout << "Enter 1 if you know the ground pressure, 0 otherwise\n";
               cin >> ipress;
               if(ipress==1){
               cout << "Enter ground pressure\n";
               cin >> groundpress;}
               }
        }


     cout << "Enter 1 for aeroplane remote sensing, enter 0 otherwise\n";
     cin >> aeroplane;

     if(aeroplane){
         cout << "Enter height of aeroplane\n";
         cin >> heightplane;
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
        
        cout << "Enter order of Gauss-Legendre quadrature rule\n";
        cin >> ngaussL;
        cout << "Enter order of Gauss-Hermite rule\n";
        cin >> ngaussH;

     fprintf(fp0, "%d \n", imode);
     fprintf(fp0, "%d %lf %lf\n", vis, vist, visb);
     fprintf(fp0, "%lf %lf \n", lambda1, lambda2);
     fprintf(fp0, "%d %d %d %d\n", nsplit1, nsplit2, nsplit3, nsplit4);
     fprintf(fp0, "%d\n", iatm);
     fprintf(fp0, "%d %d\n", switchR, switchA);
     fprintf(fp0, "%d %d %d %d %d\n", itypeu, itypes, itypet, itypeb, ocean);
     fprintf(fp0, "%d %d %lf %lf\n", ihumid, know, groundtemp, groundpress);
     fprintf(fp0, "%d %lf\n", aeroplane, heightplane);
     fprintf(fp0, "%d %lf %lf %lf %lf %lf\n", default_pause, HG, HB, HT, HS, HU);
     fprintf(fp0, "%d %d\n", ngaussL, ngaussH);
     fclose(fp0);}
     else{
        cout << "imode !=1 not coded yet\n";
         exit(1);
     }   //endif for imode=1

     }  //end for else verbose
        

}




