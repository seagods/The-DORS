#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

// slatec fortran cubic spline routines
extern "C" {

void dbint4_(double*, double*, int*, int*, int*,  double*, double*,
              int*, double*, double*, int*, int*, double*);
             

double dbvalu_(double*, double*, int* ,int*,int* ,double*, int*, double*);

//toms library Gauss quadrature
void weightcoeff_(int*, double*, double*, double*, double*,
                       double*, double*);
}

void Questioner( bool&,
                int&, int&, int&, int&, bool&, double&,
                double&, double&, double&, bool&,  int&, bool&, bool&, int&,
                int&, int&, int&, bool&, bool&, bool&, 
                int&, double&, double&, bool&, double&,
                bool&, double&, double&, double&, double&, double&, int&, 
                bool&, bool &, int&, double&, double&,
                const char*);

void WaterVap(double* , double* , double* ,
              double* , double &, double &, double &, int &, int &, int &);
              

int main(int argc, char* argv[]){

  bool verbose;   //verbose mode for input, read from file if false;
  bool verbose_out=true;  // output loads of info to standard output
  bool verbose_height=false;   //output height integration info to std::out
/**************************************************************************************
   Stage 1:  declare and define various data to define the atmosphere
             read response function, read atmospheric definitions
             get gauss quadrature rule
             declare and initialise physical constants
             read molparam, gasmask, isotopolgues, and set PPM files to read later
             read in basic atmospheric data
****************************************************************************************
   Stage 2:  determine numbers of slabs in the 4 regions
             Reverse the presssure profile so top down             
             interpolate x axis is pressure, y axis is height
             determine the four sets of equal pressures and the heights of all the "floors"
             Find the fine grid points for gauss quadrature w.r.t height
             Interpolate data onto fine grid
****************************************************************************************
   Stage 3:  Do the integrations w.r.t. height for masses, number densities, average
             temperature pressure and pressure for all the nsplit1+nsplit4+nsplit3+nsplit4 sub-slabs.
****************************************************************************************
   Stage 4:  Loop over Hitran par files, read in line centres, calculate the spectrum
****************************************************************************************/


/************************   BEGIN STAGE 1  *********************************************/
  // Four Aerosol layers, each to be split into nsplit slabs
  // Four layers are boundary layer, troposphere, stratosphere, and upper atmosphere
  int nsplit1,nsplit2,nsplit3,nsplit4;    
  // boundary layer visibility , troposphere visiblility, wavelength1, wavelength2
  double visb, vist, lambda1, lambda2;
  bool rfun=true;  //true if we are using myresp.ddat or myrespX.dat as a response function file
  bool vis;  //true if entered as visibility, false if entered as optical depth.
  // default atmosphere 1-6 or user defined
  int iatm;    //atmosphere ID
  bool switchR;    //Switch Rayleigh scattering on/off
  bool switchA;    //Switch aerosol scattering on/off
  int itypeu, itypes, itypet, itypeb;    //types for upper atmosphere, stratosphere, troposphere, bdry layer
  bool ocean=false, groundT=false, groundP=false;  //ocean or maritime, know ground temperature, know ground pressure
  int  ihumid;      // use default humid or user defined
  double  groundtemp, groundpress;  //values for ground pressure and temperature (negative if unknown)
  bool cloud; // do we have clouds
  bool aeroplane;  //do we want to simulate aeroplane remote sensing
  double heightplane;  //negative if no aeroplane

  bool default_pause;   //use default stratopause, topopause etc
  double HG, HB, HT, HS, HU;  //height of ground, bdry layer, tropopause, stratosphere
  int ngauss_height, ngauss_correlkay;       //order of Gauss-Legendre quadrature height and cumulative distribution
  bool calcspec, outspec;              //if false, don't bother with calculating the 
  double nu_cut, nu_Delta;
  const char* ReadInput;      //name of input file
  const char* RespFile;       //name of response function file

   // First off we read in the user input, whether 
   // from command line (Dors1 -v) or an input file (Dors1)



   /*****************************************************************/
   // NOTE -NO SAFETY BELTS, THE USER MUST MAKE SURE THE INPUT IS O.K.
   /*****************************************************************/
   //Questions the user for atmospheric types and so on
   //or reads data from file if verbose is false

    if(argc==1)
     verbose=false;

     ReadInput="Dors1in.dat";
     RespFile="myresp.dat";

     if(argc !=1)
     {
       if(argc==2){
         //it was either Dors1 -v or Dors1 Dors1inX.dat
         if(argv[1][0]=='-' && argv[1][1]=='v'){
           verbose=true;
           ReadInput="Dors1in.dat"; //Not used, but initialise anyway 
           RespFile="myresp.dat";
          }
         else
         {
             ReadInput=argv[1];
             RespFile="myresp.dat";
             cout << "Reading input from file=" << ReadInput << endl;
         }
       }  //endif argc==2

         if(argc==3){
          // It should have been Dors1 Dors1inX.dat myrespX.dat
             ReadInput=argv[1];
             cout << "Reading input from file=" << ReadInput << endl;
             RespFile=argv[2];
             cout << "Response function from file=" << RespFile<< endl;
         }

         if(argc>4){
            // presumably some garbage or other crept in
            cout << "Usage Dors1 or Dors1 -v" << endl;
            cout << "argc=" << argc << endl;
            exit(2);
         }
    }  //endif argc !=1


   Questioner(verbose,
              nsplit1, nsplit2, nsplit3, nsplit4, vis, visb,
              vist, lambda1, lambda2, rfun,  iatm, switchR, switchA, itypeu,
              itypes, itypet, itypeb, ocean, groundP, groundT,
              ihumid, groundtemp, groundpress, aeroplane, heightplane,
              default_pause, HG, HB, HT, HS, HU, ngauss_height, 
              calcspec, outspec, ngauss_correlkay, nu_cut, nu_Delta,
              ReadInput);  

    int nresp;
    double* lambda_resp; double* respval;
    if(rfun){
       ifstream fp_in;
       fp_in.open(RespFile, ios::in);
       if(!fp_in.is_open()){ cerr << "Failed to open file (" << RespFile << ")"  << endl; exit(1); }
       fp_in >> nresp;
        lambda_resp=new double[nresp]; respval=new double[nresp];
        for(int i=0; i< nresp ; i++){
         fp_in >> lambda_resp[i] >> respval[i];
        }
        lambda1=lambda_resp[0]; lambda2=lambda_resp[nresp-1];
        fp_in.close();
        cout << "Response function file overrides wavelengths in Dors1in.dat\n";
        cout << "Wavelengths used are " << lambda1 << " and " << lambda2 << endl;
    }



   // Use Gauss quadrature routine in toms library

 
   /* *********************************************************** */
   /*        Quadrature Rule for height integration               */
   /* *********************************************************** */
    double Q[ngauss_height],E[ngauss_height],Xgauss[ngauss_height],Wgauss[ngauss_height],Work_g[9*ngauss_height+8];
    double EPS, Xtemp[ngauss_height]; //needed for Gauss quadrature Xgauss=knots, WGauss are the weights.

    EPS=1e-15;
    for(int i=2;i<=ngauss_height;i++){
       Q[i-1]=2.*i*i/(2.*i*(2.0*i-1.0));
       E[i-1]=2.*i*i/(2.*i*(2.0*i+1.0));
    }
    Q[0]=1.0;
    E[0]=1.0/3.0;
    //toms algorithm 125
    weightcoeff_(&ngauss_height,Q,E,&EPS,Wgauss,Xgauss,Work_g);
    // why?

    //for some reason the x interval is [0,2]
    //and you have to double the weights if you want [0,2]
    //not only that the X are in reverse order
    //now transform to [0,1] 

    for(int i=1; i <=ngauss_height; i++){
       Xgauss[i-1]=Xgauss[i-1]/2.0;
       Xtemp[i-1]=Xgauss[i-1];
    } 
    for(int i=0; i <ngauss_height; i++){
        Xgauss[i]=Xtemp[ngauss_height-i-1];
    }
            
   /* *********************************************************** */
   /*        Quadrature Rule for Correletad kay                   */
   /* *********************************************************** */
    double QK[ngauss_correlkay],EK[ngauss_correlkay],XgaussK[ngauss_correlkay],WgaussK[ngauss_correlkay],Work_gK[9*ngauss_correlkay+8];
    double EPSK, XtempK[ngauss_correlkay]; //needed for Gauss quadrature Xgauss=knots, WGauss are the weights.

    EPSK=1e-15;
    for(int i=2;i<=ngauss_correlkay;i++){
       QK[i-1]=2.*i*i/(2.*i*(2.0*i-1.0));
       EK[i-1]=2.*i*i/(2.*i*(2.0*i+1.0));
    }
    QK[0]=1.0;
    EK[0]=1.0/3.0;
    //toms algorithm 125
    weightcoeff_(&ngauss_correlkay,QK,EK,&EPSK,WgaussK,XgaussK,Work_gK);
    // why?

    //for some reason the x interval is [0,2]
    //and you have to double the weights if you want [0,2]
    //not only that the X are in reverse order
    //now transform to [0,1] 

    for(int i=1; i <=ngauss_correlkay; i++){
       XgaussK[i-1]=XgaussK[i-1]/2.0;
       XtempK[i-1]=XgaussK[i-1];
    } 
    for(int i=0; i <ngauss_correlkay; i++){
        XgaussK[i]=XtempK[ngauss_correlkay-i-1];
    }           
   /* *********************************************************** */
   /*        Quadrature Rule for correlate k done                 */
   /* *********************************************************** */

   //Some constants

   double pi,abszero,R,AV,MWA,delta,nucleon;

   pi=acos(-1.0);               //Well pi obviusly
   abszero=273.15;              //0 degrees C in Kelvin
   R=8.314472;                  //Universal gas constant
   AV=6.0221415e23;             //Avogadro's number;
   MWA=28.946;                  //Molecular Weight of Air
   delta=0.029;                 //depolarisation factor
   nucleon=1.66053886E-24;      //"nucleon mass" (a.m.u) in grammes


   ifstream fp_in;   //input file stream
   int gid[50];      //HITRAN gas numbers (last few are dummies)
   bool gmask[50];   //if false, ignore that gas
   fp_in.open("GasMask.dat", ios::in);
   if(!fp_in.is_open()){ cout << "Failed to open file (GasMask.dat)" << endl;  exit(1);}
   fp_in >> gid[0] >> gid[1] >> gid[2] >> gid[3] >> gid[4] >>
            gid[5] >> gid[6] >> gid[7] >> gid[8] >> gid[9];
   fp_in >> gmask[0] >> gmask[1] >> gmask[2] >> gmask[3] >> gmask[4] >>
            gmask[5] >> gmask[6] >> gmask[7] >> gmask[8] >> gmask[9];
   fp_in >> gid[10] >> gid[11] >> gid[12] >> gid[13] >> gid[14] >>
            gid[15] >> gid[16] >> gid[17] >> gid[18] >> gid[19];
   fp_in >> gmask[10] >> gmask[11] >> gmask[12] >> gmask[13] >> gmask[14] >>
            gmask[15] >> gmask[16] >> gmask[17] >> gmask[18] >> gmask[19];
   fp_in >> gid[20] >> gid[21] >> gid[22] >> gid[23] >> gid[24] >>
            gid[25] >> gid[26] >> gid[27] >> gid[28] >> gid[29];
   fp_in >> gmask[20] >> gmask[21] >> gmask[22] >> gmask[23] >> gmask[24] >>
            gmask[25] >> gmask[26] >> gmask[27] >> gmask[28] >> gmask[29];
   fp_in >> gid[30] >> gid[31] >> gid[32] >> gid[33] >> gid[34] >>
            gid[35] >> gid[36] >> gid[37] >> gid[38] >> gid[39];
   fp_in >> gmask[30] >> gmask[31] >> gmask[32] >> gmask[33] >> gmask[34] >>
            gmask[35] >> gmask[36] >> gmask[37] >> gmask[38] >> gmask[39];
   fp_in >> gid[40] >> gid[41] >> gid[42] >> gid[43] >> gid[44] >>
            gid[45] >> gid[46] >> gid[47] >> gid[48] >> gid[49];
   fp_in >> gmask[40] >> gmask[41] >> gmask[42] >> gmask[43] >> gmask[44] >>
            gmask[45] >> gmask[46] >> gmask[47] >> gmask[48] >> gmask[49];
   fp_in.close();
   const char* AllMols[50];   //File names for reading gas number density profiles
   if(iatm==1)AllMols[0]="Atmospheres/Mols1_H2O.dat";if(iatm==2)AllMols[0]="Atmospheres/Mols2_H2O.dat";
   if(iatm==3)AllMols[0]="Atmospheres/Mols3_H2O.dat";if(iatm==4)AllMols[0]="Atmospheres/Mols4_H2O.dat";
   if(iatm==5)AllMols[0]="Atmospheres/Mols5_H2O.dat";if(iatm==6)AllMols[0]="Atmospheres/Mols6_H2O.dat";
   if(iatm==1)AllMols[2]="Atmospheres/Mols1_O3.dat";if(iatm==2)AllMols[2]="Atmospheres/Mols2_O3.dat";
   if(iatm==3)AllMols[2]="Atmospheres/Mols3_O3.dat";if(iatm==4)AllMols[2]="Atmospheres/Mols4_O3.dat";
   if(iatm==5)AllMols[2]="Atmospheres/Mols5_O3.dat";if(iatm==6)AllMols[2]="Atmospheres/Mols6_O3.dat";
   if(iatm==1)AllMols[3]="Atmospheres/Mols1_N2O.dat";if(iatm==2)AllMols[3]="Atmospheres/Mols2_N2O.dat";
   if(iatm==3)AllMols[3]="Atmospheres/Mols3_N2O.dat";if(iatm==4)AllMols[3]="Atmospheres/Mols4_N2O.dat";
   if(iatm==5)AllMols[3]="Atmospheres/Mols5_N2O.dat";if(iatm==6)AllMols[3]="Atmospheres/Mols6_N2O.dat";
   if(iatm==1)AllMols[4]="Atmospheres/Mols1_CO.dat";if(iatm==2)AllMols[4]="Atmospheres/Mols2_CO.dat";
   if(iatm==3)AllMols[4]="Atmospheres/Mols3_CO.dat";if(iatm==4)AllMols[4]="Atmospheres/Mols4_CO.dat";
   if(iatm==5)AllMols[4]="Atmospheres/Mols5_CO.dat";if(iatm==6)AllMols[4]="Atmospheres/Mols6_CO.dat";
   if(iatm==1)AllMols[5]="Atmospheres/Mols1_CH4.dat";if(iatm==2)AllMols[5]="Atmospheres/Mols2_CH4.dat";
   if(iatm==3)AllMols[5]="Atmospheres/Mols3_CH4.dat";if(iatm==4)AllMols[5]="Atmospheres/Mols4_CH4.dat";
   if(iatm==5)AllMols[5]="Atmospheres/Mols5_CH4.dat";if(iatm==6)AllMols[5]="Atmospheres/Mols6_CH4.dat";
   if(iatm==7){
   AllMols[0]="Atmospheres/Mols7_H2O.dat";  AllMols[2]="Atmospheres/Mols7_O3.dat";
   AllMols[3]="Atmospheres/Mols7_N2O.dat";  AllMols[4]="Atmospheres/Mols7_CO.dat";
   AllMols[5]="Atmospheres/Mols7_CH4.dat";  
   }
   AllMols[1]="Atmospheres/Mols_CO2.dat";
   AllMols[6]="Atmospheres/Mols_O2.dat";  AllMols[7]="Atmospheres/Mols_NO.dat";
   AllMols[8]="Atmospheres/Mols_SO2.dat";  AllMols[9]="Atmospheres/Mols_NO2.dat";
   AllMols[10]="Atmospheres/Mols_NH3.dat";  AllMols[11]="Atmospheres/Mols_HNO3.dat";
   AllMols[12]="Atmospheres/Mols_OH.dat";  AllMols[13]="Atmospheres/Mols_HF.dat";
   AllMols[14]="Atmospheres/Mols_HCl.dat";  AllMols[15]="Atmospheres/Mols_HBr.dat";
   AllMols[16]="Atmospheres/Mols_HI.dat";  AllMols[17]="Atmospheres/Mols_ClO.dat";
   AllMols[18]="Atmospheres/Mols_OCS.dat";  AllMols[19]="Atmospheres/Mols_H2CO.dat";
   AllMols[20]="Atmospheres/Mols_HOCl.dat";  AllMols[21]="Atmospheres/Mols_N2.dat";
   AllMols[22]="Atmospheres/Mols_HCN.dat";  AllMols[23]="Atmospheres/Mols_CH3Cl.dat";
   AllMols[24]="Atmospheres/Mols_H2O2.dat";  AllMols[25]="Atmospheres/Mols_C2H2.dat";
   AllMols[26]="Atmospheres/Mols_C2H6.dat";  AllMols[27]="Atmospheres/Mols_PH3.dat";
   AllMols[28]="Atmospheres/Mols_COF2.dat";  AllMols[29]="Atmospheres/Mols_SF6.dat";
   AllMols[30]="Atmospheres/Mols_H2S.dat";  AllMols[31]="Atmospheres/Mols_HCOOH.dat";
   AllMols[32]="Atmospheres/Mols_HO2.dat";  AllMols[33]="Atmospheres/Mols_O.dat";
   AllMols[34]="Atmospheres/Mols_ClONO2.dat";  AllMols[35]="Atmospheres/Mols_NO+.dat";
   AllMols[36]="Atmospheres/Mols_HOBr.dat";  AllMols[37]="Atmospheres/Mols_C2H4.dat";
   AllMols[38]="Atmospheres/Mols_CH3OH.dat";  AllMols[39]="Atmospheres/Mols_CH3Br.dat";
   AllMols[40]="Atmospheres/Mols_CH3CN.dat";  AllMols[41]="Atmospheres/Mols_CF4.dat";
   AllMols[42]="Atmospheres/Mols_X1.dat";  AllMols[43]="Atmospheres/Mols_X2.dat";
   AllMols[44]="Atmospheres/Mols_X3.dat";  AllMols[45]="Atmospheres/Mols_X4.dat";
   AllMols[46]="Atmospheres/Mols_X5.dat";  AllMols[47]="Atmospheres/Mols_X6.dat";
   AllMols[48]="Atmospheres/Mols_X7.dat";  AllMols[49]="Atmospheres/Mols_X8.dat";

   fp_in.open("Atmospheres/AllIsos.dat",ios::in);
   int nice[50]; // number of isotope variations (isotopolgues) for  each gas  molecule
   fp_in >>nice[0]>>nice[1]>>nice[2]>>nice[3]>>nice[4]>>nice[5]>>nice[6]>>nice[7]>>nice[8]>>nice[9];
   fp_in >>nice[10]>>nice[11]>>nice[12]>>nice[13]>>nice[14]>>nice[15]>>nice[16]>>nice[17]>>nice[18]>>nice[19];
   fp_in >>nice[20]>>nice[21]>>nice[22]>>nice[23]>>nice[24]>>nice[25]>>nice[26]>>nice[27]>>nice[28]>>nice[29];
   fp_in >>nice[30]>>nice[31]>>nice[32]>>nice[33]>>nice[34]>>nice[35]>>nice[36]>>nice[37]>>nice[38]>>nice[39];
   fp_in >>nice[40]>>nice[41]>>nice[42]>>nice[43]>>nice[44]>>nice[45]>>nice[46]>>nice[47]>>nice[48]>>nice[49];
   fp_in.close();
   //The relative abundance of each isotopologues are held in decreasing order
   //Many of these are very small indeed and may be ignored for most purposes
   //we introduce a cutoff, so if nicemax=1 only the most common is isotopologue
   //and so on. CO2 has the most isotopolgues, equal to ten.
   int nicemax=2;
   
   int nicecode[50][10];  //HITRAN isotope code for molecule
   double AllMolW[50][10]; //molecular weights for isotopologue
   double niceabund[50][10];  //relative abundance of isotopologue
   double QTref[50][10];      //partition function at reference temperature 296K
   int geejay[50][10];        //state degeneracy factor
   //initialise
   for(int i=0; i<50; i++){
     for(int j=0; j<10; j++){
         nicecode[i][j]=-1;
         AllMolW[i][j]=-1;
         niceabund[i][j]=-1;
         geejay[i][j]=-1;
         QTref[i][j]=-1;}}
   fp_in.open("HITRAN/molparam.txt",ios::in);
   string line1;
   getline(fp_in, line1);  //get first line, no info, then one getline
                         //for each molecule (no info)
   cout << line1 << endl;
   int  imol=0;
   while(!getline(fp_in, line1).eof()){
      cout << line1 << endl;
      for(int iso=0; iso<nice[imol];iso++){
         fp_in >>  nicecode[imol][iso] >> niceabund[imol][iso] >> QTref[imol][iso]
               >> geejay[imol][iso] >> AllMolW[imol][iso];
         if(verbose_out)cout << imol+1 << " " << nice[imol] << "  " << nicecode[imol][iso] << "  " << AllMolW[imol][iso] << endl;
      }
      if(!getline(fp_in,line1)){ break;}
       // << "  getline return\n";
      cout << line1;
     // if(getline(fp_in, line1).eof()){ break;}
      imol++;
   }
   fp_in.close(); //molparam.txt finished with

   int ngas=0;  //the number of gases not masked out.
   for(int  i=0;i<50;i++){
       if(gmask[i])ngas++;}   


   int gases[ngas];  //HITRAN gas IDs
   //MolFiles are just the Atmospheres/MolsX files actually used
   //rather than (but might possibly be) the whole set 
   const char* MolFiles[ngas];      //A copy of AllMols, but only containing the gases actually taken into account
   double MWGAS[ngas];              //The MEAN molecular weight of each gas, taking into account the
                                    // isotopologues and abundances

   ngas=0;     //reset to zero now declarations have been made
   for(int i=0;i<50;i++){
       if(gmask[i]){
           gases[ngas]=gid[i];
           MolFiles[ngas]=AllMols[i];
           //Assume natural abundances from molparams.txt
           MWGAS[ngas]=0.0;
           for(int j=0; j<nice[i]; j++){
              MWGAS[ngas]=MWGAS[ngas]+niceabund[i][j]*AllMolW[i][j];
           }
           ngas++;                           //recount ngas
           }}
      cout << "iatm=" << iatm << endl;         //output as user check

      for(int i=0; i<ngas;i++){
      if(verbose_out)cout << gases[i] << "  " << MolFiles[i] << "   MW=" << MWGAS[i] << endl;
          }
             
/* ********************************************************************* */
//    Points for integrating over height - first read in the data
/* ********************************************************************* */

   fp_in.open("Atmospheres/Alts.dat", ios::in);
   if(!fp_in.is_open()){ cout << "Failed to open file (Atmospheres/Alts.dat)" << endl;  exit(1);}

   int i_altitudes;     //i_altitudes is the number of altitudes defining the temperature and pressure profiles
   fp_in >> i_altitudes;
   double altitude[i_altitudes];    //altitude[] contains these altitudes
   for(int i=0;i<i_altitudes;i++){
     fp_in >> altitude[i];
     altitude[i]=altitude[i]*1000.0;  //convert to metres
     }
    fp_in.close();

   char sub[3];
   sprintf(sub,"%d",iatm);   //we shall substitute the atmosphere number into the strings for the input files

   string Tempsname="Atmospheres/Temps .dat";
   // stringreplace(index, num1, const char*, num2)
   //index tells it where, up to  num1 characters repleced with up to num 2
   Tempsname.replace(17,1,sub,1);

   fp_in.open(Tempsname.c_str(), ios::in);
   if(!fp_in.is_open()){ cout << "Failed to open file, file= " <<  Tempsname.c_str() << endl;  exit(1);}

   int idata;    //Make sure the number of temperatures (idata) is the same as i_altitudes
   fp_in >> idata;
   if(idata !=i_altitudes){
        cout << "Temperature data does not match altitudes\n";  exit(1);}
   double T[idata];

   for(int i=0;i<idata;i++){
      fp_in >>  T[i];
   }
   fp_in.close();  // now do the same for the pressure data

   string Pressname="Atmospheres/Press .dat";
   Pressname.replace(17,1,sub,1);
   fp_in.open(Pressname.c_str(), ios::in);
   if(!fp_in.is_open()){ cout << "Failed to open file, file= " <<  Pressname.c_str() << endl;  exit(1);}

   fp_in >> idata;
   if(idata !=i_altitudes){
         cout << " Pressure data does not match altitude data\n" << endl; exit(1);}

   double P[idata];      //Pressures at the input heights
   for(int i=0;i<idata;i++){
      fp_in >>  P[i] ;
    //  cout << i << "  " << P[i] << endl;
      }
    fp_in.close();

    // If we have non default water vapour via a humidity file we
    // need to calculate the number density of the water vapour molecules, store in Gee1 and Gee2
    double Gee1[idata];  //Gee1 only needed if using Hum1_7.dat data 
    double Gee2[idata];  //Gee2 only needed if using Hum2_7.dat data
                         //We shall need humidity for aerosol hygroscopic growth
    double Hum1[idata];  //Hum1 and Humids1 needed to store humidity values over liquid water
    double Hum2[idata];  //Hum2 and Humids2 needed to store humidity values over ice
    double GasPPM[idata][ngas];  //Gas concentration of the selected gasess in ppm
    double x1;
    //now read in data for gases
    for(int i=0;i<ngas;i++){
       fp_in.open(MolFiles[i],ios::in);  //open the HITRAN molecular profiles
       if(!fp_in.is_open()){cout << "Failed to open file, file= " <<  MolFiles[i] << endl;  exit(1);}
       fp_in >> idata;
       if(idata!=i_altitudes){cout << MolFiles[i] << "data does not match altitude data\n"; exit(1);} 
       for(int j=0; j< idata; j++){
             fp_in >> x1;
             GasPPM[j][i]=x1;
             if((gases[i]==1) && (ihumid==1) ){Gee1[j]=GasPPM[j][i]; Gee2[j]=GasPPM[j][i];}
       //      cout << i << " " << j << "  " << MolFiles[i] << "  " << GasPPM[j][i] << endl;

             }
        fp_in.close();
    }
    int ig;

       ig=0;  //Water Vapour
       int iwi; // 1 for equilibrium with water, 2 for equilibrium with ice
       int idataW1, idataW2;
       if(ihumid==2){    //not using default profiles but user humidity profiles instead
         ifstream fp_humid;
         //CHANGE FOR 1-6 +User  later
         fp_humid.open("Atmospheres/Hum1_U.dat", ios::in);
         fp_humid >> idataW1;
         if(idataW1 != i_altitudes){
              cout <<"Humidity Files don't match altitude\n"; exit(1);}
         for(int i=0; i <idataW1; i++){ 
             fp_humid >> Hum1[i];}
         fp_humid.close();
         //note that humidity files with idata=1 or 2 recommended
         //then GasPPM for H20 only changed at the ground or first km


         fp_humid.open("Atmospheres/Hum2_U.dat", ios::in);
         fp_humid >> idataW2;
         if(idataW1 != idataW2){
              cout <<"Humidity Files don't match\n"; exit(1);}
         for(int i=0; i <idataW1; i++){ 
             fp_humid >> Hum2[i];}
         fp_humid.close();
          //have two humidity humidity profiles
          //since ihumid=2 convert to ppm via function WaterVap 

         iwi=1;
         WaterVap(P,T,Gee1,Hum1,AV,nucleon,MWGAS[0],idataW1,ihumid,iwi);
         iwi=2;
         WaterVap(P,T,Gee2,Hum2,AV,nucleon,MWGAS[0],idataW1,ihumid,iwi);

         for(int i=0; i <idataW1; i++){ 
             if(T[i]>abszero){  // remember 0 Celcius=abszero Kelvin
                 GasPPM[i][ig]=Gee1[i];}
             else{
                 GasPPM[i][ig]=Gee2[i];}
        }
           //no longer want or need ihumid=2
           ihumid=1;
        }  // end for using Hum1 and Hum2 files and calculating GasPPM for H20
        else{
           // we calculate Hum1 and Hum2 from GasPPM for H20 since ihumid=1
           // needed for aerosols later, 
           iwi=1;
           WaterVap(P,T,Gee1,Hum1,AV,nucleon,MWGAS[0],idata,ihumid,iwi);
           iwi=2;
           WaterVap(P,T,Gee2,Hum2,AV,nucleon,MWGAS[0],idata,ihumid,iwi);

//         UNCOMMENT THESE LINES IF YOU WANT HUMIDITY OUTPUT
//           ofstream Hum1file, Hum2file;
//            Hum1file.open("Temp_Humid1.dat",ios::out);
//            Hum2file.open("Temp_Humid2.dat",ios::out);
//           Hum1file << idata << endl; Hum2file << idata << endl;
//           for(int ihum=0; ihum< idata; ihum++){
//              Hum1file << Hum1[ihum] << endl; Hum2file << Hum2[ihum] << endl;
//          }

        }
  
      if(groundP){
            //scale temperature and pressure to ground;
            double tempP=groundpress/P[0];
            for(int i=0; i< idata; i++){
               P[i]=P[i]*tempP;  
              }}
            if(groundT){
            double tempT=groundtemp/T[0];
            T[0]=T[0]*tempT;
            T[1]=T[1]*tempT;    // only scale lower atmosphere
           }


    // We need to integrate quanties in each sublayer, get gauss
    // quadrature weights and coefficients on (0,1)
/***************************************************************************************
                          END STAGE 1
****************************************************************************************
                         BEGIN STAGE 2
              Calculate integration nodes and interpolate data  to them
****************************************************************************************/
        //need T and P data etc at fine resolution (heights)
        // Spline from T and P data etc at
        //declare P2 and T2 idata long for spline, store fine data
        // in Temp, Press, etc

        //prepare for calls to dbint4.f
        //
        double HBDRY[5], PBDRY[5];   //Heights and pressures at 4 main boundaries
        
        HBDRY[0]=HG*1000.0; HBDRY[1]=HB*1000.0; HBDRY[2]=HT*1000.0;
        HBDRY[3]=HS*1000.0; HBDRY[4]=HU*1000.0; 
          
        int IBCL,IBCR,IN,KORDER,KNOPT,ideriv,inbv;
        double FBCL,FBCR,estim;

        //input 
        //WORK- just need to size array
        
        //output
        //KORDER  (will be 4)
        //IN  (number of coefficients, will be ndata+2)
        //BCOEFF

        KORDER=4; // order of spline
        KNOPT=1;  //option for no extrapolation outside data
        IBCL=2; IBCR=2;   //natural spline
        FBCL=0.0; FBCR=0.0;  //natural spline
        IN=idata+2;
    
        double TEE[idata+6];
        double BCOEF[IN];
        //for calculating splines
        double WORK[5*(idata+2)]; double WORK2[3*KORDER];
        
        dbint4_(altitude,P,&idata,&IBCL,&IBCR,&FBCL,&FBCR,
                       &KNOPT, TEE, BCOEF,&IN, &KORDER, WORK);

        ideriv=0; inbv=1; //inbv must always be set to 1 on first call to dbvalu
        for(int i=0; i<5;i++){
           PBDRY[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv, HBDRY+i, &inbv,WORK2);  //find pressures at the boundaries
         }

        //Need to find heights at equal pressure intervals
        //
        int up_to_top, up_to_up, up_to_strat, up_to_trop;  //numbers of slabs from ground to each bdry
        up_to_top=nsplit4+nsplit3+nsplit2+nsplit1;
        up_to_up=nsplit4+nsplit3+nsplit2;
        up_to_strat=nsplit4+nsplit3;
        up_to_trop=nsplit4;
        double floorpress[up_to_top];    //pressures at the "floor" of each slab
        double HeightFloors[up_to_top];   //heights of those floors
        double ReverseP[idata];           // i_altitude pressures top down
        double ReverseH[idata];           //altitude top down
       //     cout  << "Reverse " << idata <<  endl;
        for(int i=0; i< idata; i++){
            ReverseP[i]=P[idata-1-i];
            ReverseH[i]=altitude[idata-1-i];
         }


        //ground to troposphere
        for(int i=0; i<nsplit4; i++){
           floorpress[i]=PBDRY[0]
                        +(PBDRY[1]-PBDRY[0])*((double) i)/( (double) nsplit4);
        }
        //troposphere to stratosphere
        for(int i=0; i<nsplit3; i++){
           floorpress[up_to_trop+i]=PBDRY[1]
                        +(PBDRY[2]-PBDRY[1])*((double) i)/( (double) nsplit3);
        }
        //stratosphere to upper atmosphere
        for(int i=0; i<nsplit2; i++){
           floorpress[up_to_strat+i]=PBDRY[2]
                        +(PBDRY[3]-PBDRY[2])*((double) i)/( (double) nsplit2);
        }
        //stratosphere to upper atmosphere
        for(int i=0; i<nsplit1; i++){
           floorpress[up_to_up+i]=PBDRY[3]
                        +(PBDRY[4]-PBDRY[3])*((double) i)/( (double) nsplit1);
        }
        //spline altitude as a function of pressure
        inbv=1; //inbv must always be set to 1 on first call to dbvalu
        dbint4_(ReverseP,ReverseH,&idata,&IBCL,&IBCR,&FBCL,&FBCR,
                       &KNOPT, TEE, BCOEF,&IN, &KORDER, WORK);

        //find the heights of the bottom of every sub-slab in the atmosphere
        //by sub-slab I mean slabs within the four aerosol regions
        for(int i=0; i<up_to_top;i++)
           HeightFloors[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv,floorpress+i,&inbv,WORK2);


        ideriv=0,inbv=1; //inbv must always be set to 1 on first call to dbvalu
        //spline temperature as a function of height
        dbint4_(altitude,T,&idata,&IBCL,&IBCR,&FBCL,&FBCR,
                       &KNOPT, TEE, BCOEF, &IN, &KORDER, WORK);

        double FloorTemp[up_to_top];   //Temperature at the floor of each sub-slab
        double SlabThick[up_to_top];   //thickness of each slab
        for(int i=0; i<up_to_top;i++){
           FloorTemp[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv,HeightFloors+i,&inbv,WORK2);

        //we shall do all calculations ground up
        // but output them top down  --- HeightFloors is top down
           if(i==up_to_top-1){
              SlabThick[i]=HBDRY[4]-HeightFloors[i];
                if(verbose_height)cout << "Floors " << HeightFloors[i] << " "  <<floorpress[i] << " " << FloorTemp[i] <<
               " " << SlabThick[i] << endl;
               }
              else{
              SlabThick[i]=HeightFloors[i+1]-HeightFloors[i];
                if(verbose_height)cout << "Floors " << HeightFloors[i] << " "  <<floorpress[i] << " " << FloorTemp[i] <<
                " " << SlabThick[i] << endl;
               }
         } //end i<up_to_top loop

        //Calc Gauss integration heights ground up
        int iheight;     //fine grid for temperature and pressure
        iheight=ngauss_height*up_to_top;
        double heights[iheight];
        int kount=0;
        for(int i=0; i<nsplit4; i++){
           for(int j=0; j<ngauss_height; j++){
           heights[kount]=HeightFloors[i]+SlabThick[i]*Xgauss[j];
           if(verbose_height)cout <<"kount=" << kount << "  heights=" << heights[kount] 
           kount++;
        }}
        for(int i=0; i<nsplit3; i++){
           for(int j=0; j<ngauss_height; j++){
           heights[kount]=HeightFloors[up_to_trop+i]+SlabThick[up_to_trop+i]*Xgauss[j];
           if(verbose_height)cout <<"kount=" << kount << "  heights=" << heights[kount] << endl;
           kount++;
        }}
        for(int i=0; i<nsplit2; i++){
           for(int j=0; j<ngauss_height; j++){
           heights[kount]=HeightFloors[up_to_strat+i]+SlabThick[up_to_strat+i]*Xgauss[j];
           if(verbose_height)cout <<"kount=" << kount << "  heights=" << heights[kount] << endl;
           kount++;
        }}
        for(int i=0; i<nsplit1; i++){
           for(int j=0; j<ngauss_height; j++){
           heights[kount]=HeightFloors[up_to_up+i]+SlabThick[up_to_up+i]*Xgauss[j];
           if(verbose_height)cout <<"kount=" << kount << "  heights=" << heights[kount] << endl;
           kount++;
        }}
           
        double Temp[iheight],Press[iheight];  //temperature and pressure on fine grid
        for(int i=0; i<iheight;i++){
           Temp[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv,heights+i,&inbv,WORK);}
        ideriv=0,inbv=1; //inbv must always be set to 1 on first call to dbvalu
        dbint4_(altitude,P,&idata,&IBCL,&IBCR,&FBCL,&FBCR,
                       &KNOPT, TEE, BCOEF, &IN, &KORDER, WORK2);
        for(int i=0; i<iheight;i++)
           Press[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv,heights+i,&inbv,WORK2);
         
        double ENN[iheight],Dens[iheight];  //Number of molecules per cubic metre, density Kg per m^3
        for(int i=0; i<iheight;i++){
            //multiply by 100 to convert pressure to Pascals
            ENN[i]=Press[i]/R/Temp[i]*AV*100.; 
            // multiply by MW to get mass per vol in gramms
            // divide by 1000 to get mass per vol in Kg
            // nucleon is just 1/Av 
            Dens[i]=ENN[i]*MWA*nucleon/1000.0;    
            }

        double GPPM[idata];   //coarse grid gas parts per million (spare copy needed for interpolation)
        double G[iheight];    //fine grid data gas ppm, spare copy
        double Humids1[iheight];   //fine grid data version of Hum1
        double Humids2[iheight];   //fine grid data version of Hum2
        double Gas[iheight][ngas];  //fine grid array of gas ppm for all gases not ignored

        for(int n=0;n<ngas;n++){
           for(int i=0;i<idata;i++)
                GPPM[i]=GasPPM[i][n];

           inbv=1; //inbv must always be set to 1 on first call to dbvalu
           dbint4_(altitude,GPPM,&idata,&IBCL,&IBCR,&FBCL,&FBCR,
                  &KNOPT, TEE, BCOEF,&IN, &KORDER, WORK);

           for(int i=0; i<iheight;i++){
                 G[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                          ,&ideriv,heights+i,&inbv,WORK2);
                 Gas[i][n]=G[i];
                       }
           if(n==0){
           //remember ihumid is now 1, so we now get
           //Humids1 andd Humids2 from temporary ppm profile in G.
           iwi=1;
           WaterVap(Press,Temp,G,Humids1,AV,nucleon,MWGAS[0],iheight,ihumid,iwi);
           iwi=2;
           WaterVap(Press,Temp,G,Humids2,AV,nucleon,MWGAS[0],iheight,ihumid,iwi);}
        }
        double DensGas[iheight][ngas];  // Density of each gas on fine grid.
        for( int n=0; n< ngas; n++){
           for( int i=0;i<iheight;i++){
               // divide ENN by 1E6 since Gas[i][n] is ppmv
               // divide by a 1000 and DensGas is in Kg m^{-2} (remember necleon is in grammes.
               DensGas[i][n]=ENN[i]/1.0E6*Gas[i][n]*MWGAS[n]*nucleon/1000.0;
           }
         }
/***************************************************************************************
                           END STAGE  2
****************************************************************************************
                           BEGIN STAGE 3
                           Do Integrals
****************************************************************************************/
            
       cout << "Begin Stage 3\n";
       int kountH=iheight; // count integration points top down down
       int kountS=0;       //kount slabs down
       int islab=up_to_top;
       double AvTemp[up_to_top]; double  AvPress[up_to_top];
       //mass of air (per square metre) in each layet and the total mass of all layers
       double massgas[up_to_top]; double totalmass=0.0;
       //mass of each constituent in each layer and the total mass
       double massgases[up_to_top][ngas]; double  totalmasses[ngas];

       int idown=0;
       for(int i=0; i<ngas; i++){
             totalmasses[i]=0.0;
       }
       cout << "***************************************\n";
       for(int i=nsplit1; i > 0; i--){
       cout << nsplit1-i<< "_______________________________________\n";
           AvTemp[nsplit1-i]=0.0;
           AvPress[nsplit1-i]=0.0;
           massgas[nsplit1-i]=0.0;
           for(int k=0; k<ngas; k++){
             massgases[nsplit1-i][k]=0.0;
           }
           for(int j=0; j < ngauss_height; j++){
             kountH--;
             cout << "j=" <<  j << "  " << Wgauss[j] << " " << Temp[kountH] << " " << kountH << endl;
             AvTemp[nsplit1-i]=AvTemp[nsplit1-i]+Wgauss[j]*Temp[kountH];
             AvPress[nsplit1-i]=AvPress[nsplit1-i]+Wgauss[j]*Press[kountH];
             massgas[nsplit1-i]=massgas[nsplit1-i]+Wgauss[j]*Dens[kountH];
             for(int k=0; k<ngas; k++){
               massgases[nsplit1-i][k]=massgases[nsplit1-i][k]+Wgauss[j]*DensGas[kountH][k];
             }
        }
        
        kountS++;
        islab=up_to_top-kountS;
        massgas[nsplit1-i]=massgas[nsplit1-i]*SlabThick[islab];
        for(int k=0; k<ngas; k++){
          massgases[nsplit1-i][k]=massgases[nsplit1-i][k]*SlabThick[islab];
          totalmasses[k]=totalmasses[k]+massgases[nsplit1-i][k];
        }
        totalmass=totalmass+massgas[nsplit1-i];
        if(verbose_height)cout << AvTemp[nsplit1-i] << "  " << AvPress[nsplit1-i] << "  " << totalmass <<endl;
        }
       idown=idown+nsplit1;
       cout << "***************************************\n";
       for(int i=nsplit2; i > 0; i--){
       cout << nsplit2-i+idown <<"_______________________________________\n";
           AvTemp[nsplit2-i+idown]=0.0;
           AvPress[nsplit2-i+idown]=0.0;
           massgas[nsplit2-i+idown]=0.0;
           for(int k=0; k<ngas; k++){
             massgases[nsplit2-i+idown][k]=0.0;
           }
           for(int j=0; j < ngauss_height; j++){
             kountH--;
             cout << "j=" << j << "  " << Wgauss[j] << " " << Temp[kountH] << " " << kountH << endl;
             AvTemp[nsplit2-i+idown]=AvTemp[nsplit2-i+idown]+Wgauss[j]*Temp[kountH];
             AvPress[nsplit2-i+idown]=AvPress[nsplit2-i+idown]+Wgauss[j]*Press[kountH];
             massgas[nsplit2-i+idown]=massgas[nsplit2-i+idown]+Wgauss[j]*Dens[kountH];
             for(int k=0; k<ngas; k++){
               massgases[nsplit2-i+idown][k]=massgases[nsplit2-i+idown][k]+Wgauss[j]*DensGas[kountH][k];
             }
        }
        kountS++;
        islab=up_to_top-kountS;
        massgas[nsplit2-i+idown]=massgas[nsplit2-i+idown]*SlabThick[islab];
        for(int k=0; k<ngas; k++){
          massgases[nsplit2-i+idown][k]=massgases[nsplit2-i+idown][k]*SlabThick[islab];
          totalmasses[k]=totalmasses[k]+massgases[nsplit2-i+idown][k];
        }
        totalmass=totalmass+massgas[nsplit2-i+idown];
        if(verbose_height)cout << AvTemp[nsplit2-i+idown] << "  " << AvPress[nsplit2-i+idown] << "  " << totalmass << endl;
        }
       idown=idown+nsplit2;
       cout << "***************************************\n";
       for(int i=nsplit3; i > 0; i--){
       cout << nsplit3-i+idown<< "_______________________________________\n";
           AvTemp[nsplit3-i+idown]=0.0;
           AvPress[nsplit3-i+idown]=0.0;
           massgas[nsplit3-i+idown]=0.0;
           for(int k=0; k<ngas; k++){
             massgases[nsplit3-i+idown][k]=0.0;
           }
           for(int j=0; j < ngauss_height; j++){
             kountH--;
             cout << "j=" << j << "  " << Wgauss[j] << " " << Temp[kountH] << " " << kountH << endl;
             AvTemp[nsplit3-i+idown]=AvTemp[nsplit3-i+idown]+Wgauss[j]*Temp[kountH];
             AvPress[nsplit3-i+idown]=AvPress[nsplit3-i+idown]+Wgauss[j]*Press[kountH];
             massgas[nsplit3-i+idown]=massgas[nsplit3-i+idown]+Wgauss[j]*Dens[kountH];
             for(int k=0; k<ngas; k++){
               massgases[nsplit3-i+idown][k]=massgases[nsplit3-i+idown][k]+Wgauss[j]*DensGas[kountH][k];
             }
        }
        kountS++;
        islab=up_to_top-kountS;
        massgas[nsplit3-i+idown]=massgas[nsplit3-i+idown]*SlabThick[islab];
        for(int k=0; k<ngas; k++){
          massgases[nsplit3-i+idown][k]=massgases[nsplit3-i+idown][k]*SlabThick[islab];
          totalmasses[k]=totalmasses[k]+massgases[nsplit3-i+idown][k];
        }
        totalmass=totalmass+massgas[nsplit3-i+idown];
        if(verbose_height)cout << AvTemp[nsplit3-i+idown] << "  " << AvPress[nsplit3-i+idown] << "  " << totalmass << endl;
        }
       idown=idown+nsplit3;
       cout << "***************************************\n";
       for(int i=nsplit4; i > 0; i--){
       cout << nsplit4-i+idown<< "_______________________________________\n";
           AvTemp[nsplit4-i+idown]=0.0;
           AvPress[nsplit4-i+idown]=0.0;
           massgas[nsplit4-i+idown]=0.0;
           for(int k=0; k<ngas; k++){
             massgases[nsplit4-i+idown][k]=0.0;
           }
           for(int j=0; j < ngauss_height; j++){
             kountH--;
             cout << "j=" << j << "  " << Wgauss[j] << " " << Temp[kountH] << " " << kountH << endl;
             AvTemp[nsplit4-i+idown]=AvTemp[nsplit4-i+idown]+Wgauss[j]*Temp[kountH];
             AvPress[nsplit4-i+idown]=AvPress[nsplit4-i+idown]+Wgauss[j]*Press[kountH];
             massgas[nsplit4-i+idown]=massgas[nsplit4-i+idown]+Wgauss[j]*Dens[kountH];
             for(int k=0; k<ngas; k++){
               massgases[nsplit4-i+idown][k]=massgases[nsplit4-i+idown][k]+Wgauss[j]*DensGas[kountH][k];
             }
        }
        kountS++;
        islab=up_to_top-kountS;
        massgas[nsplit4-i+idown]=massgas[nsplit4-i+idown]*SlabThick[islab];
        for(int k=0; k<ngas; k++){
          massgases[nsplit4-i+idown][k]=massgases[nsplit4-i+idown][k]*SlabThick[islab];
          totalmasses[k]=totalmasses[k]+massgases[nsplit4-i+idown][k];
        }
        totalmass=totalmass+massgas[nsplit4-i+idown];
        if(verbose_height)cout << AvTemp[i+idown] << "  " << AvPress[i+idown] << "  " << totalmass << endl;
        }

        cout << "mass of gas per square metre * g=" << 9.81*totalmass << "  Pascals="
              << 9.81*totalmass/100.0 << "mbar\n";
        cout << "Ground pressure was given as " << P[0] << "mbar:  Discrepancy due to interpolation error\n";
        cout << "We now scale masses to regain original ground pressure and convert to g per cm^2\n";
        cout << "amd calculate number densities per cm^2 to calculate line strengths\n";

        double correction=9.81*totalmass/P[0];

        for(int i=0; i< iheight; i++){
               Press[i]=Press[i]/correction; }
        for(int i=0; i<5; i++){
               PBDRY[i]=PBDRY[i]/correction; }
        for(int i=0; i<up_to_top; i++){
               PBDRY[i]=PBDRY[i]/correction; }
               

        //convert from Kg per square metre to g per square cm.
        // Kg to g times=1000; metre square to cm square /=10000

        for(int k=0; k<ngas; k++){
          for(int i=0; i<up_to_top; i++){
             massgases[i][k]=massgases[i][k]/10.0/correction;
          }
          totalmasses[k]=totalmasses[k]/10.0/correction;
        }
        totalmass=totalmass/10.0/correction;

        cout << "gases in each layer --- top down --- g per cm^2\n";
        for(int i=0; i<up_to_top; i++){ //top down actually
            for(int k=0; k < ngas ; k++){
                  cout << massgases[i][k] << " ";}
                  cout << endl;}
        cout << "totals\n";
        for(int k=0; k < ngas ; k++){
             cout << totalmasses[k] << " ";}

        cout <<  "  totalmass=" << totalmass << "g per cm^2" << "  ground pressure now " << 9.81*10*totalmass << endl; 

/***************************************************************************************
                          END STAGE 3
****************************************************************************************
                         BEGIN STAGE 4
                         HITRAN DATA 
****************************************************************************************/

//    declarations for HITRAN input
      cout <<" Now read HITRAN par files\n";    
      int Imol;                                    //molecule number
      int Iso;                                     //Isotopologue
      double nu;                                   //Wave number
      double S_intense;                            //line intensity
      double Einstein_A;                           //Einstein A coefficient
      double gamma_air, gamma_self;                //broadened half-widths
      double LSE;                                  //lower state energy
      double Tdep_gamma_air;                       //For temperature dependence gammma_air
      double Press_shift   ;                       //Shift due to pressure
      int ierr[6];                                 //6 error parameters
      int iref[6];                                 //6 reference parameters
      char flag;                                   //flag for possibility of line mixing
      double gp;                                   //Stat. Weight Upper 
      double gpp;                                  //Stat. Weight Lower

      //wave number interval based on lambda1 and lambda2
      //look at wavenumbers nu1-cutoff to nu2+cutoff;
      double nu1, nu2, cutoff, lambda1X,lambda2X;
      nu1=1e4/lambda2; nu2=1e4/lambda1; cutoff=0.1;
      nu1=nu1-nu1*cutoff; nu2=nu2+nu2*cutoff;
      lambda1X=1.0e4/nu2; lambda2X=1.0e4/nu1;
     
      cout << "nu1 nu2 " << nu1 << "  " << nu2  <<endl;
      cout << "lambda1X  lambda2X  " << lambda1X << "  " << lambda2X << endl;

      string  PARFILE="HITRAN/   hit08.par";
      double lambda;   // a wavelength
      
      string line;
      string entry[19];
      istringstream iss[19]; //19 input stringstreams
      for(int ig=0; ig<ngas; ig++){ //loop to read from HITRAN files.
        int molecule=gases[ig];
        int ireplace1=molecule/10; int  ireplace2=molecule%10;
        ireplace1=ireplace1+48; ireplace2=ireplace2+48; //ascci '0' is 48
        ostringstream replacewith;
        replacewith << (char)ireplace1 << (char)ireplace2 << '_';
        PARFILE.replace(7,3,replacewith.str());
        cout <<"PARFILE=" << PARFILE << endl;
        fp_in.open(PARFILE.c_str(),ios::in);
        if(!fp_in.is_open()){cout <<" can't open HITRAN file \n";
                             cout << PARFILE.c_str() << endl;  exit(1);}
        int kountlines=0;
        int vecsize=nice[ig];
        if(vecsize>nicemax)vecsize=nicemax;

        vector<double> linecentres[vecsize];
        vector<double> linestrengths[vecsize];
        vector<double> linegammaAir[vecsize];     // air broadening
        vector<double> linegammaSelf[vecsize];    // self broadening broadening
        vector<double> lineLSE[vecsize];          // Lower state Energy
        vector<double> line_gp[vecsize];          // Stat weight upper
        vector<double> line_gpp[vecsize];         // Stat weight lower     

        for(int k=0; k<vecsize; k++){
            linecentres[k].clear();
            linestrengths[k].clear();
            linegammaAir[k].clear();
            linegammaSelf[k].clear();
            lineLSE[k].clear();    
            line_gp[k].clear();    
            line_gpp[k].clear();            
         }

        while(!getline(fp_in, line).eof()){
          entry[1]=line.substr(2,1); //
          iss[1].str(entry[1]);
          iss[1] >> Iso;


          if(Iso<=nicemax){ // number of isotopologues considered
          entry[2]=line.substr(3,12); //get wave number first
          iss[2].str(entry[2]);
          iss[2] >> nu;
          if(nu > nu1 && nu < nu2){
          cout << line << endl;
          kountlines++;
          entry[0]=line.substr(0,2); 
          entry[3]=line.substr(15,10); entry[4]=line.substr(25,10); entry[5]=line.substr(35,5);
          entry[6]=line.substr(40,5); entry[7]=line.substr(45,10); entry[8]=line.substr(55,4);
          entry[9]=line.substr(59,8); entry[10]=line.substr(67,15); entry[11]=line.substr(82,15);
          entry[12]=line.substr(97,15); entry[13]=line.substr(112,15); entry[14]=line.substr(127,6);
          entry[15]=line.substr(133,12); entry[16]=line.substr(145,1); entry[17]=line.substr(146,7);
          entry[18]=line.substr(153,7);

          iss[0].str(entry[0]);  
          iss[0] >> Imol;
          cout << "Imol="  << Imol << endl;
          cout << "Iso=" <<Iso << "  " << nicecode[Imol-1][Iso-1]  <<endl;
          cout << "nu="  << nu <<  "  "  << nu1 << "  " << nu2 << endl;
          cout << 10000/nu  << " " << lambda1X << " " << lambda2X<<  endl;

          linecentres[Iso-1].push_back(nu);
          linestrengths[Iso-1].push_back(S_intense);
          linegammaAir[Iso-1].push_back(gamma_air);
          linegammaSelf[Iso-1].push_back(gamma_self);
          lineLSE[Iso-1].push_back(LSE);
          line_gp[Iso-1].push_back(gp);
          line_gpp[Iso-1].push_back(gpp);

          iss[3].str(entry[3]);
          iss[3] >> S_intense;
           cout << "S_intense="  << S_intense<< endl;
          iss[4].str(entry[4]);
          iss[4] >> Einstein_A;
          cout << "Einstein_A="  << Einstein_A<< endl;
          iss[5].str(entry[5]);
          iss[5] >> gamma_air;
          cout << "gamma_air="  << gamma_air<< endl;
          iss[6].str(entry[6]);
          iss[6] >> gamma_self;
          cout << "gamma_self="  << gamma_self<< endl;
          iss[7].str(entry[7]);
          iss[7] >> LSE;
          cout << "LSE="  << LSE << endl;
          iss[8].str(entry[8]);
          iss[8] >> Tdep_gamma_air;
          cout << "Tdep_gamma_air="  << Tdep_gamma_air<< endl;
          iss[9].str(entry[9]);
          iss[9] >> Press_shift;
          cout << "Press_shift="  << Press_shift << endl;
          // entry[10]=Upper state global parse according to classes 1 to 10 Table 3 HITRAN 2004
          // entry[11]=Lower state global parse according to classes 1 to 10 Table 3 HITRAN 2004
          // entry[12]=Upper state local  parse according to groups 1 to 6 Table 4 HITRAN 2004
          // entry[13]=Lower state local  parse according to groups 1 to 6 Table 4 HITRAN 2004
          // not used in this code
          // skip to entry[14]
          //ascii character '0'=48  'space'=32
          ierr[0]=(int)entry[14][0]-48; ierr[1]=(int)entry[14][1]-48;
          ierr[2]=(int)entry[14][2]-48; ierr[3]=(int)entry[14][3]-48;
          ierr[4]=(int)entry[14][4]-48; ierr[5]=(int)entry[14][5]-48;
       
 //         cout << "ierr[6]="  << ierr[0] << " " << ierr[1] << " "<<  ierr[2] << " "
  //            << ierr[3] << " " << ierr[4] << " "  <<  ierr[5] << endl;
          int idum1,idum2;
          idum1=(int)entry[15][0]; idum2=(int)entry[15][1];
          if(idum1==32){idum1=0;} else {idum1=idum1-48;}  if(idum2==32){idum2=0;} else {idum2=idum2-48;}
          iref[0]=10*idum1+idum2;    
          idum1=(int)entry[15][2]; idum2=(int)entry[15][3];
          if(idum1==32){idum1=0;} else {idum1=idum1-48;}  if(idum2==32){idum2=0;} else {idum2=idum2-48;}
          iref[1]=10*idum1+idum2;    
          idum1=(int)entry[15][4]; idum2=(int)entry[15][5];
          if(idum1==32){idum1=0;} else {idum1=idum1-48;}  if(idum2==32){idum2=0;} else {idum2=idum2-48;}
          iref[2]=10*idum1+idum2;    
          idum1=(int)entry[15][6]; idum2=(int)entry[15][7];
          if(idum1==32){idum1=0;} else {idum1=idum1-48;}
          if(idum2==32){idum2=0;} else {idum2=idum2-48;}
          iref[3]=10*idum1+idum2;    
          idum1=(int)entry[15][8]; idum2=(int)entry[15][9];
          if(idum1==32){idum1=0;} else {idum1=idum1-48;}
          if(idum2==32){idum2=0;} else {idum2=idum2-48;}
          iref[4]=10*idum1+idum2;    
          idum1=(int)entry[15][10]; idum2=(int)entry[15][11];
          if(idum1==32){idum1=0;} else {idum1=idum1-48;}
          if(idum2==32){idum2=0;} else {idum2=idum2-48;}
          iref[5]=10*idum1+idum2;    
  //        cout << "iref[6]="<< iref[0] << " " << iref[1] << " " << iref[2] <<  " "
   //         << iref[3] << " " << iref[4] << " " << iref[5] << endl;
          flag=(char)entry[16][0];
          iss[17].str(entry[17]); iss[17] >> gp;   iss[18].str(entry[18]); iss[18] >> gpp;
   //       cout << "gp and gpp (stat weights upper and lower) are " << gp << " and "<< gpp << endl;
   //
          iss[0].clear();iss[0].clear();iss[3].clear();iss[4].clear();
          iss[0].clear();iss[5].clear();iss[6].clear();iss[7].clear();iss[8].clear();iss[9].clear();
          iss[10].clear();iss[11].clear();iss[12].clear();iss[13].clear();iss[14].clear();iss[15].clear();
          iss[16].clear();iss[17].clear();iss[18].clear();
          }  //endif for lambda within range
          iss[2].clear();
          } //endif for isotope less than nicemax
          iss[1].clear();

          } //end while loop
        fp_in.close();
        //At last we have the line centres for this particular gas.
        
/*****************************************************************************************
        Calculate Spectrum - Work out Distributions
******************************************************************************************/


/*****************************************************************************************
       Write Input Files for Dors3
******************************************************************************************/

/****************************************************************************************/
        for(int i=0;i<vecsize;i++){
            cout << "gas humber=" << ig+1 << " Isotopologue=" << i+1 << "  code=" <<nicecode[ig][i]
                 << "  number of lines=" << linecentres[i].size()<< endl;
                 linecentres[i].erase(linecentres[i].begin(), linecentres[i].end()); 
                 linestrengths[i].erase(linestrengths[i].begin(), linestrengths[i].end()); 
                 linegammaAir[i].erase(linestrengths[i].begin(), linestrengths[i].end()); 
                 linegammaSelf[i].erase(linestrengths[i].begin(), linestrengths[i].end());  
                 lineLSE[i].erase(linestrengths[i].begin(), linestrengths[i].end());       
                 line_gp[i].erase(linestrengths[i].begin(), linestrengths[i].end());     
                 line_gpp[i].erase(linestrengths[i].begin(), linestrengths[i].end());            
            }
      }
      //end loop over HITRAN par files ig=0 to ngas


  return 0;
}


