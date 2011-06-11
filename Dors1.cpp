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

//HITRAN Partition Function  --- note calling from c++ means all lower case
void bd_tips_2003_(int&, double&, int&, double&, double&);

void humlik_(int&, double*, double&, double*);
}

void Questioner( bool&,
                int&, int&, int&, int&, bool&, double&,
                double&, double&, double&, bool&,  int&, bool&, bool&, int&,
                int&, int&, int&, bool&, bool&, bool&, 
                int&, double&, double&,
                bool&, double&, double&, double&, double&, double&, int&, 
                bool&, bool &, int&, int&, int&, int&,
                const char*);

void WaterVap(double* , double* , double* ,
              double* , double &, double &, double &, int &, int &, int &);

              

int main(int argc, char* argv[]){

  bool verbose;   //verbose mode for input, read from file if false;
  //verbose output bools
  bool verbose_out=true;  // output loads of info to standard output
  bool verbose_height=false;   //output height integration info to std::out
  bool verbose_spect=false;   // output spectrum details to std out
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
  double nu_cut; int icutL, icutD, istep;   //wave number cutoffs and spectral resolution
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
   //  verbose=true;   // comment out if using ddd or other display debugger

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
              ihumid, groundtemp, groundpress,
              default_pause, HG, HB, HT, HS, HU, ngauss_height, 
              calcspec, outspec, ngauss_correlkay, icutL,icutD,istep,
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
   /*        Quadrature Rule for Correlated kay                   */
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

   //constants for spectral calculations using HITRAN data
   double RefTemp=296.0;   //Temperature at which linestrength data are taken
   double RefPress=1013*100.0; //Converts mbar to Pascals    
   double C2=1.4387752;   //=hc/Boltzman constant=second radiation constant. N.B.  units are Kelvin per cm
                          //so C2*wavenumber/T per cm is (photon energy)/kT
   double log2=log(2.0);
   double sqrtlog2=sqrt(log2);
   double sqrtlog2opi=sqrt(log2/pi);
   double Boltz=1.3806503e-16;  //Boltzmann's constant (cgs units)
   //speed of light in cm per sec.
   double Speedlight=2.99792458e10;



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
   if(verbose_spect)cout << line1 << endl;
   int  imol=0;
   while(!getline(fp_in, line1).eof()){
      if(verbose_spect)cout << line1 << endl;
      for(int iso=0; iso<nice[imol];iso++){
         fp_in >>  nicecode[imol][iso] >> niceabund[imol][iso] >> QTref[imol][iso]
               >> geejay[imol][iso] >> AllMolW[imol][iso];
         if(verbose_spect)cout << imol+1 << " " << nice[imol] << "  " << nicecode[imol][iso] << "  " << AllMolW[imol][iso] << endl;
      }
      if(!getline(fp_in,line1)){ break;}
       // << "  getline return\n";
      if(verbose_spect)cout << line1;
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
           if(verbose_height)cout <<"kount=" << kount << "  heights=" << heights[kount] << endl;
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
            ENN[i]=Press[i]*100.0/R/Temp[i]*AV; 
            // multiply by MW to get mass per vol in grammes
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
       if(verbose_out)cout << "***************************************\n";
       for(int i=nsplit1; i > 0; i--){
           if(verbose_out)cout << nsplit1-i<< "_______________________________________\n";
           AvTemp[nsplit1-i]=0.0;
           AvPress[nsplit1-i]=0.0;
           massgas[nsplit1-i]=0.0;
           for(int k=0; k<ngas; k++){
             massgases[nsplit1-i][k]=0.0;
           }
           for(int j=0; j < ngauss_height; j++){
             kountH--;
             if(verbose_out)cout << "j=" <<  j << "  " << Wgauss[j] << " " << Temp[kountH] << " " << kountH << endl;
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
       if(verbose_out)cout << "***************************************\n";
       for(int i=nsplit2; i > 0; i--){
           if(verbose_out)cout << nsplit2-i+idown <<"_______________________________________\n";
           AvTemp[nsplit2-i+idown]=0.0;
           AvPress[nsplit2-i+idown]=0.0;
           massgas[nsplit2-i+idown]=0.0;
           for(int k=0; k<ngas; k++){
             massgases[nsplit2-i+idown][k]=0.0;
           }
           for(int j=0; j < ngauss_height; j++){
             kountH--;
             if(verbose_out)cout << "j=" << j << "  " << Wgauss[j] << " " << Temp[kountH] << " " << kountH << endl;
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
       if(verbose_out)cout << "***************************************\n";
       for(int i=nsplit3; i > 0; i--){
          if(verbose_out)cout << nsplit3-i+idown<< "_______________________________________\n";
           AvTemp[nsplit3-i+idown]=0.0;
           AvPress[nsplit3-i+idown]=0.0;
           massgas[nsplit3-i+idown]=0.0;
           for(int k=0; k<ngas; k++){
             massgases[nsplit3-i+idown][k]=0.0;
           }
           for(int j=0; j < ngauss_height; j++){
             kountH--;
             if(verbose_out)cout << "j=" << j << "  " << Wgauss[j] << " " << Temp[kountH] << " " << kountH << endl;
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
       if(verbose_out)cout << "***************************************\n";
       for(int i=nsplit4; i > 0; i--){
           if(verbose_out)cout << nsplit4-i+idown<< "_______________________________________\n";
           AvTemp[nsplit4-i+idown]=0.0;
           AvPress[nsplit4-i+idown]=0.0;
           massgas[nsplit4-i+idown]=0.0;
           for(int k=0; k<ngas; k++){
             massgases[nsplit4-i+idown][k]=0.0;
           }
           for(int j=0; j < ngauss_height; j++){
             kountH--;
             if(verbose_out)cout << "j=" << j << "  " << Wgauss[j] << " " << Temp[kountH] << " " << kountH << endl;
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
              << "  " <<9.81*totalmass/100.0 << "mbar\n";
        cout << "Ground pressure was given as " << P[0] << "mbar:  Discrepancy due to interpolation+integration error\n";
        cout << "We now scale masses to regain original ground pressure and convert to g per cm^2\n";
        cout << "and calculate number densities per cm^2 to calculate line strengths\n";

        //convert P[0] mbar to Pascals,  1mbar=100  Pascals - otherwise masses are all wrong!
       
        double correction=9.81*totalmass/(P[0]*100);   

        for(int i=0; i< iheight; i++){
               Press[i]=Press[i]/correction; }
        for(int i=0; i<5; i++){
               PBDRY[i]=PBDRY[i]/correction; }
        for(int i=0; i<up_to_top; i++){
               PBDRY[i]=PBDRY[i]/correction; }
               

        //convert from Kg per square metre to g per square cm.
        // Kg to g times=1000; metre square to cm square /=10000
        // g  per square cm = Kg per square metre /1000*10000 
        // We have calculated Kg per square metre---divide by 10 to get g per cm^2.
        //

        for(int k=0; k<ngas; k++){
          for(int i=0; i<up_to_top; i++){
             massgases[i][k]=massgases[i][k]/10.0/correction;
          }
          totalmasses[k]=totalmasses[k]/10.0/correction;
        }
        totalmass=totalmass/10.0/correction;

        if(verbose_out)cout << "gases in each layer --- top down --- g per cm^2\n";
        for(int i=0; i<up_to_top; i++){ //top down actually
            for(int k=0; k < ngas ; k++){
                  if(verbose_out)cout << massgases[i][k] << " ";}
                  if(verbose_out)cout << endl;}
        if(verbose_out)cout << "totals\n";
        for(int k=0; k < ngas ; k++){
             if(verbose_out)cout << totalmasses[k] << " ";}

        if(verbose_out)cout <<  "  totalmass=" << totalmass << "g per cm^2" << "  ground pressure now "
                            << 9.81*10*totalmass << " Pascals" <<endl; 

/***************************************************************************************
                          END STAGE 3
****************************************************************************************
                         BEGIN STAGE 4
                         HITRAN DATA 
****************************************************************************************/
      if(calcspec){    //there are reasons we might not want to bother!
//    declarations for HITRAN input
      cout <<" Now read HITRAN par files\n";    
      int Imol;                                    //molecule number
      int Iso;                                     //Isotopologue
      double nu=0.0;                                   //Wave number
      double S_intense;                            //line intensity
      double Einstein_A;                           //Einstein A coefficient
      double gamma_air, gamma_self;                //broadened half-widths
      double LSE;                                  //lower state energy
      double T_depLorentz;                       //For temperature dependence gammma_air
      double Press_shift   ;                       //Shift due to pressure
      int ierr[6];                                 //6 error parameters
      int iref[6];                                 //6 reference parameters
      char flag;                                   //flag for possibility of line mixing
      double gp;                                   //Stat. Weight Upper 
      double gpp;                                  //Stat. Weight Lower

      //wave number interval based on lambda1 and lambda2
      //look at wavenumbers nu1-cutoff to nu2+cutoff;
      double nu1, nu2, nu1X, nu2X, lambda1X, lambda2X;
      nu1=1e4/lambda2; nu2=1e4/lambda1;
      nu1X=nu1-nu_cut; nu2X=nu2+nu_cut;
      lambda1X=1.0e4/nu2X; lambda2X=1.0e4/nu1X;
     
      cout << "nu1X nu2X " << nu1X << "  " << nu2X  <<endl;
      cout << "lambda1X  lambda2X  " << lambda1X << "  " << lambda2X << endl;
      bool singlewave=false;

      if(lambda2<lambda1){
          cout << "second wavelength smaller than first, Dors1 Failed \n"; exit(1);}

      if(lambda2 > lambda1){
          singlewave=false;
      }
      else{singlewave==true;} 

      int n_waves;  // number of wave numbers for fine detail spectrum
      double* wavearray;  // array for fine detail wave numbers
      
      
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
        int vecsize=nice[ig];
        if(vecsize>nicemax)vecsize=nicemax;
      
        // for picking out max line strength for each Isotopologue
        int istr_max[vecsize];    // which line has strength str_max?
        int ikount[vecsize];      // number of lines for isotopolgue
        double str_max[vecsize];  // max line strength for Isotopolgue
        for(int i=0; i<vecsize; i++){
           istr_max[i]=0;
           ikount[i]=-1;
           str_max[i]=-1.0;
        }
        //vectors for each isotopologue, vecsize <=nicemax depending on nice[ig]
        vector<double> linecentres[vecsize];
        vector<double> linestrengths[vecsize];
        vector<double> linegammaAir[vecsize];     // air broadening
        vector<double> linegammaSelf[vecsize];    // self broadening broadening
        vector<double> lineLSE[vecsize];          // Lower state Energy
        vector<double> lineTdep[vecsize];          // Stat weight upper
        vector<double> lineShift[vecsize];         // Stat weight lower     

        for(int k=0; k<vecsize; k++){
            linecentres[k].clear();
            linestrengths[k].clear();
            linegammaAir[k].clear();
            linegammaSelf[k].clear();
            lineLSE[k].clear();    
            lineTdep[k].clear();    
            lineShift[k].clear();            
         }

        // Loop that reads in HITRAN data
        while(!getline(fp_in, line).eof()){  // get the line "line" --- what a jazzy name!
          entry[1]=line.substr(2,1); // first of 19 substrings called entry  (start pos=2, length=1)
          iss[1].str(entry[1]);      // input to first of 19 input string streams called iss
          iss[1] >> Iso;             // put content of second iss to Iso (iss[0] is lower down)


          if(Iso<=nicemax){ // number of isotopologues considered
          entry[2]=line.substr(3,11); //get wave number first
          iss[2].str(entry[2]);
          iss[2] >> fixed;
          // setprecision fixes the (maximum) number of digits
          // but fortran F12.6 is the width of the field and 6 digits after the decimal
          // width of field includes the decimal point -hence setprecision with 11.
          nu=0.0;
          iss[2] >> setprecision(11) >> nu;
          if(nu > nu1X && nu < nu2X){
          if(verbose_spect)cout << line << endl;
          entry[0]=line.substr(0,2); 
          entry[3]=line.substr(15,10); entry[4]=line.substr(25,10); entry[5]=line.substr(35,5);
          entry[6]=line.substr(40,5); entry[7]=line.substr(45,10); entry[8]=line.substr(55,4);
          entry[9]=line.substr(59,8); entry[10]=line.substr(67,15); entry[11]=line.substr(82,15);
          entry[12]=line.substr(97,15); entry[13]=line.substr(112,15); entry[14]=line.substr(127,6);
          entry[15]=line.substr(133,12); entry[16]=line.substr(145,1); entry[17]=line.substr(146,7);
          entry[18]=line.substr(153,7);

          iss[0].str(entry[0]);  
          iss[0] >> Imol;
          if(verbose_spect)cout << "Imol="  << Imol << endl;
          if(verbose_spect)cout << "Iso=" <<Iso << "  " << nicecode[Imol-1][Iso-1]  <<endl;
          if(verbose_spect)cout << "nu="  << nu <<  "  "  << nu1 << "  " << nu2 << endl;
          if(verbose_spect)cout << 10000/nu  << " " << lambda1X << " " << lambda2X<<  endl;

          linecentres[Iso-1].push_back(nu);
          
          iss[3].str(entry[3]);
          iss[3] >> S_intense;
          linestrengths[Iso-1].push_back(S_intense);
          if(verbose_spect)cout << "S_intense="  << S_intense<< endl;

          // pick out the maximum line strength for each isotopologue
          // and its index
          ikount[Iso-1]++;
          if(linestrengths[Iso-1][ikount[Iso-1]]>str_max[Iso-1]){
                    str_max[Iso-1]=linestrengths[Iso-1][ikount[Iso-1]];
                    istr_max[Iso-1]=ikount[Iso-1];
          }

          iss[4].str(entry[4]);
          iss[4] >> Einstein_A;
          if(verbose_spect)cout << "Einstein_A="  << Einstein_A<< endl; //not used for local thermodynamic equilibrium

          iss[5].str(entry[5]);
          iss[5] >> gamma_air;
          linegammaAir[Iso-1].push_back(gamma_air);
          if(verbose_spect)cout << "gamma_air="  << gamma_air<< endl;

          iss[6].str(entry[6]);
          iss[6] >> gamma_self;
          linegammaSelf[Iso-1].push_back(gamma_self);
          if(verbose_spect)cout << "gamma_self="  << gamma_self<< endl;

          iss[7].str(entry[7]);
          iss[7] >> LSE;
          lineLSE[Iso-1].push_back(LSE);
          if(verbose_spect)cout << "LSE="  << LSE << endl;

          iss[8].str(entry[8]);
          iss[8] >> T_depLorentz;
          lineTdep[Iso-1].push_back(T_depLorentz);
          if(verbose_spect)cout << "T_depLorentz="  << T_depLorentz<< endl;

          iss[9].str(entry[9]);
          iss[9] >> Press_shift;
          lineShift[Iso-1].push_back(Press_shift);
          if(verbose_spect)cout << "Press_shift="  << Press_shift << endl;

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

          } //end  if not end of file while loop
        fp_in.close();
        //Finished reading in HITRAN data
        //At last we have the line centres for this particular gas.

        cout << "Finished with PARFILE=" << PARFILE << endl;
        



       for(int k=0; k<vecsize; k++){   //loop over isotopologues

        int ig_fort=ig+1;  //gas and isotopoluge for BD_TIPS_2003.f
        int k_fort=k+1;

        int nlines=linecentres[k].size();

     // we can  examine only the spectrum from nu1X to nu2X
     // remember only line centres in this interval were pushed onto the linecentres[k]  vector!
     //   
     // Spectra for each layer top down --- see up_to_top and kountS
        double TempNlay,PressNlay;  //temperature of layerm Nlay
        double gee_i=0.0, QT=0.0;  // state independent degeneracy factor, Total internal partition sum
        double Refgee_i;
         //The line strengths are for a reference temperature of 296K
         //BD_TIPS_2003 Routine provided by HITRAN.

        bd_tips_2003_(ig_fort, RefTemp, k_fort, Refgee_i, QT);

        for(int nlay=0; nlay < up_to_top; nlay++){
            TempNlay=AvTemp[nlay];
            PressNlay=AvPress[nlay];
          //BD_TIPS_2003 Routine provided by HITRAN.

            bd_tips_2003_(ig_fort, TempNlay, k_fort, gee_i, QT);
           // cout << "gi=" << gee_i << "  QT=" << QT <<  endl;
 
            double wavenumlines[nlines];    // temporary copies
            double strengthlines[nlines];
            double knaughtlines[nlines];
            double alphaLorentzA[nlines];  //for air -add partial self and partial pressures later;
            double alphaDop[nlines];
            double widthlines[nlines];
 
            

            //First we find all the wave numbers at which to calculate
            //the entire spectrum over nu1X to nu2X and stuff them in SpectrumWaves;
            int kountwaves=0;  //wavelengths which are line centres
            for(int l=0; l<nlines;l++){
              if(QT>0){
              wavenumlines[l]=linecentres[k][l];
              //Voigt function needs half widths-not 1/e-hence log2
              alphaDop[l]=(wavenumlines[l]-lineShift[k][l]*AvPress[nlay]/RefPress)
                          /Speedlight*sqrt(2.0*Boltz*TempNlay*log2/(AllMolW[ig][k]*nucleon));
              alphaLorentzA[l]=linegammaAir[k][l]*AvPress[nlay]/RefPress*pow(RefTemp/TempNlay,lineTdep[k][l]);
              strengthlines[l]=linestrengths[k][l]*QTref[ig][k]/QT
                      *exp(-C2*lineLSE[k][l]/TempNlay)/exp(-C2*lineLSE[k][l]/RefTemp)
                      *( 1.0-exp(-C2*linecentres[k][l]/TempNlay) )*( 1.0-exp(-C2*linecentres[k][l]/RefTemp) );
              knaughtlines[l]=strengthlines[l]/alphaDop[l]*sqrtlog2opi;
              }
              else{
                 cout << "Error: BD_TIPS_2003 returned QT=" << QT << endl;
                 exit(1);
              }
              //now we have the x and y vectors for all the lines
            }  //end loop over nlines


            vector<double> SpectrumWaves; vector<double>SpectrumKays;
            vector<int> CentreWaves;  //which values of nu in SpectrumWaves are at line centres?
            vector<int> CentreLines;
            SpectrumWaves.clear();
            CentreWaves.clear();
            CentreLines.clear();
            double startnu,stepnu,tempnu,currentnu0,
            currentwidth,currentwidthL,currentwidthD,nextnu0,last_stopnu, stopnu, Wtol;


            Wtol=0.000001;  //because of fortran F12.6  for wave number in par file

            int icurrent=0;
            last_stopnu=nu1X;
            bool skipline=false;
            double midpoint;

            while(icurrent < nlines){  //we get some "double lines" -skip the next line if it it's the same nu0
              if(skipline){
                 skipline=false;
                 icurrent++;
                 //we are just calculating where on the nu axis our fine detail spectrum is calculated
                 //the two transitions with the same frequency have different widths and strengths
                 //so it must appear in the CentreWaves list. It will be the last centre wave number.
                 CentreWaves.push_back(CentreWaves[CentreWaves.size()-1]);
                 CentreLines.push_back(CentreLines[CentreLines.size()-1]);
              }
              else{
                currentnu0=wavenumlines[icurrent];
                if(icurrent<nlines-1){               
                  nextnu0=wavenumlines[icurrent+1];
                   if(currentnu0==nextnu0){
                     skipline=true;
                 }
              }
              
              currentnu0=currentnu0-lineShift[k][icurrent]*AvPress[nlay]/RefPress;
              nextnu0=nextnu0-lineShift[k][icurrent]*AvPress[nlay]/RefPress;
              //set current width
              currentwidthD=alphaDop[icurrent]*icutD;
              currentwidthL=alphaLorentzA[icurrent]*icutL;
              if(currentwidthL<currentwidthD){currentwidth=currentwidthD;}
              else{currentwidth=currentwidthL;}
          
              stepnu=currentwidth/istep;
              nu_cut=currentwidth;

              widthlines[icurrent]=currentwidth;

              if(currentnu0-nu_cut<last_stopnu)
                   {
                   startnu=last_stopnu+stepnu;
                   }

              else {startnu=currentnu0-nu_cut;}

              tempnu=startnu;

              while(tempnu<currentnu0-Wtol){   //step up to left of line which is first stopping point
                SpectrumWaves.push_back(tempnu);
                tempnu+=stepnu; kountwaves++;
              }

              //we have reached the line centre. The last +=stepnu in the while
              //loop may have just stepped over so we don't actually want it. 
              if(tempnu<currentnu0-Wtol){
                   SpectrumWaves.push_back(currentnu0);
                   CentreWaves.push_back(kountwaves);
                   CentreLines.push_back(icurrent);
                   kountwaves++;}
              else{
                   SpectrumWaves.pop_back();
                   kountwaves--;
                   SpectrumWaves.push_back(currentnu0);
                   CentreWaves.push_back(kountwaves);
                   CentreLines.push_back(icurrent);
                   kountwaves++;}

              midpoint=(currentnu0+nextnu0)/2.0;
              if(icurrent==nlines-1){               
                  stopnu=currentnu0+nu_cut;
              }//last line -- nothing to the right
              else{
                if(midpoint > currentnu0+nu_cut){
                  stopnu=currentnu0+nu_cut;}
                  else{
                  stopnu=midpoint;}
              }

              tempnu=currentnu0; //start at line centre and step up
              bool past_it=false;

              while(tempnu<stopnu-Wtol && tempnu < nextnu0-Wtol){

               tempnu+=stepnu;
               if(tempnu>=nextnu0)past_it=true;
               if(!past_it){
                  SpectrumWaves.push_back(tempnu); 
                  kountwaves++;}
 
                        int idebug=kountwaves;
                        if(idebug==26395){
                             cout << "debug point\n";}
           }  //got to next half way point or end of spectrum or lines so close
              //that half way point+stepnu is past nextnu0
           if(past_it){
                  tempnu=tempnu-stepnu;
                  past_it=false;
                  SpectrumWaves.push_back(stopnu); 
                  kountwaves++;
           }
           last_stopnu=tempnu;
           icurrent++;
              } //end of if else for skipline
           }  //end of while icurrent < nlines


           int templength=SpectrumWaves.size();
           double myzero=0.0;
           for(int l=0; l<templength-1; l++){
              SpectrumKays.push_back(myzero);
           //chop this next bit of loop out once all the above thoroughly tested
              if(l==0)
              cout << "start of " << templength << " wave numbers--wave number=" << SpectrumWaves[l] << endl;
              if(SpectrumWaves[l+1] <= SpectrumWaves[l]){
                    cout << "Error at l=" << l << "  layer=" << nlay << "  Iso=" << k << endl;
                    exit(0);
              }
              if(l==templength-1)
              cout << "end of " << templength << " wave numbers--wave number=" << SpectrumWaves[l] << endl;
           }
           SpectrumKays.push_back(myzero);        

           double xfac, yfac, xtemp, width, alphaD, alphaL;
           int leftwaves, rightwaves, iwave, SpectrumSize;
           int nulength=CentreWaves.size();
           int speclength=SpectrumWaves.size();
           
           for(int l=0; l<nulength-1; l++){
              currentnu0=wavenumlines[l]-lineShift[k][icurrent]*AvPress[nlay]/RefPress;
              width=widthlines[l];   //cutoff from previous loop
              alphaD=alphaDop[l];
              alphaL=alphaLorentzA[l];
              xfac=log2/alphaD;
              yfac=alphaL/alphaD*log2; 
              //Need to calculate Xvector for HUMLIK
              //But it needs to be ordered. We will have a normal array (ordered)
              // and a C++ vector (disordered)
              vector <double> Xvector;
              vector <int>  ISpect;
              Xvector.clear();
              ISpect.clear();
              leftwaves=0; rightwaves=0; iwave=CentreWaves[l];

              //linecentre then everything within range on left
              xtemp=SpectrumWaves[iwave-leftwaves];
           //   cout << setprecision(11) << xtemp << endl;
              Xvector.push_back(xtemp);
              ISpect.push_back(iwave-leftwaves);
              leftwaves++;              

              while(xtemp>currentnu0-width && iwave-leftwaves>0){
                   xtemp=SpectrumWaves[iwave-leftwaves];
                   // cout << setprecision(11) << xtemp << endl;
                   Xvector.push_back(xtemp); 
                   ISpect.push_back(iwave-leftwaves);
                   leftwaves++;
              } 

              rightwaves=1;
              //linecentre then everything within range on left

              if(iwave+rightwaves<speclength-1){
              xtemp=SpectrumWaves[iwave+rightwaves];
              Xvector.push_back(xtemp); 
              ISpect.push_back(iwave+rightwaves);
              rightwaves++;
              while(xtemp<currentnu0+width && iwave+rightwaves<speclength){
                   xtemp=SpectrumWaves[iwave+rightwaves];
                  // cout << setprecision(11) << xtemp << endl;
                   Xvector.push_back(xtemp); 
                   ISpect.push_back(iwave+rightwaves);
                   rightwaves++;
              } 
              } //endif for going past spectrum limit
              // the rightwaves block looks like a 1 to N vector, not a 0 to N-1 so shave 1 off rightwaves
              rightwaves--;
              int Nwaves=leftwaves+rightwaves;

              double XVEC[Nwaves],KVEC[Nwaves];
              for(int j=leftwaves-1; j>=0; j--){
                  XVEC[leftwaves-1-j]=(Xvector[j]-Xvector[leftwaves-1])*xfac;
                  KVEC[leftwaves-1-j]=0.0;
              }
             for(int j=leftwaves; j< Nwaves; j++){
                  XVEC[j]=(Xvector[j]-Xvector[leftwaves-1])*xfac;
                  KVEC[j]=0.0;
              }
              //calculate voigt profile
              humlik_(Nwaves, XVEC, yfac, KVEC);

              for(int j=leftwaves-1; j>=0; j--){
                  SpectrumKays[ISpect[leftwaves-1-j]]=SpectrumKays[ISpect[leftwaves-1-j]]
                                                     +KVEC[leftwaves]*knaughtlines[l];
                  cout << setprecision(11) << XVEC[leftwaves-1-j] << " K= " << KVEC[leftwaves-1-j] << " I="
                       << leftwaves-1-j << "  J=" << ISpect[leftwaves-1-j] << endl;
              }
             for(int j=leftwaves; j< Nwaves; j++){
                   SpectrumKays[ISpect[j]]=SpectrumKays[ISpect[j]]+KVEC[j]*knaughtlines[l];
                   cout << setprecision(11) << XVEC[j] << " K= " << KVEC[j] << " I= " 
                        << j <<  "  J=" << ISpect[j]  <<endl;
              }


              Xvector.erase(Xvector.begin(), Xvector.end());
              ISpect.erase(ISpect.begin(), ISpect.end());                     
           }  //end for loop over centre wavelengths
 



          if(outspec){
/*****************************************************************************************
        Output Spectra if needed
******************************************************************************************/
        // for each isotopologue output files to MolSpect directory
        // first --- the file names
        // current gas=ig  --- Yep we are still in the ig loop!

        string SpectFile;
        SpectFile=AllMols[ig];
        const string SpectDir("MolSpect");
        SpectFile.replace(0, 11, SpectDir);
        cout << SpectFile << endl;
        ostringstream oss_outcode; //output file has isotopologue number and  layer number
        int isorep1, isorep2;
        if(Iso<10){
           isorep1=Iso+48; //ASCII 0=48
           oss_outcode <<  "_Iso_"  << (char)isorep1;
          // cout << "Iso Number=" << Iso << "  character="<<  (char)isorep1 << endl;
        }
        if(Iso>=10){
           int i1=Iso-10,i2=Iso%10;
           isorep1=i1+48;
           isorep2=i2+48;
           oss_outcode <<  "_Iso_"  << (char)isorep1 << (char)isorep2;
          // cout << "Iso Number=" << Iso << endl;
         }

         string modstring;
         modstring=oss_outcode.str();
         cout << modstring << endl;
         int stringsize,ipos;
         stringsize=SpectFile.size();
         ipos=stringsize-4; //".dat" at the end
         SpectFile.insert(ipos,modstring);
         cout << SpectFile << endl;    
         } //endif for outspec




/*****************************************************************************************
        Work out k Distributions
******************************************************************************************/




           

                          
           SpectrumWaves.erase(SpectrumWaves.begin(), SpectrumWaves.end());
           CentreWaves.erase(CentreWaves.begin(), CentreWaves.end());
           CentreLines.erase(CentreLines.begin(), CentreLines.end());
           }  //end loop over layers nlay
          } // end loop over isotopologues k=0, k< nicemax




/****************************************************************************************
        Clear memory
******************************************************************************************/
        for(int i=0;i<vecsize;i++){
            cout << "gas humber=" << ig+1 << " Isotopologue=" << i+1 << "  code=" <<nicecode[ig][i]
                 << "  number of lines=" << linecentres[i].size()<< endl;
                 linecentres[i].erase(linecentres[i].begin(), linecentres[i].end()); 
                 linestrengths[i].erase(linestrengths[i].begin(), linestrengths[i].end()); 
                 linegammaAir[i].erase(linestrengths[i].begin(), linestrengths[i].end()); 
                 linegammaSelf[i].erase(linestrengths[i].begin(), linestrengths[i].end());  
                 lineLSE[i].erase(linestrengths[i].begin(), linestrengths[i].end());       
                 lineTdep[i].erase(linestrengths[i].begin(), linestrengths[i].end());     
                 lineShift[i].erase(linestrengths[i].begin(), linestrengths[i].end());            
            }

      }
      //end loop over HITRAN par files ig=0 to ngas

         } //endif for calcspec

     //delete respval and lambdaresp - add later



  return 0;
}


