//-------------------------------------------------------------------------------
// Copyright 2011 Christopher Godsalve.
// All Rights Reserved.
//
// Permission to use, copy, modify and distribute this software (if not modified) and its
// documentation for educational, research and non-profit purposes, without fee,
// and without a written agreement is hereby granted, provided that the above
// copyright notice, this paragraph and the following three paragraphs appear in all copies.
// 
//
// To request permission to incorporate this software into commercial products
// contact Dr C. Godsalve, 42 Swainstone Road, Reading, Berks, UK or by email at
// seagods@btinternet.com or seagods@hotmail.com.
//
// IN NO EVENT SHALL CHRISTOPHER GODSALVE BE LIABLE TO ANY PARTY FOR
// DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING 
// LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, 
// EVEN IF CHRISTOPHER GODSALVE HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
//
// CHRISTOPHER GODSALVE SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN `AS IS' BASIS, AND CHRISTOPHER 
// GODSALVE HAS NO OBLIGATIONS TO PROVIDE MAINTAINANCE, SUPPORT, UPDATES, 
// ENHANCEMENTS, OR MODIFICATIONS IF HE CHOOSES NOT TO DO SO.
//--------------------------------------------------------------------------------
#include <stdio.h>  //standard input output
#include <iostream> // input output (cin cout etc)
#include <fstream>  //file stream
#include <stdlib.h> //string conversion + malloc and calloc
#include <math.h>
#include <iomanip>  //setprecision
#include <sstream>  //string stream
#include <string>   // nice string operations
#include <vector>   // a bit like a souped up stack
#include <sys/stat.h>  //POSIX stuff including file queries
#include <new>   // use nothrow
using namespace std;

// fortran programs, use nm (name mangle) on object files and .so files
// to get calling name. Usually the same name with all lower case letters
// and a trailing underscore at the end. Also, fortran always passes data to
// functions and subroutines by address, and square arrays A_ij in fortran are A_ji in C.

extern "C" {

// slatec fortran B-spline routines, k=4 is cubic
// dbint4 calculates the spline representation
// SUBROUTINE DBINT4 (X, Y, NDATA, IBCL, IBCR, FBCL, FBCR, KNTOPT, T, BCOEF, N, K, W)
void dbint4_(double*, double*, int*, int*, int*,  double*, double*,
              int*, double*, double*, int*, int*, double*);

// dbvalu evalutes the spline or a derivative of the function
// dbspev similar, but calculates a vector of derivatives
// DOUBLE PRECISION FUNCTION DBVALU (T, A, N, K, IDERIV, X, INBV, WORK)
double dbvalu_(double*, double*, int* ,int*, int* ,double*, int*, double*);


// toms library Gauss quadrature
// SUBROUTINE WEIGHTCOEFF(N,Q,E,EPS,W,X,WORK)
void weightcoeff_(int*, double*, double*, double*, double*,
                       double*, double*);

//HITRAN Partition Function
void bd_tips_2003_(int&, double&, int&, double&, double&);

//calculates Voigt profile.
void humlik_(int&, double*, double&, double*);
}

// Either prompts user for input values, or reads an input file.
void Questioner( bool&,
                int&, int&, int&, int&, bool&, double&,
                double&, double&, double&, double&, bool&,  int&, bool&, bool&, int&,
                int&, int&, int&, bool&, bool&, bool&, 
                int&, double&, double&,
                bool&, double&, double&, double&, double&, double&, int&, 
                bool&, bool&, bool&, bool&, bool&, bool&,
                int&, int&, int&, int&,
                const char*);

// Either calculates humidity from number density or number density from humidity.
// The values can be over liquid water or ice.
void WaterVap(double* , double* , double* ,
              double* , double &, double &, double &, int &, int &, int &);

int main(int argc, char* argv[]){

  bool verbose;   //verbose mode for input, read from file if false;
  //verbose output bools

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
   Stage 4:  4A: Loop over Hitran par files, read in line centres.
             4B: Work out where to calculate the spectrum.
             4C: calculate the spectrum.
***************************************************************************************

/************************   BEGIN STAGE 1  *********************************************/
  // Four Aerosol layers, each to be split into nsplit slabs
  // Four layers are boundary layer, troposphere, stratosphere, and upper atmosphere
  int nsplit1,nsplit2,nsplit3,nsplit4;    
  // boundary layer visibility , troposphere visiblility, wavelength1, wavelength2
  double visb, vist, lambda1, lambda2, lambdacut;
  bool rfun;  //true if we are using myresp.dat or myrespX.dat as a response function file
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


  bool PlotIt; //if we are using PlotIt need to format the output for this plotting routine
  bool logplot; //if we are outputing the log of the cross section
  bool calcRef; //calculate refractive index spectrum
  bool outRef;  // output refractive spectrum files
  

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
              vist, lambda1, lambda2, lambdacut, rfun,  iatm, switchR, switchA, itypeu,
              itypes, itypet, itypeb, ocean, groundP, groundT,
              ihumid, groundtemp, groundpress,
              default_pause, HG, HB, HT, HS, HU, ngauss_height, 
              calcspec, outspec, PlotIt, logplot, calcRef, outRef,
              ngauss_correlkay, icutL,icutD,istep,
              ReadInput);  

    int nresp;
    double* lambda_resp; double* respval;
    if(rfun){
       ifstream fp_in;
       fp_in.open(RespFile, ios::in);
       if(!fp_in.is_open()){ cerr << "Failed to open file (" << RespFile << ")"  << endl; exit(1); }
       fp_in >> nresp;
        lambda_resp=new (nothrow) double[nresp]; respval=new (nothrow) double[nresp];
        if(lambda_resp==0){cout << "memory for lambda_resp failed\n" << endl; exit(10);}
        if(respval==0){cout << "memory for respval failed\n" << endl; exit(10);}
        for(int i=0; i< nresp ; i++){
         fp_in >> lambda_resp[i] >> respval[i];
        }
        lambda1=lambda_resp[0]; lambda2=lambda_resp[nresp-1];
        fp_in.close();
        cout << "Response function file overrides wavelengths in Dors1in.dat\n";
        cout << "Wavelengths used are " << lambda1 << " and " << lambda2 << endl;
    }

/*   NAUGHTY STUFF   --- Hacks to do one off stuff and/or over-ride Questioner */

     bool verbose_out=true;  // output loads of info to standard output
     bool verbose_height=true;   //output height integration info to std::out
     bool verbose_spect=false;   // output spectrum details to std out

     bool  ShiftLine=0;  //switch pressure dependent line shifts on and off
                         // Make 1 or 0 depending on whether we bother with line shifts
     bool LayerChange=true;  // We will pratt about with layers
     int NLStart=9,NLStop=10;  // eg start=9, stop=10 1 layer

     int nicemax=1;  // max number of isotopes considered --- should go to Questioner!
     bool LEGEND=false;

     int ichap=5;  // 5 temperature files for Chappuis

     int iozHH=6;  // 5 temperature sets for Ozone HH

     int iHHdata=581; //581 data points in ozone HH files

     double startXHH=29164.0; double stopXHH=40798.0;
     //roughly 343nm-245nm

     double HHTemp[6]; //Temperatures for Ozone HH Bands
     HHTemp[0]=200.0; HHTemp[1]=220.0; HHTemp[2]=240.0; HHTemp[3]=260.0;
     HHTemp[4]=280.0; HHTemp[5]=300.0;


     double ChapTemp[5]; //Temperatures for Chappuis Bands
     ChapTemp[0]=218.0; ChapTemp[1]=228.0; ChapTemp[2]=243.0; ChapTemp[3]=273.0;
     ChapTemp[4]=295.0;

     double startXChap[5];  //Lowest wave number in each Chappuis file
     startXChap[0]=15384.379; startXChap[1]=19230.399; startXChap[2]=19267.451; 
     startXChap[3]=19230.399; startXChap[4]=12048.193; 

     vector<double> OzWavesHH; vector<double> OzXHH[6];  //Huggins and Hartley Ozone bands, 6 spectra
     vector<double> OzWavesChap[5]; vector<double> OzXChap[5];  //Huggins and Hartley Ozone bands, 6 spectra
     bool ozhh=false; bool ozchap=false;

     int HHlines=851;
     int ChapLines[5]; //number of lines in file - 2 (first 2 not data)
     ChapLines[0]=45552; ChapLines[1]=32552; ChapLines[2]=32452;
     ChapLines[3]=22052; ChapLines[4]=63501;

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

   double pi,abszero,R,AV,MWA,delta,nucleon,rootpi,twopi;

   pi=acos(-1.0);               //Well pi obviusly
   twopi=2.0*pi;
   rootpi=sqrt(pi);
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
   double Boltz=1.3806503e-16;  //Boltzmann's constant (cgs units)
   //speed of light in cm per sec.
   double Speedlight=2.99792458e10;



   ifstream fp_in;   //input file stream
   int gid[50];      //HITRAN gas numbers (last few are dummies)
   bool gmask[50];   //if false, ignore that gas
   fp_in.open("GasMask.dat", ios::in);


   if(!fp_in.is_open()){ cout << "Failed to open file (GasMask.dat)" << endl;  exit(1);}

   fp_in >> gid[0] >> gid[1] >> gid[2] >> gid[3] >> gid[4] 
         >> gid[5] >> gid[6] >> gid[7] >> gid[8] >> gid[9];
   fp_in >> gmask[0] >> gmask[1] >> gmask[2] >> gmask[3] >> gmask[4] 
         >> gmask[5] >> gmask[6] >> gmask[7] >> gmask[8] >> gmask[9];
   fp_in >> gid[10] >> gid[11] >> gid[12] >> gid[13] >> gid[14] 
         >> gid[15] >> gid[16] >> gid[17] >> gid[18] >> gid[19];
   fp_in >> gmask[10] >> gmask[11] >> gmask[12] >> gmask[13] >> gmask[14] 
         >> gmask[15] >> gmask[16] >> gmask[17] >> gmask[18] >> gmask[19];
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

   double PartPress=0;

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
    double GasPPM[idata][ngas];  //Gas concentration of the selected gases in ppm
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
           if(i==up_to_top-1){
              SlabThick[i]=HBDRY[4]-HeightFloors[i];
                if(verbose_height)cout << "Floors " << HeightFloors[i] << " "  <<floorpress[i] << " " << FloorTemp[i] <<
               " " << SlabThick[i] << endl;
               }
              else{
              SlabThick[i]=HeightFloors[i+1]-HeightFloors[i];
                if(verbose_height){
                  cout << "Floors " << HeightFloors[i] << " "  <<floorpress[i] << " " << FloorTemp[i] <<
                  " " << SlabThick[i] << endl;
                //  cout << i << " " << HeightFloors[i] << endl;
                  }
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
        double Gas[iheight][ngas];  //fine grid array of gas ppm for all gases which are not ignored via gasmask

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
               // divide by a 1000 and DensGas is in Kg m^{-2} (remember nucleon is in grammes.
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
       //mass of air (per square metre) in each layer and the total mass of all layers
       double massgas[up_to_top]; double totalmass=0.0;
       //mass of each constituent in each layer and the total mass
       double massgases[up_to_top][ngas]; double  totalmasses[ngas];
       double densgases[up_to_top][ngas]; 

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
        if(verbose_height)cout << "Averages " << AvTemp[nsplit1-i] << "  " << AvPress[nsplit1-i] << "  " << totalmass <<endl;
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
        if(verbose_height)cout << "Averages " << AvTemp[nsplit2-i+idown] << "  " << AvPress[nsplit2-i+idown] << "  " << totalmass << endl;
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
        if(verbose_height)cout << "Averages " << AvTemp[nsplit3-i+idown] << "  " << AvPress[nsplit3-i+idown] << "  " << totalmass << endl;
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
        if(verbose_height)cout <<"Averages " <<  AvTemp[nsplit4-i+idown] << "  " << AvPress[nsplit4-i+idown] << "  " << totalmass <<   "  layer=" << nsplit4-i+idown << "  idown=" << idown <<endl;
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
        //

        for(int k=0; k<ngas; k++){
          for(int i=0; i<up_to_top; i++){
             massgases[i][k]=massgases[i][k]/10.0/correction;
          }
          totalmasses[k]=totalmasses[k]/10.0/correction;
        }
        totalmass=totalmass/10.0/correction;



        if(verbose_out)cout << "gases in each layer --- top down --- g per cm^2\n";
        // Will need gas density in g cm^{-3} for refractive index

        for(int i=0; i<up_to_top; i++){ //top down actually
            for(int k=0; k < ngas ; k++){
                  densgases[i][k]=massgases[i][k]/SlabThick[up_to_top-1-i]/100.0;
                  //divided by 100 because slabthick in metres - want cm^{-3}
             }
        }
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

/******************************* BEGIN STAGE 4A     ************************************/
/********************  FOR EACH GAS -- FOR EACH ISO -- FOR EACH LAYER  *****************/


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
      
      // this nu_cut just widens the wavelength range to include lines from outside.

      nu_cut=lambdacut/lambda1*nu1;
      nu1X=nu1-nu_cut;
      nu_cut=lambdacut/lambda2*nu2;
      nu2X=nu2+nu_cut;
      lambda1X=1.0e4/nu2X; lambda2X=1.0e4/nu1X;

      //later on nu_cut is used as line dependent cut off

      cout << "nu1 nu2 " << nu1 << "  " << nu2  <<endl;
      cout << "lambda1  lambda2  " << lambda1 << "  " << lambda2 << endl;
     
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
      
      string  PARFILE="HITRAN/   hit08.par";
      double lambda;   // a wavelength
      
      string line;
      string entry[19];
      istringstream iss[19]; //19 input stringstreams
      
      // FOR EACH GAS
      // Note igx is NOT the HITRAN gas number
      // gases[igx] IS the HITRAN molecule number
      // more confusion, Hitran numbers start at 1, C numbers start at zero
      // hence int molecule and int molek

      for(int igx=0; igx<ngas; igx++){ //loop to read from HITRAN files.
        
        int molecule=gases[igx];
        int molek=molecule-1;

        int ireplace1=molecule/10; int  ireplace2=molecule%10;
        ireplace1=ireplace1+48; ireplace2=ireplace2+48; //ascci '0' is 48
        ostringstream replacewith;
        replacewith << (char)ireplace1 << (char)ireplace2 << '_';
        PARFILE.replace(7,3,replacewith.str());
        cout <<"PARFILE=" << PARFILE << endl;
        fp_in.open(PARFILE.c_str(),ios::in);
        if(!fp_in.is_open()){
                             cout << igx << "  ngas=" << ngas <<endl;
                             cout <<" can't open HITRAN file \n";
                             cout << PARFILE.c_str() << endl;  exit(1);}
        int vecsize=nice[molek];
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
        //vectors for each isotopologue, vecsize <=nicemax depending on nice[molek]
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
       // We assume we are using cross section files and NOT line by line
       // HITRAN .par files for gases 30,35, and 42
       // Ozone has both a .par file and a cross section file
        bool UseParFile=true;
        if( (molecule==30) || (molecule==35) || (molecule==42) )UseParFile=false;

        if(UseParFile){

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

          if(nu > nu1X && nu < nu2X){  // if nu is in range
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
          linestrengths[Iso-1].push_back(S_intense/rootpi);  //rootpi from HUMLIK
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
          if(verbose_spect)cout << "Einstein_A="  << Einstein_A<< endl; 
          //not used for local thermodynamic equilibrium

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
          }  //endif for lambda within range (nu>nu1X nu < nu2X)
          iss[2].clear();
          } //endif for isotope less than nicemax
          iss[1].clear();

          } //end  if not end of file while loop
        fp_in.close();
        //Finished reading in HITRAN data
        //At last we have the line centres for this particular gas.

        cout << "Finished with PARFILE=" << PARFILE << endl;

        } //endif for UseParFile





        cout <<" gid=" << molecule << "  gmask= " << gmask[molek] << endl;
             cout <<"nu1="  << nu1 << " nu2=" << nu2 <<endl;

//**************************** Now for X section files  *****************************************************/


        if(molecule==3){

          double hello=-1.0; double goodbye=-1.0;

//         startXHH=29164.0 (343nm) and stopXHH=40798.0 (245nm) 


          if(nu1 <= startXHH  && nu2 > startXHH){
               hello=startXHH;
               ozhh=true;
               if(nu2 < stopXHH)goodbye=nu2; 
                         else goodbye=stopXHH; }
          if(nu1 > startXHH  && nu1 < stopXHH){
               hello=nu1;
               ozhh=true;
               if(nu2 < stopXHH)goodbye=nu2; 
                         else goodbye=stopXHH; }

          if(ozhh){



          if(goodbye<0){cout << "Ozone UV gone wrong\n"; exit(0);}

          fp_in.open("HITRAN/UV_X/O3-UV04.xsc", ios::in);

          if(!fp_in.is_open()){ cerr << "Failed to open file O3-UV04,xsc"  << endl; exit(1); }



          //OK - We know the data format
          string ozdata;  // We have 6 temperatures for the ozone Hartley and Huggins spectra
          string entry;
          istringstream input_ozo;

          double OzoX;


          for(int ioz=0; ioz<iozHH; ioz++){
             getline(fp_in,ozdata);   // 6 headers
             cout << ozdata << endl;
             double OzW1=startXHH;

             for(int id=0; id<iHHdata; id++){
                getline(fp_in,ozdata);

                for(int idat=0; idat<10; idat++){  //581 lines have 10 data entries
                  entry=ozdata.substr(idat*10+1,9);
                  input_ozo.str(entry);
                  input_ozo >> OzoX;  input_ozo.clear();
                  if(hello <= OzW1 && OzW1<=goodbye){
                     if(ioz==0)OzWavesHH.push_back(OzW1);
                     OzXHH[ioz].push_back(OzoX);                     
                  }
                  OzW1+=2.0;
                } //end idat loop
             }  //  have read in 581 lines of 10 entries


               getline(fp_in,ozdata);  // last line is only 8 numbers
               for(int idat=0; idat<8; idat++){  //581 lines have 10 data entries
                  entry=ozdata.substr(idat*10+1,9);
                  input_ozo.str(entry);
                  input_ozo >> OzoX;  input_ozo.clear();
                  if(hello <= OzW1 && OzW1<=goodbye){
                     if(ioz==0)OzWavesHH.push_back(OzW1);
                     OzXHH[ioz].push_back(OzoX);                     
                  }
                  OzW1+=2.0;
                } //end idat loop for last line of 8
                
          }  // end 6 ioz temperatures loop


          fp_in.close();

          } // endif ozhh true





         // Chappuis Bands!
         double stopXChap=startXHH; //Use the HITRAN data for high wavenumbers where possible


/*   startXChap already declared and initiated as

     startXChap[0]=15384.379; startXChap[1]=19230.399; startXChap[2]=19267.451; 
     startXChap[3]=19230.399; startXChap[4]=12048.193; 

     corresponding (roughly) to 650nm, 520nm, 519nm, 520nm,  and 830nm
*/ 

         double hellochap[ichap]; double goodbyechap=-1.0;

         hellochap[0]=-1.0; hellochap[1]=-1.0; hellochap[2]=-1.0;
         hellochap[3]=-1.0; hellochap[4]=-1.0;
         
         for(int ic=0; ic< ichap; ic++){
          if(nu1 <= startXChap[ic]  && nu2 > startXChap[ic]){
               hellochap[ic]=startXChap[ic];
               ozchap=true;
               if(nu2 < stopXChap)goodbyechap=nu2; 
                         else goodbyechap=stopXChap; }
          if(nu1 > startXChap[ic]  && nu1 < stopXChap){
               hellochap[ic]=nu1;
               ozchap=true;
               if(nu2 < stopXChap)goodbyechap=nu2; 
                         else goodbyechap=stopXChap; }
         }

         if(ozchap){


         string  ChapFile[5];
         ChapFile[0]="HITRAN/UV_X/SMPO/SMPO1.cs"; ChapFile[1]="HITRAN/UV_X/SMPO/SMPO2.cs";
         ChapFile[2]="HITRAN/UV_X/SMPO/SMPO3.cs"; ChapFile[3]="HITRAN/UV_X/SMPO/SMPO4.cs";       
         ChapFile[4]="HITRAN/UV_X/SMPO/SMPO5.cs"; 

         string chapline;
         double wavechap, xchap;
         

         for(int ifile=0; ifile<ichap;ifile++){
            fp_in.open(ChapFile[ifile].c_str(), ios::in);

            if(!fp_in.is_open()){ cerr << "Failed to open file Chappuis file " 
                                       << ifile << endl; exit(1); }
            getline(fp_in,chapline);
            cout << chapline << endl;
            getline(fp_in,chapline);
            cout << chapline << endl;

            for(int id=0; id< ChapLines[ifile]; id++){
               fp_in >> wavechap >> xchap;  
               if( (hellochap[ifile] <= wavechap) && wavechap <=goodbyechap){
                      OzWavesChap[ifile].push_back(wavechap);
                      OzXChap[ifile].push_back(xchap);
               } //endif for in range
            } // end read in data file  


            fp_in.close();
          }  //end ifile < ichap loop
          }  //endif ozchap



          if(ozhh){
             cout <<"Molecule=3 OZONE HH\n";
             cout <<"startXHH=" << startXHH << "  stopXHH=" << stopXHH <<endl;
             cout <<"hello="  << hello << " goodbye=" << goodbye <<endl;
             cout <<"nu1="  << nu1 << " nu2=" << nu2 <<endl;
       
            cout << "Size of Ozwaves=" << OzWavesHH.size() << endl;
            for(int ioz=0; ioz<iozHH; ioz++){
              cout << "Size of OzXHH=" << OzXHH[ioz].size() << endl; }
          }

          if(ozchap){
             cout <<"Molecule=3 OZONE Chap\n";
             for(int ic=0; ic< ichap; ic++){
             cout <<"startXChap=" << startXChap[ic] << "  stopXChap=" << stopXChap <<endl;
             cout <<"hellochap="  << hellochap[ic] << " goodbyechap=" << goodbyechap <<endl;
             }
             cout <<"nu1="  << nu1 << " nu2=" << nu2 <<endl;
      
            for(int ic=0; ic<5; ic++){
              cout << "Size of OzWChap OzXChap=" << OzWavesChap[ic].size() 
                       << "  " << OzXChap[ic].size() << endl; }
      
          }

          cout <<"Exit at at THIS LINE ARSE!\n";  exit(0);

          } // endif molecule=3


        if(molecule==30){
          //Read in X-section for IR Sulphur Hexaflouride (SF6) if wavelength range requires it
          //wn range 925-955
        }
        if(molecule==35){
          //Read in X-section for IR Chlorine Nitrate (ClONO2) if wavelength range requires it
          //wn ranges 750-830, 1260-1320, 1680-1790

        }

        if(molecule==42){
          //Read in X-section for IR Carbon Tetraflouride (CF4) if wavelength range requires it
          //wn range 1250-1290

        }
/*******************************END STAGE 4A *****************************************/

/*******************************BEGIN STAGE 4B ***************************************/
/************************ 4B BEGINS INSIDE LOOP INSIDE  LOOP ****************************/
/********************     FOR EACH ISO -- FOR EACH LAYER            *****************/

       if(UseParFile){

        //FOR EACH ISO
       for(int k=0; k<vecsize; k++){   //loop over isotopologues

        //gas and isotopologue for BD_TIPS_2003.f
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

        bd_tips_2003_(molecule, RefTemp, k_fort, Refgee_i, QT);
        ofstream OutSpect; ofstream RefSpect; ofstream TauSpect; ofstream SigSpect;

        if(!LayerChange){
           NLStart=0;
           NLStop=up_to_top;
        }
        if(NLStart <0 || NLStop>up_to_top){
             cout << "Layers Wrong" << endl; exit(1);
         }


    
        //FOR EACH LAYER
        for(int nlay=NLStart; nlay < NLStop; nlay++){

            TempNlay=AvTemp[nlay];
            PressNlay=AvPress[nlay]*100.;

            //AvPress in mbar, Ref Pressure in Pascal



          //BD_TIPS_2003 Routine provided by HITRAN.

            bd_tips_2003_(molecule, TempNlay, k_fort, gee_i, QT);
           // cout << "gi=" << gee_i << "  QT=" << QT <<  endl;

           //get seg faults if declared on stack

            double* wavenumlines; double* strengthlines; double* alphaLorentzA;
            double* alphaDop; double* widthlines;

        //     cout << "about to declare wave\n"; 
            wavenumlines=new (nothrow) double[nlines];    // temporary copies
        //     cout << "about to declare strength\n";
            strengthlines=new (nothrow) double[nlines];
       //      cout << "about to declare Lor\n";
            alphaLorentzA=new(nothrow) double[nlines];  //for air -add partial self and partial pressures later;
       //      cout << "about to declare Dop\n";
            alphaDop=new (nothrow)  double[nlines];
      //       cout << "about to declare widthlines\n";
            widthlines=new (nothrow) double[nlines];
      //       cout << " Done! \n";
            if(wavenumlines==0){cout << "memory failed for wavenumlines\n"; exit(10);}
            if(strengthlines==0){cout << "memory failed for strengthlines\n"; exit(10);}
            if(alphaLorentzA==0){cout << "memory failed for alphaLorentz\n"; exit(10);}
            if(alphaDop==0){cout << "memory failed for alphaDop\n"; exit(10);}
            if(widthlines==0){cout << "memory failed for widthlines\n"; exit(10);}



 
            //First we find all the wave numbers at which to calculate
            //the entire spectrum over nu1X to nu2X and stuff them in SpectrumWaves;
            int kountwaves=0;  //wavelengths which are line centres

            for(int l=0; l<nlines;l++){
              if(QT>0){


  //            cout <<  "  k and l " << k << "  " << l << endl; 

              wavenumlines[l]=linecentres[k][l];


           //   cout << "new centre at " << linecentres[k][l] <<  "  k and l " << k << "  " << l << endl; 
              //This is a half width -- dont want it                
//              alphaDop[l]=(wavenumlines[l]-lineShift[k][l]*PressNlay/RefPress)
//                          /Speedlight*sqrt(2.0*Boltz*TempNlay*log2/(AllMolW[molek][k]*nucleon));
               //Want 1/e width =half width divided by sqrt(log(2))
              //  awkward - but for now multiply <P> by 100 since RefPress in Pascals

               if(ShiftLine){
               alphaDop[l]=(wavenumlines[l]-lineShift[k][l]*PressNlay/RefPress)
                          /Speedlight*sqrt(2.0*Boltz*TempNlay/(AllMolW[molek][k]*nucleon));}
               else{
               alphaDop[l]=wavenumlines[l]
                          /Speedlight*sqrt(2.0*Boltz*TempNlay/(AllMolW[molek][k]*nucleon));
               }
 
              // HITRAN 96 Appendix A. "Hitran Parameters: Definitions and Usage"
              // Ref Paper, HITRAN 96 eqn.A.11 and A.12
              // Partial Pressure
              //
              //    ARSE!
              //
              //  We assume there is only air broadening for k>0

              if(k==0){
                 PartPress=massgases[nlay][igx]/massgas[nlay]*MWA/AllMolW[molek][k];

                 alphaLorentzA[l]=PressNlay/RefPress*(
                 linegammaAir[k][l]*(1.0-PartPress)+linegammaSelf[k][l]*PartPress
                   )*pow(RefTemp/TempNlay,lineTdep[k][l]);
               }
               else{
                 //don't bother with self broadening for other isotopes [k]
                 alphaLorentzA[l]=linegammaAir[k][l]*PressNlay/RefPress
                                *pow(RefTemp/TempNlay,lineTdep[k][l]);
               }
/*            // diagnostic write
              cout << alphaLorentzA[l] << "  P=" << PressNlay  
                    << "  " << alphaDop[l]  << " Part Pres=  " << PartPress << 
                   "  " <<  pow(RefTemp/TempNlay,lineTdep[k][l]) <<
                   "  " <<  linegammaAir[k][l] << "  " << linegammaSelf[k][l] << 
                   "  " << linestrengths[k][l] << endl;
*/

              // Ref Paper HITRAN 96 eqn A.10
              strengthlines[l]=linestrengths[k][l]*QTref[molek][k]/QT
                      *exp(-C2*lineLSE[k][l]/TempNlay)/exp(-C2*lineLSE[k][l]/RefTemp)
                      *( 1.0-exp(-C2*linecentres[k][l]/TempNlay) )*( 1.0-exp(-C2*linecentres[k][l]/RefTemp) );


            //    cout << strengthlines[l] << "=strength\n";  //compare with JavaHAWKS sometime!
              }
              else{
                 cout << "Error: BD_TIPS_2003 returned QT=" << QT << endl;
                 exit(1);
              }
              //now we have the x and y vectors for all the lines
            }  //end loop over nlines
 
/*******************************BEGIN STAGE 4B PROPER - WORK OUT SPECTRUMWAVES **********************/
            vector<double> SpectrumWaves; 
            vector<double>SpectrumKays;  //extinction
            vector<double>SpectrumDeltaN;  //Refractive index minus one

            vector<int> CentreWaves;  //which values of nu in SpectrumWaves are at line centres?

            SpectrumWaves.clear();
            CentreWaves.clear();

            SpectrumKays.clear();
            SpectrumDeltaN.clear();
 
            double startnu,stepnu,tempnu,currentnu0,
            currentwidth,currentwidthL,currentwidthD,nextnu0,last_stopnu, stopnu, Wtol;


            Wtol=0.000001;  //because of fortran F12.6  for wave number in par file

            int icurrent=0;
            last_stopnu=nu1X;
            bool skipline=false;
            double midpoint;

            while(icurrent < nlines){  
              //we get some "double lines" -skip the next line if it it's the same nu0
              if(skipline){
                 skipline=false;
                 //icurrent has already been incremented

                 //we are just calculating where on the nu axis our fine detail spectrum is calculated
                 //the two transitions with the same frequency have different widths and strengths
                 //so it must appear in the CentreWaves list. It will be the last centre wave number.
                 CentreWaves.push_back(CentreWaves[CentreWaves.size()-1]);
                 //current nu0 is the same as last nu0's
                 currentnu0=wavenumlines[icurrent];
                 if(icurrent<nlines){
                    nextnu0=wavenumlines[icurrent+1];
                 }

                 if(ShiftLine){
                 currentnu0=currentnu0-lineShift[k][icurrent]*PressNlay/RefPress;}
                 else{
                 currentnu0=currentnu0;}


                 //set current width
                 currentwidthD=alphaDop[icurrent]*icutD;
                 currentwidthL=alphaLorentzA[icurrent]*icutL;
                 int whichcut;
                //widthlines[icurrent] has not been initialised.
                //It could contain any garbage imaginable.
                 if(currentwidthL<currentwidthD){
                       currentwidth=currentwidthD;
                       whichcut=icutD;}
                 else{
                     currentwidth=currentwidthL;
                     whichcut=icutL;
                     }

                     widthlines[icurrent]=currentwidth;
                     icurrent++;
              }   // end if skipline is true
              else{ 
                  //skipline is false
                  currentnu0=wavenumlines[icurrent];
                  if(icurrent<nlines-1){               
                    nextnu0=wavenumlines[icurrent+1];
                    if(currentnu0==nextnu0){
                      skipline=true;
                    }
                   }

                   //skipline is mostly still false, but occasionally nextnu0 is same
                   //as currentnu0  --- but occasionally skipline will be reset to true
              if(ShiftLine){
                 currentnu0=currentnu0-lineShift[k][icurrent]*PressNlay/RefPress;
              }

              //set current width
              currentwidthD=alphaDop[icurrent]*icutD;
              currentwidthL=alphaLorentzA[icurrent]*icutL;
              int whichcut;
              if(currentwidthL<currentwidthD){
                    currentwidth=currentwidthD;
                    whichcut=icutD;}
              else{
                    currentwidth=currentwidthL;
                    whichcut=icutL;}
          
              stepnu=currentwidth/istep;
              stepnu=stepnu/whichcut;
              nu_cut=currentwidth;

              widthlines[icurrent]=currentwidth;

              if(currentnu0-nu_cut<last_stopnu){
                   startnu=last_stopnu+stepnu; } else {startnu=currentnu0-nu_cut;}

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
                   kountwaves++;}
              else{
                   if(SpectrumWaves.size()>0){
                   SpectrumWaves.pop_back();
                   kountwaves--;
                   SpectrumWaves.push_back(currentnu0);
                   CentreWaves.push_back(kountwaves);
                   kountwaves++;
                   } //endif 
                   }


              bool double_last_line=false;

              //Start If Else
              if(!skipline){
                 midpoint=(currentnu0+nextnu0)/2.0;}
              else{

                 if(icurrent+2<nlines){
                   nextnu0=wavenumlines[icurrent+2];}
                 else{
                   double_last_line=true;}

                midpoint=(currentnu0+nextnu0)/2.0;
              } //End If Else

              if(icurrent==nlines-1){               
                  stopnu=currentnu0+nu_cut;
              }//last line -- nothing to the right

              //Start If Else
              if(double_last_line){
                  stopnu=currentnu0+nu_cut;
              }
              else{
                if(midpoint > currentnu0+nu_cut){
                  stopnu=currentnu0+nu_cut;}
                  else{
                  stopnu=midpoint;
                }
              } //End If Else

              tempnu=currentnu0; //start at line centre and step up
              bool past_it=false;

              while(tempnu<stopnu-Wtol && tempnu < nextnu0-Wtol){

               tempnu+=stepnu;
               if(tempnu>=nextnu0)past_it=true;
               if(!past_it){
                  SpectrumWaves.push_back(tempnu); 
                  kountwaves++;}
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
           } //end of while icurrent << nlines



/*******************************END STAGE 4B  -- WE HAVE SPECTRUMWAVES ***************/
/*******************************BEGIN STAGE 4C  KAYWAVES *****************************/    


           int templength=SpectrumWaves.size();

           double myzero=0.0;
           for(int l=0; l<templength-1; l++){
              SpectrumKays.push_back(myzero);
              SpectrumDeltaN.push_back(myzero);
           }    

           double yfac;
           double   xtemp, width, alphaD, alphaL;
           int leftwaves, rightwaves, iwave, SpectrumSize;
           int nulength=CentreWaves.size();
           int speclength=SpectrumWaves.size();

           
           for(int l=0; l<nulength-1; l++){
          
              if(ShiftLine){
              currentnu0=wavenumlines[l]-lineShift[k][l]*PressNlay/RefPress;}
              else{
              currentnu0=wavenumlines[l];}
            
              width=widthlines[l];   //cutoff from previous loop
              alphaD=alphaDop[l];
              alphaL=alphaLorentzA[l];
              yfac=alphaL/alphaD; 
              //Need to calculate Xvector for HUMLIK
              //But it needs to be ordered. We will have a normal array (ordered)
              // and a C++ vector (disordered)
              vector<double> Xvector;
              vector<int>  ISpect;
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
                   if(xtemp>currentnu0-width){
                   Xvector.push_back(xtemp); 
                   ISpect.push_back(iwave-leftwaves);
                   leftwaves++;}
              } 

              rightwaves=1;
              //linecentre then everything within range on right

              if(iwave+rightwaves<speclength-1){
              xtemp=SpectrumWaves[iwave+rightwaves];
              Xvector.push_back(xtemp); 
              ISpect.push_back(iwave+rightwaves);
              rightwaves++;

              while(xtemp<currentnu0+width && iwave+rightwaves<speclength){
                   xtemp=SpectrumWaves[iwave+rightwaves];
                  // cout << setprecision(11) << xtemp << endl;
                   if(xtemp<currentnu0+width){
                   Xvector.push_back(xtemp); 
                   ISpect.push_back(iwave+rightwaves);
                   rightwaves++;}
              } 
              } //endif for going past spectrum limit
              // the rightwaves block looks like a 1 to N vector, not a 0 to N-1 so shave 1 off rightwaves
              rightwaves--;
              int Nwaves=leftwaves+rightwaves;

              //remember the first entry in Xvector is the line centre

              double* XVEC=new double[Nwaves];
              double* KVEC=new double[Nwaves];
    //          double XVEC[Nwaves],KVEC[Nwaves];
              for(int j=leftwaves-1; j>=0; j--){
              //    cout << leftwaves-1-j << "  " << j << endl;
                  XVEC[leftwaves-1-j]=(Xvector[j]-Xvector[0])/alphaD;
                  KVEC[leftwaves-1-j]=0.0;
              }
             for(int j=leftwaves; j< Nwaves; j++){
                  XVEC[j]=(Xvector[j]-Xvector[0])/alphaD;
                  KVEC[j]=0.0;
              }
              //calculate voigt profile
              humlik_(Nwaves, XVEC, yfac, KVEC);



              for(int j=0; j<leftwaves; j++){
                  SpectrumKays[ISpect[leftwaves-1-j]]=SpectrumKays[ISpect[leftwaves-1-j]]
                                                     +(KVEC[j]/alphaD)*strengthlines[l];

                  if(outRef)SpectrumDeltaN[ISpect[leftwaves-1-j]]=SpectrumDeltaN[ISpect[leftwaves-1-j]]
                           +(currentnu0-SpectrumWaves[ISpect[leftwaves-1-j]])
                           /SpectrumWaves[ISpect[leftwaves-1-j]]/alphaL
                           *(KVEC[j]/alphaD)*strengthlines[l]/twopi;

            /*      cout << setprecision(11) << XVEC[j] << " K= " << KVEC[j] << " Left-1-j  "
                       << leftwaves-1-j << "  ISpect=" << ISpect[leftwaves-1-j] << "  "
                       << Xvector[j] <<  "  " << Xvector[0]  << "  " << l << endl; */
              }

             for(int j=leftwaves; j< Nwaves; j++){
                   SpectrumKays[ISpect[j]]=SpectrumKays[ISpect[j]]+(KVEC[j]/alphaD)*strengthlines[l];

                   if(outRef)SpectrumDeltaN[ISpect[j]]=SpectrumDeltaN[ISpect[j]]
                           +(currentnu0-SpectrumWaves[ISpect[j]])
                           /SpectrumWaves[ISpect[j]]/alphaL
                           *(KVEC[j]/alphaD)*strengthlines[l]/twopi;
              } 

              //cout << currentnu0 << "  " << wavenumlines[l] << "  " << leftwaves << "  " << rightwaves << endl;
              delete[] XVEC; delete[]KVEC;
              Xvector.erase(Xvector.begin(), Xvector.end());
              ISpect.erase(ISpect.begin(), ISpect.end());                 
           }  //end for loop l<nulength over centre wavelengths


            delete[] wavenumlines;
            delete[] strengthlines;
            delete[] alphaLorentzA;
            delete[] alphaDop;
            delete[] widthlines;
 
          if(outspec){

/*****************************************************************************************
        Output absorption cross section Spectra if needed
******************************************************************************************/
        // for each isotopologue output files to MolSpect directory
        // first --- the file names
        // current gas=gases[igx]  --- Yep we are still in the igx loop!
        bool firstwrite=false;
        if(nlay==NLStart)firstwrite=true;
        string SpectFile;
        string SpectFileT;
        string SpectFileN;
        string SpectFileS;
        SpectFile=AllMols[molek];
        SpectFileT=AllMols[molek];
        SpectFileN=AllMols[molek];
        SpectFileS=AllMols[molek];
      //  cout << AllMols[0] << "  " << AllMols[1] <<  "  " << AllMols[2]  << "  molecule=  " << molecule 
      //  << " gases ig= " << gases[igx] << endl;
     //   cout << SpectFile << endl; exit(0);
        const string SpectDir("MolSpect");
        const string SpectDirT("TauSpect");
        const string SpectDirN("RefSpect");
        const string SpectDirS("SigSpect");
        SpectFile.replace(0, 11, SpectDir);
        SpectFileT.replace(0, 11, SpectDirT);
        SpectFileN.replace(0, 11, SpectDirN);
        SpectFileS.replace(0, 11, SpectDirS);  //Sigma spect will contain Rayleigh Depth from H2O scattering.
        cout << SpectFile << endl;
        cout << SpectFileT << endl;
        cout << SpectFileN << endl;
        cout << SpectFileS << endl;
        ostringstream oss_outcode; //output file has isotopologue number and  layer number
        int isorep1, isorep2;

        if(k+1<10){
           isorep1=k+1+48; //ASCII 0=48
           oss_outcode <<  "_Iso_"  << (char)isorep1;
          // cout << "Iso Number=" << k+1 << "  character="<<  (char)isorep1 << endl;
        }

        if(k+1>=10){
           int i1=k+1-10,i2=(k+1)%10;
           isorep1=i1+48;
           isorep2=i2+48;
           oss_outcode <<  "_Iso_"  << (char)isorep1 << (char)isorep2;
          // cout << "Iso Number=" << k+1 << endl;
         }

         string modstring;
         modstring=oss_outcode.str();
         cout << modstring << endl;

         int stringsize,ipos;
         stringsize=SpectFile.size();
         ipos=stringsize-4; //".dat" at the end
         SpectFile.insert(ipos,modstring);
         SpectFileT.insert(ipos,modstring);
         SpectFileN.insert(ipos,modstring);
         SpectFileS.insert(ipos,modstring);
         cout << SpectFile << "  " << nlay << "  "  << k << endl;  
         cout << SpectFileT << "  " << nlay << "  "  << k << endl;  
         cout << SpectFileN << "  " << nlay << "  "  << k << endl;  
         cout << SpectFileS << "  " << nlay << "  "  << k << endl;  

         if(firstwrite){
           int Status=-1000;
           struct stat FileInfo;
           Status=stat(SpectFile.c_str(),&FileInfo);
           //returns 0 if successfully get File Infor, -1 if it fails to find the file
           if(Status==0){
             cout << "File " << SpectFile << " already exists-We will not overwrite it or append to it!\n";
             exit(0);}
          } //endif firstwrite

         firstwrite=false;
         int kountzero=0;
         int isize=SpectrumWaves.size();

         int isize2;

         double molspercm2, molspercm3;

         molspercm2=massgases[igx][nlay]*niceabund[molek][k]/(nucleon*AllMolW[molek][k]);
         molspercm3=molspercm2/(SlabThick[up_to_top-1-nlay]*100.0);


         // HITRAN line strength is in cm^2 per molecule
         // we need to convert again  --- SpectrumKays will now contain optical depth of each layer
         // No of molecules per cm^2=mass per cm^2/(mass 1 molecule)
         for(int iz=0; iz < isize; iz++){
            //convert to optical depth
            SpectrumKays[iz]=SpectrumKays[iz]*molspercm2;
             //massgases in g per cm^2, nucleon in grammes.
             // 100 slabthick = slab in centimetres
            if(outRef)SpectrumDeltaN[iz]=SpectrumDeltaN[iz]*molspercm3;
         }

         if(logplot){
         //strip out the zeros!
           for(int iz=0; iz < isize; iz++){
              if(SpectrumKays[iz]<=0.0)kountzero++;
         }}
         if(logplot){isize2=isize-kountzero;} else {isize2=isize;}

         // for PlotIt        
         int ntype=0;   //points=1   lines=0
         int ncol=1+nlay-NLStart;    //Give each line a different material colour
         int nstyle=0;  // solid or dippled
         int npoint=1;  //point size


         OutSpect.open(SpectFile.c_str(),ios::out | ios::app);
         TauSpect.open(SpectFileT.c_str(),ios::out | ios::app);
         if(outRef)RefSpect.open(SpectFileN.c_str(),ios::out | ios::app);
         if(outRef)SigSpect.open(SpectFileS.c_str(),ios::out | ios::app);
         if(!LayerChange){
         if(nlay==0)OutSpect << up_to_top << endl;
         if(nlay==0)TauSpect << up_to_top << endl;
         if(nlay==0)RefSpect << up_to_top << endl;
         if(nlay==0)SigSpect << up_to_top << endl;
         }
         else{
         if(nlay==NLStart)OutSpect << NLStop-NLStart << endl;
         if(nlay==NLStart)TauSpect << NLStop-NLStart << endl;
         if(nlay==NLStart)RefSpect << NLStop-NLStart << endl;
         if(nlay==NLStart)SigSpect << NLStop-NLStart << endl;
         }


          if(PlotIt){
            OutSpect << isize2 <<  "  " << ntype << "  " << ncol << "  " << nstyle << "  " << npoint << endl;
            TauSpect << isize2 <<  "  " << ntype << "  " << ncol << "  " << nstyle << "  " << npoint << endl;
            if(outRef)
            RefSpect << isize2 <<  "  " << ntype << "  " << ncol << "  " << nstyle << "  " << npoint << endl;
            if(outRef)
            SigSpect << isize2 <<  "  " << ntype << "  " << ncol << "  " << nstyle << "  " << npoint << endl;
          } 
          else{       
            OutSpect << isize2 << endl;
            TauSpect << isize2 << endl;
            if(outRef)RefSpect << isize << endl;
            if(outRef)SigSpect << isize << endl;
          } 

         double sigmatau; //optical depth per cm based on molecular absorption of gas
         double xwave;  //wavelength



         if(logplot){

            for(int i=0; i<isize; i++){
              if(SpectrumKays[i]>0){
                  //remember spectrum kays has been multiplied by mass per cm^2/mass per molecule
                  //convert back to cm^2 per molecule for output
                  OutSpect << setprecision(11) << SpectrumWaves[i] <<  "   " 
                           << log10(SpectrumKays[i]/molspercm2) << endl;
                  TauSpect << setprecision(11) << SpectrumWaves[i] <<  "   " 
                           << log10(SpectrumKays[i])<< endl;
                  if(outRef){

                  xwave=1./SpectrumWaves[i]; //wavelength in cm (cgs)

                  sigmatau=SpectrumDeltaN[i]/2.0/pi/molspercm3;
                           //so far we have (n^2-1)/2 pi N = atomic polarisability of 1 molecule
                           // next convert to scattering cross section in cm^2 per molecule
                  sigmatau=sigmatau*sigmatau*128*pow(pi,5.0)/
                           (3.0*pow(xwave,4.0));
                          //now convert to  scattering cross section for 1cm^3 of gas
                  sigmatau=sigmatau*molspercm3;


 
                  if(sigmatau>0){
                  SigSpect <<  setprecision(11) << SpectrumWaves[i] <<  "   " 
                           << log10(sigmatau) << endl;}
                 
                  }

              }  //endif  >0
             } //end for loop}
         }
         else{
           // not log plot
           for(int i=0; i<isize; i++){
              OutSpect << setprecision(11) << SpectrumWaves[i] <<  "   " 
                       << SpectrumKays[i]/molspercm2 << endl;
              TauSpect << setprecision(11) << SpectrumWaves[i] <<  "   " 
                       << SpectrumKays[i] << endl;
              if(outRef){

                  xwave=1./SpectrumWaves[i]; //wavelength in cm (cgs)

                  sigmatau=SpectrumDeltaN[i]/2.0/pi/molspercm3;
                           //so far we have (n^2-1)/2 pi N = atomic polarisability of 1 molecule
                           // next convert to scattering cross section in cm^2 per molecule
                  sigmatau=sigmatau*sigmatau*128*pow(pi,5.0)/
                           (3.0*pow(xwave,4.0));
                          //now convert to  scattering cross section for 1cm^3 of gas
                  sigmatau=sigmatau*molspercm3;

 
                  SigSpect <<  setprecision(11) << SpectrumWaves[i] <<  "   " 
                           << sigmatau << endl;}


           }

         }  //end for if logplot else
         //Ref index always has negatives, no log plot
           if(outRef){
           for(int i=0; i<isize; i++){
          //    if(fabs(SpectrumDeltaN[i])<2e-13){
              RefSpect << setprecision(11) << SpectrumWaves[i] <<  "   " << Speedlight*SpectrumDeltaN[i] << endl;
          //    }
           }}


           if(nlay==NLStop-1){

            if(PlotIt){
  
              OutSpect << 1 << "  " <<  1 << endl;
              OutSpect << "Wave@Number   @per_cm" << endl;
              OutSpect << 1 << "  " <<  1 << endl;  // this line for vertical axis

              TauSpect << 1 << "  " <<  1 << endl;
              TauSpect << "Wave@Number   @per_cm" << endl;
              TauSpect << 1 << "  " <<  0 << endl;  // this line for vertical axis

              if(outRef){
                RefSpect << 1 << "  " <<  1 << endl;
                RefSpect << "Wave@Number   @per_cm" << endl;
                RefSpect << 1 << "  " <<  0 << endl;  // this linefor vertical axis
                SigSpect << 1 << "  " <<  1 << endl;
                SigSpect << "Wave@Number   @per_cm" << endl;
                SigSpect << 1 << "  " <<  0 << endl;  // this linefor vertical axis
              }
              if(logplot){
                 OutSpect << "log10(Cross@Section)    @cm^2@Per@Mol" << endl;}
              else{
                 OutSpect << "Cross@Section   @cm^2@Per@Mol" << endl;
              }
              if(logplot){
                 TauSpect << "log10(Optical@Depth)" << endl;
                 if(outRef)SigSpect << "log10(scattering@optical@depth)    ." << endl;}
              else{
                 TauSpect << "Optical@Depth " << endl;
                 if(outRef)SigSpect << "scattering@optical@depth    ." << endl;
              }
              if(outRef)RefSpect << "(n-1)" << endl;

               if(LEGEND){
                 OutSpect << "1" << endl;
                 TauSpect << "1" << endl;
                 if(outRef)RefSpect << "1" << endl;
                 if(outRef)SigSpect << "1" << endl;

             for(int laynumber=0; laynumber< up_to_top; laynumber++){
                OutSpect << "Layer=@" << laynumber << endl;
                TauSpect << "Layer=@" << laynumber << endl;
                if(outRef)RefSpect << "Layer=@" << laynumber << endl;
                if(outRef)SigSpect << "Layer=@" << laynumber << endl;
             }} else {
                 OutSpect << "0" << endl;
                 TauSpect << "0" << endl;
                 if(outRef)RefSpect << "0" << endl;
                 if(outRef)SigSpect << "0" << endl;
             }
           
           } //endif for PlotIt}

         

           } //endif for nlay=up_to_top

           OutSpect.close();  
           TauSpect.close();    
           RefSpect.close();   
           SigSpect.close(); 
         
         } //endif for outspec

/*****************************************************************************************
        Work out k Distributions before releasing memory
******************************************************************************************/    

/*****************Release Memory**********************************************************/       
           SpectrumWaves.erase(SpectrumWaves.begin(), SpectrumWaves.end());
           SpectrumKays.erase(SpectrumKays.begin(), SpectrumKays.end());                         
           SpectrumDeltaN.erase(SpectrumDeltaN.begin(), SpectrumDeltaN.end());
           CentreWaves.erase(CentreWaves.begin(), CentreWaves.end());

           }  //end loop over layers nlay


          } // end loop over isotopologues k=0, k< nicemax




/****************************************************************************************
        Clear memory
******************************************************************************************/
        for(int i=0;i<vecsize;i++){
            cout << "gas number=" << molecule << " Isotopologue=" << i+1 << "  code=" <<nicecode[molek][i]
                 << "  number of lines=" << linecentres[i].size()<< endl;
                 linecentres[i].erase(linecentres[i].begin(), linecentres[i].end()); 
                 linestrengths[i].erase(linestrengths[i].begin(), linestrengths[i].end()); 
                 linegammaAir[i].erase(linegammaAir[i].begin(), linegammaAir[i].end()); 
                 linegammaSelf[i].erase(linegammaSelf[i].begin(), linegammaSelf[i].end());  
                 lineLSE[i].erase(lineLSE[i].begin(), lineLSE[i].end());       
                 lineTdep[i].erase(lineTdep[i].begin(), lineTdep[i].end());     
                 lineShift[i].erase(lineShift[i].begin(), lineShift[i].end());            
            }

          } // endif for UseParFile

      } //end loop over HITRAN par files igx=0 to ngas



         } //endif for calcspec

       cout << "The Value of ShiftLine was " << ShiftLine << endl;
     if(ozhh){
          OzWavesHH.erase( OzWavesHH.begin(), OzWavesHH.end() );
          for(int ioz=0; ioz <iozHH; ioz++){
             OzXHH[ioz].erase( OzXHH[ioz].begin(), OzXHH[ioz].end() );
          }
     }
     if(ozchap){

          for(int ioz=0; ioz < ichap; ioz++){
             OzWavesChap[ioz].erase( OzWavesChap[ioz].begin(), OzWavesChap[ioz].end() );
             OzXChap[ioz].erase( OzXChap[ioz].begin(), OzXChap[ioz].end() );
          }
     }


     if(rfun){delete[] lambda_resp; delete[] respval;}

  return 0;
}


