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

//toms library
void weightcoeff_(int*, double*, double*, double*, double*,
                       double*, double*);


}

void Questioner(int&, bool&, int&, int&, int&, int&, bool&, double&, double&,
                double&, double&, int&, bool&, bool&, int&, int&, int&, int&, bool&,
                int&, bool&, double&, double&, bool&, double&,
                bool&, double&, double&, double&, double&, double&, int&, int&,
                const char*);

void WaterVap(double* , double* , double* ,
              double* , double &, double &, double &, int &, int &, int &);
              

int main(int argc, char* argv[]){
/*
Stage 1: Define basic physical data, read in basic molecular atmosperic data
Stage 2: Calculate integration grid and interpolate data to integration grid
Stage 3: Integrate data to get slabs for all molecular masses, air density, and humidity
*/

  bool verbose;
/***************************************************************************************
   Begin Stage 1
****************************************************************************************
   Define basic physical data, read in basic atmosperic data
****************************************************************************************/

  // Four Aerosol layers, each to be split into nsplit slabs
  // Four layers are boundary layer, troposphere, stratosphere, and upper atmosphere
  int nsplit1,nsplit2,nsplit3,nsplit4;    
  // boundary layer visibility (or optical depth), toposphere visiblility
  // (or optical depth),  wavelength
  double visb,vist,lambda1, lambda2;
  bool vis;  //true if entered as visibility, false if entered as optical depth.
  // default atmosphere 1-6 or user defined
  int imode; //1 if Earth atmosphere, othere alternative not coded for yet.
  int iatm;
  bool switchR;
  bool switchA;
  int itypeu, itypes, itypet, itypeb;
  bool ocean,ground;
  int  ihumid;
  double  groundtemp, groundpress;
  bool cloud, know;
  bool aeroplane;
  double heightplane;

  bool default_pause;
  double HG, HB, HT, HS, HU;
  int ngaussL, ngaussH;
  const char* ReadInput;
  const char* RespFile;

   // First off we read in the user input, whether 
   // from command line (Dors1 -v) or an input file (Dors1)

    if(argc==1)
     verbose=false;
     ReadInput="Dors1in.dat";
     if(argc !=1)
     {
       if(argc==2){
         RespFile="myresp.dat";
         if(argv[1][0]=='-' && argv[1][1]=='v'){
           verbose=true;
           ReadInput="Dors1in.dat"; //Not used, but initialise anyway 
          }
         else
         {
             ReadInput=argv[1];
             cout << "Reading from file=" << ReadInput << endl;
         }
       }
        else
       {
         if(argc==3){
          RespFile=argv[2];
             cout << "Response function from file=" << RespFile<< endl;
         }
         else
         {
            cout << "Usage Dors1 or Dors1 -v" << endl;
            cout << "argc=" << argc << endl;
            exit(2);
         }
        }
    }

   /*****************************************************************/
   // NOTE -NO SAFETY BELTS, THE USER MUST MAKE SURE THE INPUT IS O.K.
   /*****************************************************************/

   Questioner(imode,verbose, nsplit1, nsplit2, nsplit3, nsplit4, vis, visb, vist,
   lambda1, lambda2, iatm, switchR, switchA,itypeu, itypes, itypet, itypeb, ocean,
   ihumid, know, groundtemp, groundpress, aeroplane, heightplane,
   default_pause, HG, HB, HT, HS, HU, ngaussL, ngaussH, ReadInput);  

   // for Gauss quadrature routine in toms library
    double Q[ngaussL],E[ngaussL],Xgauss[ngaussL],Wgauss[ngaussL],Work_g[9*ngaussL+8];
    double EPS, Xtemp[ngaussL];

    EPS=1e-15;
    for(int i=2;i<=ngaussL;i++){
       Q[i-1]=2.*i*i/(2.*i*(2.0*i-1.0));
       E[i-1]=2.*i*i/(2.*i*(2.0*i+1.0));
    }
    Q[0]=1.0;
    E[0]=1.0/3.0;
    //toms algorithm 125
    weightcoeff_(&ngaussL,Q,E,&EPS,Wgauss,Xgauss,Work_g);
    // why?

    //for some reason the x interval is [0,2]
   // and you have to double the weights if you want [0,2]
    //not only that the X are in reverse order
    // W not affected since it is symmetric anyway
    //now transform to [0,1] 

    for(int i=1; i <=ngaussL; i++){
       Xgauss[i-1]=Xgauss[i-1]/2.0;
       Xtemp[i-1]=Xgauss[i-1];
    } 
    for(int i=0; i <ngaussL; i++){
        Xgauss[i]=Xtemp[ngaussL-i-1];
    }
            
   /* *********************************************************** */
   /*        Quadrature Rule Done                                 */
   /* *********************************************************** */

   //Some constants

   double pi,abszero,R,AV,MWA,delta,nucleon;

   pi=acos(-1.0);               //Well pi obviusly
   abszero=273.15;              //absolte zero =-abszero celcius
   R=8.314472;                  //Universal gas constant
   AV=6.0221415e23;             //Avogadro's number;
   MWA=28.946;                  //Molecular Weight of Air
   delta=0.029;                 //depolarisation factor
   nucleon=1.660588E-24;         // "nucleon mass" in g


   ifstream fp_in;   //input file stream
   int gid[50];
   bool gmask[50];
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
   const char* AllMols[50];
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
   int nice[50]; // number of isotope variations (isotopolgues) in each  molecule
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
   double QTref[50][10];
   int geejay[50][10];
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
         cout << imol+1 << " " << nice[imol] << "  " << nicecode[imol][iso] << "  " << AllMolW[imol][iso] << endl;
      }
      if(!getline(fp_in,line1)){ break;}
       // << "  getline return\n";
      cout << line1;
     // if(getline(fp_in, line1).eof()){ break;}
      imol++;
   }
   fp_in.close(); 

   int ngas=0;
   for(int  i=0;i<50;i++){
       if(gmask[i])ngas++;}


   int gases[ngas];
   //MolFiles are just the Atmospheres/MolsX files actually used
   //rather than (but might possibly be) the whole set 
   const char* MolFiles[ngas];
   double MWGAS[ngas];
   ngas=0;
   for(int i=0;i<50;i++){
       if(gmask[i]){
           gases[ngas]=gid[i];
           MolFiles[ngas]=AllMols[i];
           //Assume natural abundances from molparams.txt
           MWGAS[ngas]=0.0;
           for(int j=0; j<nice[i]; j++){
              MWGAS[ngas]=MWGAS[ngas]+niceabund[i][j]*AllMolW[i][j];
           }
           ngas++;
           }}
    cout << "iatm=" << iatm << endl;
    for(int i=0; i<ngas;i++){
      cout << gases[i] << "  " << MolFiles[i] << "   MW=" << MWGAS[i] << endl;
          }
             
/* ********************************************************************* */
//            Points for quadrature
/* ********************************************************************* */


    double zstep1,zstep2,zstep3,zstep4,thik1,thik2,thik3,thik4;

    thik1=(HU-HS)*1000.0/((double)(nsplit1));
    thik2=(HS-HT)*1000.0/((double)(nsplit2));
    thik3=(HT-HB)*1000.0/((double)(nsplit3));
    thik4=(HB-HG)*1000.0/((double)(nsplit4));

   int nmix=3;      //for mixing aerosols
   int nprofs=6;    //Number of profiles stored

   fp_in.open("Atmospheres/Alts.dat", ios::in);
   if(!fp_in.is_open()){ cout << "Failed to open file (Atmospheres/Alts.dat)" << endl;  exit(1);}

   int i_altitudes;
   fp_in >> i_altitudes;
   double altitude[i_altitudes];
   for(int i=0;i<i_altitudes;i++){
     fp_in >> altitude[i];
     }
     //convert to metres
     for(int i=0; i<i_altitudes; i++)
        altitude[i]=altitude[i]*1000.0;
   fp_in.close();



   char sub[3];
   sprintf(sub,"%d",iatm);

   string Tempsname="Atmospheres/Temps .dat";
   // stringreplace(index, num1, const char*, num2)
   //index tells it where, up to  num1 characters repleced with up to num 2
   Tempsname.replace(17,1,sub,1);

   fp_in.open(Tempsname.c_str(), ios::in);
   if(!fp_in.is_open()){ cout << "Failed to open file, file= " <<  Tempsname.c_str() << endl;  exit(1);}

/*
   cout << Tempsname.c_str() << "  " <<  fp_in.good() << endl;
   cout << Tempsname.c_str() << "  " << fp_in.bad() << endl;
   cout << Tempsname.c_str() << "  " << fp_in.eof() << endl;
   cout << Tempsname.c_str() << "  " << fp_in.fail() << endl;  
*/
   int idata;    //LOWTRAN T-P-H20 heights
   fp_in >> idata;
   if(idata !=i_altitudes){
        cout << "Temperature data does not match altitudes\n";  exit(1);}
   double T[idata];

   for(int i=0;i<idata;i++){
      fp_in >>  T[i];
   }

   fp_in.close();

   string Pressname="Atmospheres/Press .dat";
   Pressname.replace(17,1,sub,1);
   fp_in.open(Pressname.c_str(), ios::in);
   if(!fp_in.is_open()){ cout << "Failed to open file, file= " <<  Pressname.c_str() << endl;  exit(1);}

   fp_in >> idata;
   if(idata !=i_altitudes){
         cout << " Pressure data does not match altitude data\n" << endl; exit(1);}

   double P[idata];
   for(int i=0;i<idata;i++){
      fp_in >>  P[i] ;
    //  cout << i << "  " << P[i] << endl;
      }
    fp_in.close();

    double Gee1[idata];  //Gee1 only needed if using Hum1_7.dat data
    double Gee2[idata];  //Gee2 only needed if using Hum2_7.dat data
    double Hum1[idata];  //Hum1 and Humids1 needed to store humidity values over liquid water
    double Hum2[idata];  //Hum2 and Humids2 needed to store humidity values over ice
    double GasPPM[idata][ngas];
    double x1;
    //now read in data for gases
    for(int i=0;i<ngas;i++){
       fp_in.open(MolFiles[i],ios::in);
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
       if(ihumid==2){
         ifstream fp_humid;
         //CHANGE FOR 1-6 +User  later
         fp_humid.open("Atmospheres/Hum1_U.dat", ios::in);
         fp_humid >> idataW1;
         if(idataW1 > i_altitudes){
              cout <<"Humidity Files don't match altitude\n"; exit(1);}
         //note that humidity files with idata=1 or 2 recommended
         //then GasPPM for H20 only changed at the ground or first km

         for(int i=0; i <idataW1; i++){ 
             fp_humid >> Hum1[i];}
         fp_humid.close();
         fp_humid.open("Atmospheres/Hum2_U.dat", ios::in);
         fp_humid >> idataW2;
         if(idataW1 != idataW2){
              cout <<"Humidity Files don't match\n"; exit(1);}
         for(int i=0; i <idataW1; i++){ 
             fp_humid >> Hum2[i];}
         fp_humid.close();
         iwi=1;
         WaterVap(P,T,Gee1,Hum1,AV,nucleon,MWGAS[0],idataW1,ihumid,iwi);
         iwi=2;
         WaterVap(P,T,Gee2,Hum2,AV,nucleon,MWGAS[0],idataW1,ihumid,iwi);

         for(int i=0; i <idataW1; i++){ 
             if(T[i]>abszero){  // remember abszero is the negative of absolute zero
                 GasPPM[i][ig]=Gee1[i];}
             else{
                 GasPPM[i][ig]=Gee2[i];}
        }
           //no longer want or need ihumid=2
           ihumid=1;
        }  // end for using Hum1 and Hum2 files and calculating GasPPM for H20
        else{
           // we calculate Hum1 and Hum2 from GassPPM for H20
           iwi=1;
           WaterVap(P,T,Gee1,Hum1,AV,nucleon,MWGAS[0],idata,ihumid,iwi);
           iwi=2;
           WaterVap(P,T,Gee2,Hum2,AV,nucleon,MWGAS[0],idata,ihumid,iwi);
        }
  
      if(ground==1){
            //scale temperature and pressure to ground;
            double tempP=groundpress/P[0];
            for(int i=0; i< idata; i++){
               P[i]=P[i]*tempP;  
              }
            double tempT=groundtemp/T[0];
            T[0]=T[0]*tempT;
            T[1]=T[1]*tempT;    // only scale lower atmosphere
           }


    // We need to integrate quanties in each sublayer, get gauss
    // quadrature weights and coefficients on (0,1)
/***************************************************************************************
   End Stage 1
****************************************************************************************
   Begin Stage 2, calc integration nodes and interpolate data  to them
****************************************************************************************/
        //need T and P data etc at fine resolution (heights)
        // Spline from T and P data etc at
        //declare P2 and T2 idata long for spline, store fine data
        // in Temp, Press, etc

        //prepare for calls to dbint4.f
        //
        double HBDRY[5], PBDRY[5];
        
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
                                 ,&ideriv, HBDRY+i, &inbv,WORK2);
         }

        //Need to find heights at equal pressure intervals
        //
        int up_to_top, up_to_up, up_to_strat, up_to_trop;
        up_to_top=nsplit4+nsplit3+nsplit2+nsplit1;
        up_to_up=nsplit4+nsplit3+nsplit2;
        up_to_strat=nsplit4+nsplit3;
        up_to_trop=nsplit4;
        double floorpress[up_to_top];
        double HeightFloors[up_to_top];
        double ReverseP[idata];
        double ReverseH[idata];
       //     cout  << "Reverse " << idata <<  endl;
        for(int i=0; i< idata; i++){
            ReverseP[i]=P[idata-1-i];
            ReverseH[i]=altitude[idata-1-i];}


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

        for(int i=0; i<up_to_top;i++)
           HeightFloors[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv,floorpress+i,&inbv,WORK2);


        ideriv=0,inbv=1; //inbv must always be set to 1 on first call to dbvalu
        dbint4_(altitude,T,&idata,&IBCL,&IBCR,&FBCL,&FBCR,
                       &KNOPT, TEE, BCOEF, &IN, &KORDER, WORK);

        double FloorTemp[up_to_top];
        double SlabThick[up_to_top];
        for(int i=0; i<up_to_top;i++){
           FloorTemp[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv,HeightFloors+i,&inbv,WORK2);

        //we shall do all calculations ground up
        // but output them top down
           if(i==up_to_top-1){
              SlabThick[i]=HBDRY[4]-HeightFloors[i];
            //  cout << "Floors " << HeightFloors[i] << " "  <<floorpress[i] << " " << FloorTemp[i] <<
            //   " " << SlabThick[i] << endl;
               }
              else{
              SlabThick[i]=HeightFloors[i+1]-HeightFloors[i];
            //  cout << "Floors " << HeightFloors[i] << " "  <<floorpress[i] << " " << FloorTemp[i] <<
            //   " " << SlabThick[i] << endl;
               }
         } //end i<up_to_top loop

        //Calc Gauss integration heights ground up
        int iheight;
        iheight=ngaussL*up_to_top;
        double heights[iheight];
        int kount=0;
        for(int i=0; i<nsplit4; i++){
           for(int j=0; j<ngaussL; j++){
           heights[kount]=HeightFloors[i]+SlabThick[i]*Xgauss[j];
   //        cout <<"kount=" << kount << "  heights=" << heights[kount] 
           kount++;
        }}
        for(int i=0; i<nsplit3; i++){
           for(int j=0; j<ngaussL; j++){
           heights[kount]=HeightFloors[up_to_trop+i]+SlabThick[up_to_trop+i]*Xgauss[j];
        //   cout <<"kount=" << kount << "  heights=" << heights[kount] << endl;
           kount++;
        }}
        for(int i=0; i<nsplit2; i++){
           for(int j=0; j<ngaussL; j++){
           heights[kount]=HeightFloors[up_to_strat+i]+SlabThick[up_to_strat+i]*Xgauss[j];
        //   cout <<"kount=" << kount << "  heights=" << heights[kount] << endl;
           kount++;
        }}
        for(int i=0; i<nsplit1; i++){
           for(int j=0; j<ngaussL; j++){
           heights[kount]=HeightFloors[up_to_up+i]+SlabThick[up_to_up+i]*Xgauss[j];
        //   cout <<"kount=" << kount << "  heights=" << heights[kount] << endl;
           kount++;
        }}
           
        double Temp[iheight],Press[iheight];
        for(int i=0; i<iheight;i++){
           Temp[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv,heights+i,&inbv,WORK);}
        ideriv=0,inbv=1; //inbv must always be set to 1 on first call to dbvalu
        dbint4_(altitude,P,&idata,&IBCL,&IBCR,&FBCL,&FBCR,
                       &KNOPT, TEE, BCOEF, &IN, &KORDER, WORK2);
        for(int i=0; i<iheight;i++)
           Press[i]=dbvalu_(TEE,BCOEF,&IN,&KORDER
                                 ,&ideriv,heights+i,&inbv,WORK2); 
        double ENN[iheight],Dens[iheight];
        for(int i=0; i<iheight;i++){
            //multiply by 100 to get pressure in Pascals
            ENN[i]=Press[i]/R/Temp[i]*AV*100.; 
            // multiply by MW to get mass per vol in gramms
            // divide by 1000 to get mass per vol in Kg
            // nucleon is just 1/Av 
            Dens[i]=ENN[i]*MWA*nucleon/1000.0;    
            }

        double GPPM[idata];
        double G[iheight];
        double Humids1[iheight];
        double Humids2[iheight];
        double Gas[iheight][ngas];   

        for(int n=0;n<ngas;n++){

           iwi=1;   // 1 if over water, 2 if over ice

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
           iwi=1;
           WaterVap(Press,Temp,G,Humids1,AV,nucleon,MWGAS[0],iheight,ihumid,iwi);
           iwi=2;
           WaterVap(Press,Temp,G,Humids2,AV,nucleon,MWGAS[0],iheight,ihumid,iwi);}
        }
        double DensGas[iheight][ngas];
        for( int n=0; n< ngas; n++){
           for( int i=0;i<iheight;i++){
               // divide ENN by 1E6 since Gas[i][n] is ppmv
               DensGas[i][n]=ENN[i]/1.0E6*Gas[i][n]*MWGAS[n]*nucleon/1000.0;
           }
         }
/***************************************************************************************
   End Stage 2
****************************************************************************************
   Begin Stage 3, Do Integrals
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
           for(int j=0; j < ngaussL; j++){
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
        cout << AvTemp[nsplit1-i] << "  " << AvPress[nsplit1-i] << "  " << totalmass <<endl;
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
           for(int j=0; j < ngaussL; j++){
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
        cout << AvTemp[nsplit2-i+idown] << "  " << AvPress[nsplit2-i+idown] << "  " << totalmass << endl;
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
           for(int j=0; j < ngaussL; j++){
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
        cout << AvTemp[nsplit3-i+idown] << "  " << AvPress[nsplit3-i+idown] << "  " << totalmass << endl;
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
           for(int j=0; j < ngaussL; j++){
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
        cout << AvTemp[i+idown] << "  " << AvPress[i+idown] << "  " << totalmass << endl;
        }

        cout << "mass of gas per square metre * g=" << 9.81*totalmass << "  Pascals="
              << 9.81*totalmass/100.0 << "mbar\n";
        cout << "Discrepency due (mostly) to interpolation\n";

        //convert from Kg per square metre to g per square cm.
        // Kg to g times=1000; metre square to cm square /=10000

        for(int k=0; k<ngas; k++){
          for(int i=0; i<up_to_top; i++){
             massgases[i][k]=massgases[i][k]/10.0;
          }
          totalmasses[k]=totalmasses[k]/10.0;
        }

        cout << "gases in each layer --- top down --- g per cm^2\n";
        for(int i=0; i<up_to_top; i++){ //top down actually
            for(int k=0; k < ngas ; k++){
                  cout << massgases[i][k] << " ";}
                  cout << endl;}
        cout << "totals\n";
        for(int k=0; k < ngas ; k++){
             cout << totalmasses[k] << " ";}
             cout <<  "  totalmass=" << totalmass <<endl; 

//    declaration for HITRAN input
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
      char flag;                                   //flag for poss, of line mixing
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
      double lambda;

      vector<double> linecentres[ngas]; //line centres in wave number units.
      

      string line;
      string entry[19];
      istringstream iss[19]; //19 input stringstreams
      for(int ig=0; ig<ngas; ig++){ //loop to read from HITRAN files.
        linecentres[ig].clear();
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
          linecentres[ig].push_back(nu);
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
      }
      //end loop over HITRAN par files
      for(int ig=0;ig<ngas;ig++){
         cout << linecentres[ig].size()<< endl;
      }

  return 0;
}


