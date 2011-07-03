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

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
using namespace::std;
//  replace this with hard wired version and use only with 
//  5 points or less --- want the combinatorics to check
//  out the hard wired stuff to be written later.
int Pickem(int** picks, int& m, int& N, int & ifacs, int &mcrit);

    double LagrangeInt( double* X, double* Y, int &N){

    double Integral;

    if(N==4){
     //four point Lagrange
     double coeff0, coeff1, coeff2, coeff3;
     double factor0, factor1, factor2, factor3;

     factor0=1.0/(X[0]-X[1])/(X[0]-X[2])/(X[0]-X[3]);
     factor1=1.0/(X[1]-X[0])/(X[1]-X[2])/(X[1]-X[3]);
     factor2=1.0/(X[2]-X[0])/(X[2]-X[1])/(X[2]-X[3]);
     factor3=1.0/(X[3]-X[0])/(X[3]-X[1])/(X[3]-X[2]);

     coeff3=Y[0]*factor0
           +Y[1]*factor1
           +Y[2]*factor2
           +Y[3]*factor3;  //coeff x^3


     coeff2=-Y[0]*factor0*(X[1]+X[2]+X[3])
            -Y[1]*factor1*(X[0]+X[2]+X[3])
            -Y[2]*factor2*(X[0]+X[1]+X[3])
            -Y[3]*factor3*(X[0]+X[1]+X[2]);     //coeff x^2

     coeff1=Y[0]*factor0*(X[1]*X[2]+X[1]*X[3]+X[2]*X[3])
           +Y[1]*factor1*(X[0]*X[2]+X[0]*X[3]+X[2]*X[3])
           +Y[2]*factor2*(X[0]*X[1]+X[0]*X[3]+X[1]*X[3])
           +Y[3]*factor3*(X[0]*X[1]+X[0]*X[2]+X[1]*X[2]); //coeff x


     coeff0=-Y[0]*factor0*X[1]*X[2]*X[3]
            -Y[1]*factor1*X[0]*X[2]*X[3]
            -Y[2]*factor2*X[0]*X[1]*X[3]
            -Y[3]*factor3*X[0]*X[1]*X[2];    //coeff x^0


  
     double X2, X0;
     double X2sq, X0sq;
     double X2cube, X0cube;
     double X2quart, X0quart;


     X2=X[N-1];               X0=X[0];
     X2sq=X2*X2;            X0sq=X0*X0;
     X2cube=X2*X2sq;        X0cube=X0*X0sq;
     X2quart=X2*X2cube;     X0quart=X0*X0cube;


     Integral=(X2quart-X0quart)/4.0*coeff3+(X2cube-X0cube)/3.0*coeff2
             +(X2sq-X0sq)/2.0*coeff1+(X2-X0)*coeff0;

    return Integral;
    }



     if(N==3){
     //three point Lagrange
     double coeff0, coeff1, coeff2, factor0, factor1, factor2;
     factor0=1.0/(X[0]-X[1])/(X[0]-X[2]);
     factor1=1.0/(X[1]-X[0])/(X[1]-X[2]);
     factor2=1.0/(X[2]-X[1])/(X[2]-X[0]);

     coeff2=Y[0]*factor0
           +Y[1]*factor1
           +Y[2]*factor2;  //coeff x^2


     coeff1=-Y[0]*factor0*(X[1]+X[2])
            -Y[1]*factor1*(X[0]+X[2])
            -Y[2]*factor2*(X[0]+X[1]);     //coeff x


     coeff0=Y[0]*factor0*X[1]*X[2]
             +Y[1]*factor1*X[0]*X[2]
             +Y[2]*factor2*X[1]*X[0];    //coeff x^0

     

     double X2, X0, X2sq, X0sq, X2cube, X0cube;
     X2=X[N-1]; X0=X[0];
     X2sq=X2*X2; X0sq=X0*X0; X2cube=X2*X2sq; X0cube=X0*X0sq;

     Integral=(X2cube-X0cube)/3.0*coeff2+(X2sq-X0sq)/2.0*coeff1+(X2-X0)*coeff0;
     return Integral;
     }


    if(N==2){
       //Two point Lagrange

       double coeff0, coeff1, factor0, factor1;
       factor0=1.0/(X[0]-X[1]);
       factor1=-factor0;

       coeff1=Y[0]*factor0
             +Y[1]*factor1;     //coeff x

       coeff0=-Y[0]*factor0*X[1]
              -Y[1]*factor1*X[0];    //coeff x^0

        double X2, X0, X2sq, X0sq;
        X2=X[1]; X0=X[0];
        X2sq=X2*X2; X0sq=X0*X0;

        Integral=(X2sq-X0sq)/2.0*coeff1+(X2-X0)*coeff0;

       return Integral;
    }    

    if(N>4){
       cout << "Use automatically generated coefficients to debug a case for N=" << N << endl;
       cout << "Don't try to use it for practical integration \n";
       cout << "Abandoning process! \n";
       exit(0);
    }
                  

    bool verbose=true;

    double factor[N];

     for(int i=0; i<N; i++){
       factor[i]=1.0;
       for(int j=0; j<N; j++){
         if(i !=j){
             factor[i]=factor[i]/(X[i]-X[j]);         
         }
       }
     } 

    double sumI[N][N];   //list of integrals for each Y;

    double mult=1.0;

    double coeffs[N+1];
   
    for(int i=0; i<N; i++){
      sumI[i][0]=1.0;   //coeff of x^N=1.0;
    }

    int ifacs;

    int iverbose=-1;    //get verbose output from Pickem and recursive calls
                        // to Pickem for a particular value of nobj --- set negative for not verbose

    int Nend;   // we only have to do half the work
    Nend=N/2;   // the second part of the list is the first half - except the numbers are
                // to be regarded as missing instead of present 
    int iremain=N%2;

    int factorials[Nend];

    for(int i=0; i<Nend; i++){
       ifacs=1;
       int nobj=i+1;
       // pick i+1 out of N objects
        
       for(int j=0; j<nobj; j++){
      //    simple to prove j+1 always goes exactly...
            ifacs=ifacs*(N-j)/(j+1);
       }
       if(verbose)cout << "i=" << i <<  "  ifacs=" << ifacs << endl;

       factorials[i]=ifacs;
 
       if(verbose)cout <<" declare picks " <<  ifacs << "  times " << nobj << endl;

       int** picks=new int*[ifacs];


       for(int k=0; k< ifacs; k++){
           picks[k]=new int[nobj];
       }

        // we only want <= m/2
       //if >m/2 have what we picked as missing from rather than present!
       int ireturn=Pickem(picks, nobj, N, ifacs, iverbose);

       bool mask[ifacs][N];
       int complement[ifacs][N];
       for(int k=0; k< ifacs; k++){
         for(int l=0; l<N; l++){
            complement[k][l]=l;
            mask[k][l]=true;
         }
       } 

      if(verbose)cout << " for n=" << nobj << " I can pick\n";

      bool addmult[N]; 

      for(int m=0; m< N; m++){
           sumI[m][nobj]=0.0;
          } 


      for(int k=0; k< ifacs; k++){
         mult=1.0;
         for(int m=0; m< N; m++){
           addmult[m]=true;
          } 
       
         for(int l=0; l<nobj; l++){

             if(verbose)cout << picks[k][l] << "  ";
             mask[k][picks[k][l]]=false;
      
             
             mult=mult*X[picks[k][l]];
             addmult[picks[k][l]]=false;

             } //end loop over l
         if(verbose)cout <<   endl;
         for(int m=0; m< N; m++){
            if(addmult[m])sumI[m][nobj]=sumI[m][nobj]+mult;
          } 
       } // end loop over k

       if(i<Nend-1){

          for(int m=0; m< N; m++){
             sumI[m][N-nobj]=0.0;
          } 


       if(verbose)cout << " for n=" << N-nobj << " I can pick\n";
       for(int k=0; k< ifacs; k++){
         mult=1.0;
         for(int m=0; m< N; m++){
           addmult[m]=true;
          } 

         for(int l=0; l<N; l++){
            if(mask[k][l]){
                   if(verbose)cout << complement[k][l] << "  ";
                   mult=mult*X[complement[k][l]];
                   addmult[complement[k][l]]=false;
            }

         } //end loop over l
         if(verbose)cout << endl;
         for(int m=0; m< N; m++){
            if(addmult[m])sumI[m][N-nobj]=sumI[m][N-nobj]+mult;
          } 
         } //end loop over k

       }else{


           if(iremain ==1){
             if(verbose)cout << "iremain for n=" << N-nobj << " I can pick\n";
               for(int m=0; m< N; m++){
                 sumI[m][N-nobj]=0.0;
                }


             for(int k=0; k< ifacs; k++){
                 mult=1.0;
                 for(int m=0; m< N; m++){
                     addmult[m]=true;
                 } 

               for(int l=0; l<N; l++){
                  if(mask[k][l]){
                    if(verbose)cout << complement[k][l] << "  ";

                      mult=mult*X[complement[k][l]];
                      addmult[complement[k][l]]=false;
                   } //endif  mask
                 } //end loop over l
             if(verbose)cout << endl;
                for(int m=0; m< N; m++){
                   if(addmult[m])sumI[m][N-nobj]=sumI[m][N-nobj]+mult;
                } 
             } //end loop over k



            }  //endif iremain
       }



       for(int k=0; k< ifacs; k++){
          delete[] picks[k];
       } 
       delete[] picks; 
       
   
     }  

     // integral from a to b
     double a=X[0]; double b=X[N-1];

     double limpow_a[N]; double limpow_b[N];
     for(int i=0; i<N; i++){
        limpow_a[i]=1.0; limpow_b[i]=1.0;
        int idiv=0;
        for(int j=i; j<N; j++){
          idiv++;
          limpow_a[i]=limpow_a[i]*a;
          limpow_b[i]=limpow_b[i]*b;
     //     cout << i << "  a= a* \n";
        }
    //    cout << "Divide by " << idiv << endl;
        limpow_a[i]=limpow_a[i]/((double)idiv);
        limpow_b[i]=limpow_b[i]/((double)idiv);
     } 
     
    

     double Integrals[N];
     
     for(int i=0; i<N; i++){
       Integrals[i]=0;
       for(int j=0; j< N; j++){       
        if(j%2==1){
          Integrals[i]=Integrals[i]-sumI[i][j]*(limpow_b[j]-limpow_a[j]);

        }
        else{
          Integrals[i]=Integrals[i]+sumI[i][j]*(limpow_b[j]-limpow_a[j]);
        }
      }
      Integrals[i]=Integrals[i]*factor[i];
     } 

     Integral=0;
     for(int i=0; i< N; i++){
        Integral=Integral+Integrals[i]*Y[i];
     }

     if(verbose)cout << "INTEGRAL=" << Integral << endl;

     if(N==4){
     //four point Lagrange
     double coeff0, coeff1, coeff2, coeff3;
     double factor0, factor1, factor2, factor3;

     factor0=1.0/(X[0]-X[1])/(X[0]-X[2])/(X[0]-X[3]);
     factor1=1.0/(X[1]-X[0])/(X[1]-X[2])/(X[1]-X[3]);
     factor2=1.0/(X[2]-X[0])/(X[2]-X[1])/(X[2]-X[3]);
     factor3=1.0/(X[3]-X[0])/(X[3]-X[1])/(X[3]-X[2]);

     coeff3=Y[0]*factor0
           +Y[1]*factor1
           +Y[2]*factor2
           +Y[3]*factor3;  //coeff x^3


     coeff2=-Y[0]*factor0*(X[1]+X[2]+X[3])
            -Y[1]*factor1*(X[0]+X[2]+X[3])
            -Y[2]*factor2*(X[0]+X[1]+X[3])
            -Y[3]*factor3*(X[0]+X[1]+X[2]);     //coeff x^2

     coeff1=Y[0]*factor0*(X[1]*X[2]+X[1]*X[3]+X[2]*X[3])
           +Y[1]*factor1*(X[0]*X[2]+X[0]*X[3]+X[2]*X[3])
           +Y[2]*factor2*(X[0]*X[1]+X[0]*X[3]+X[1]*X[3])
           +Y[3]*factor3*(X[0]*X[1]+X[0]*X[2]+X[1]*X[2]); //coeff x


     coeff0=-Y[0]*factor0*X[1]*X[2]*X[3]
            -Y[1]*factor1*X[0]*X[2]*X[3]
            -Y[2]*factor2*X[0]*X[1]*X[3]
            -Y[3]*factor3*X[0]*X[1]*X[2];    //coeff x^0


  
     double X2, X0;
     double X2sq, X0sq;
     double X2cube, X0cube;
     double X2quart, X0quart;


     X2=X[N-1];               X0=X[0];
     X2sq=X2*X2;            X0sq=X0*X0;
     X2cube=X2*X2sq;        X0cube=X0*X0sq;
     X2quart=X2*X2cube;     X0quart=X0*X0cube;


     Integral=(X2quart-X0quart)/4.0*coeff3+(X2cube-X0cube)/3.0*coeff2
             +(X2sq-X0sq)/2.0*coeff1+(X2-X0)*coeff0;
     }

    return Integral;

}
int Pickem(int** picks, int& m, int & N, int & ifacs, int &mcrit){

    //m is the nuber of objects to picked out of N obects

    int j=m-1;


    if(m==1){
    for(int i1=0; i1<ifacs; i1++){
       for(int j=0; j<m; j++){
         picks[i1][j]=i1;
         if(mcrit==4)cout << picks[i1][j] << endl;
    }}}
    if(m==1)return 0;

    //we are picking more than one object out of N




    int nextN=N-1;
    int nextm=m-1;
    int kount=0;
    int start=0;
    int lastkount;

    while(kount<ifacs){

       for(int j1=0; j1<nextm; j1++){


       lastkount=kount;

       int nextNend=nextN/2;
       int nextiremain=nextNend%2;
       nextNend=nextNend+nextiremain;


       int nextifacs=1;
        
       for(int j1=0; j1<nextm; j1++){
            nextifacs=nextifacs*(nextN-j1)/(j1+1);
       }


       int** nextpicks=new int*[nextifacs];

       if(mcrit==4)cout <<" declare picks " <<  nextifacs << "  times " << nextm << endl;
       for(int k=0; k< nextifacs; k++){
           nextpicks[k]=new int[nextm];
       }


       int ireturn=Pickem(nextpicks, nextm, nextN, nextifacs, mcrit);

       if(mcrit==4)cout << "nextpicks picked\n";
       for(int k=0; k< nextifacs; k++){
         for(int l=0; l<nextm; l++){
            if(mcrit==4)cout << nextpicks[k][l] << "  ";
         }
        if(mcrit==4)cout << "Then " << endl;
       }


    
       for(int k=0; k<nextifacs; k++){
          

          if(mcrit==4)cout << "picks picked\n";


          picks[k+lastkount][0]=start;
          if(mcrit==4)cout << picks[k+lastkount][0] << " ";

          
          for(int l=0; l<nextm; l++){
             picks[k+lastkount][l+1]=nextpicks[k][l]+start+1;
             if(mcrit==4)cout << picks[k+lastkount][l+1] << " ";

          } 
          kount++;     
          if(mcrit==4)cout << endl;
       }
       start++;

       if(mcrit==4)cout << " Then picked " << start << "  as the new start\n";

       nextN--;

       if(mcrit==4){
          cout << " debug point " << endl;
       }

       lastkount=kount;


       

       for(int k=0; k< nextifacs; k++){
          delete[] nextpicks[k];
       } 
       delete[] nextpicks; 
       
       }

   } // end while kount
    
}

