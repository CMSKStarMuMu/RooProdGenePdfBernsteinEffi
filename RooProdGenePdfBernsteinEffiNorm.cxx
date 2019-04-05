/***************************************************************************** 
 * Project: RooProdGenePdfBernsteinEffiNorm                                                 * 
 *                                                                           * 
 * P.Dini fecit, Anno Domini MMXVIII                                         *
 * "Est modus in rebus."                                                     * 
 *                                                                           * 
 * Class to describe 3D angular efficiency in B0->K*MuMU Analysis            * 
 *                                                                           * 
 *****************************************************************************/ 


#include "Riostream.h" 

#include "RooProdGenePdfBernsteinEffiNorm.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h"
#include "RooRealProxy.h" 
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include <math.h> 
#include "TMath.h"  
#include <algorithm>  

ClassImp(RooProdGenePdfBernsteinEffiNorm); 
   
  
   RooProdGenePdfBernsteinEffiNorm::RooProdGenePdfBernsteinEffiNorm(const char *name, const char *title, 
                        RooAbsReal& x, RooAbsReal& y, RooAbsReal& z, 
			const RooArgList& coefList,
			int maxDegree1, int maxDegree2, int maxDegree3,int ilNorm
			) :
   RooAbsReal(name,title), 
   _x("_x","_x",this,x),
   _y("_y","_y",this,y),
   _z("_z","_z",this,z), 
   _coefList("_coefList","_coefList",this),
   _maxDegree1(maxDegree1),
   _maxDegree2(maxDegree2),
   _maxDegree3(maxDegree3),
   _ilNorm(ilNorm)
{ 
    TIterator* coefIter1 = coefList.createIterator() ;
     RooAbsArg* coef1 ;
     while((coef1 = (RooAbsArg*)coefIter1->Next())) {
       if (!dynamic_cast<RooAbsReal*>(coef1)) {
    	 std::cout << "RooProdGenePdfBernsteinEffiNorm::ctor(" << GetName() << ") ERROR: coefficient " << coef1->GetName()
    	      << " is not of type RooAbsReal" << std::endl ;
    	 assert(0) ;
       }
       _coefList.add(*coef1) ;
//   	 std::cout << "RooProdGenePdfBernsteinEffiNorm::ctor(" << GetName() << ") coefficient " << coef->GetName()<<std::endl ;
      }
     delete coefIter1;
     _maxDegreeV = std::max(maxDegree3,std::max(maxDegree1,maxDegree2));
     printf("RooBernstein: Max(Numbers of degree) = %d\n",_maxDegreeV);
     printf("RooBernstein: Integral Number = %d\n",ilNorm);
 } 
//==============================================================================================
RooProdGenePdfBernsteinEffiNorm::RooProdGenePdfBernsteinEffiNorm(const RooProdGenePdfBernsteinEffiNorm& other, const char* name) :  
   RooAbsReal(other,name), 
   _x("_x",this,other._x),
   _y("_y",this,other._y),
   _z("_z",this,other._z), 
   _coefList("_coefList",this,other._coefList),
   _maxDegree1(other._maxDegree1),
   _maxDegree2(other._maxDegree2),
   _maxDegree3(other._maxDegree3),
   _maxDegreeV(other._maxDegreeV),
   _ilNorm(other._ilNorm)
{ 
} 
//==============================================================================================
Double_t RooProdGenePdfBernsteinEffiNorm::evaluate() const 
{ 
       double x    = _x; // x, y, z...
       double y    = _y; // x, y, z...
       double z    = _z; // x, y, z...

       double cosxq =	  x*x;
       double cosyq =	  y*y;
       double sinxq = (1.-cosxq);
       double sinyq = (1.-cosyq);
       double cos2x = (2.*cosxq-1.);
       double sin2x =  2.*x*sqrt(sinxq);
       double sin2y =  2.*y*sqrt(sinyq);
       double sinx  =  sqrt(sinxq);
//
//
       fptype   pdf  = 9.0/(32.0*TMath::Pi());
      if     (_ilNorm==0){
       pdf *= (3./4.*sinyq);
      }else if(_ilNorm==1){ 
       pdf *= cosyq;
      }else if(_ilNorm==2){  
       pdf *= -1./4.*sinyq*cos2x;
      }else if(_ilNorm==3){  
       pdf *= -cosyq*cos2x;
       
      }else if(_ilNorm==4){  
       pdf *= 0.5*sinyq*sinxq*cos(2.*z);
      }else if(_ilNorm==5){ 
       pdf *= -cos(z)*(0.5*sin2y*sin2x);
      }else if(_ilNorm==6){ 
       pdf *= -cos(z)*sin2y*sinx;
      }else if(_ilNorm==7){ 
       pdf *= -sin(z)*sin2y*sinx;
      }else if(_ilNorm==8){
       pdf *= -0.5*sin(z)*sin2y*sin2x ;
      }else if(_ilNorm==9){
       pdf *= 2*sinyq*x;
      }else if(_ilNorm==10){
       pdf *= -sinyq*sinxq*sin(2.*z);
      }else{
        return -10000.;
       }

//        fptype   pdf  = 9.0/(32.0*TMath::Pi());
//        if     (_ilNorm==0){
//         pdf *= 3./4.*(1.0-y)*(1.0+y);
//        }else if(_ilNorm==1){
//         pdf *= y*y;
//        }else if(_ilNorm==2){
//         pdf *= -1./4.*(1.0-y)*(1.0+y)*(2.*x*x-1.);
//        }else if(_ilNorm==3){
//         pdf *= -y*y*(2.*x*x-1.);
//        }else if(_ilNorm==4){
//         pdf *= 0.5*(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Cos(2.0*z);
//        }else if(_ilNorm==5){
//         pdf *= -TMath::Cos(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
//        }else if(_ilNorm==6){
//         pdf *= -TMath::Cos(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
//        }else if(_ilNorm==7){
//         pdf *= -TMath::Sin(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
//        }else if(_ilNorm==8){
//         pdf *= -TMath::Sin(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)) ;
//        }else if(_ilNorm==9){
//         pdf *= 2.0*(1.-y*y)*x;
//        }else if(_ilNorm==10){
//         pdf *= -(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Sin(2.0*z);
//        }else{
//         return -10000.;
//        }
//






/*  
      double F_L = 0.5;
     double P4p = 0.0000e+00;
     double P5p = 0.00000e+00;
     double P6p = 0.00000e+00;
     double P8p = 1.00000e+00;
     double P_1 = 0.;
     double P_2 =0.;
     double P_3 = 0.00000e+00;


  double F_T   = (1.-F_L);

  double cosxq =     x*x;
  double cosyq =     y*y;
  double sinxq = (1.-cosxq);
  double sinyq = (1.-cosyq);
  double cos2x = (2.*cosxq-1.);
     double sin2x =  2.*x*sqrt(sinxq);
     double sin2y =  2.*y*sqrt(sinyq);
//double sin2x =  2.*x*TMath::Sin(TMath::ACos(x));
//double sin2y =  2.*y*TMath::Sin(TMath::ACos(y));
//  double sinx  =  TMath::Sin(TMath::ACos(x));
  double sinx  =  sqrt(sinxq);
  

//folded  
//   double pdf = 9./(8.*TMath::Pi())*(2./3.*(F_S+A_S*y)*(1-x*x)+A5S*sqrt(1-y*y)*sqrt(1-x*x)*TMath::Cos(z))+ 
//   (1.-F_S)*(2.*F_L*y*y*(1-x*x)+0.5*(1.-F_L)*(1-y*y)*(1.+x*x)+0.5*P_1*(1.-F_L)*(1-y*y)*(1-x*x)*TMath::Cos(2*z)+ 
//    2*P5p*y*sqrt(F_L*(1.-F_L))*sqrt(1-y*y)*sqrt(1-x*x)*TMath::Cos(z));
  
  double pdf = 9./(32.*TMath::Pi())*(3./4.*F_T*sinyq+F_L*cosyq+\
  (1./4.*F_T*sinyq-F_L*cosyq)*cos2x+0.5*P_1*F_T*sinyq*sinxq*TMath::Cos(2.*z)+\
  sqrt(F_L*F_T)*TMath::Cos(z)*(0.5*P4p*sin2y*sin2x+P5p*sin2y*sinx)-\
  sqrt(F_L*F_T)*TMath::Sin(z)*(P6p*sin2y*sinx-0.5*P8p*sin2y*sin2x)+\
  2.*P_2*F_T*sinyq*x-P_3*F_T*sinyq*sinxq*TMath::Sin(2.*z));
 */ 
//        double Norm[11];
//        
//        fptype   pdfCost  = 9.0/(32.0*TMath::Pi());
// 
 //        Norm[0] = 3./4.*(1.0-y)*(1.0+y);
//         Norm[1] = y*y;
//         Norm[2] = -1./4.*(1.0-y)*(1.0+y)*(2.*x*x-1.);
//         Norm[3] = -y*y*(2.*x*x-1.);
//         Norm[4] = 0.5*(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Cos(2.0*z);
//         Norm[5] = -TMath::Cos(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
//         Norm[6] = -TMath::Cos(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
//         Norm[7] = -TMath::Sin(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
//         Norm[8] = -TMath::Sin(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)) ;
//         Norm[9] = 2.0*(1.-y*y)*x;
//         Norm[10] = -(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Sin(2.0*z);

// Norm[0]  =    (3./4.*sinyq);
// Norm[1]  =   cosyq;
// Norm[2]  =   -1./4.*sinyq*cos2x;
// Norm[3]  =   -cosyq*cos2x;
// Norm[4]  =   0.5*sinyq*sinxq*cos(2.*z);
// Norm[5]  =   -cos(z)*(0.5*sin2y*sin2x);
// Norm[6]  =   -cos(z)*sin2y*sinx;
// Norm[7]  =   -sin(z)*sin2y*sinx;
// Norm[8]  =   -0.5*sin(z)*sin2y*sin2x ;
// Norm[9]  =   2*sinyq*x;
// Norm[10] =   -sinyq*sinxq*sin(2.*z);
//        
//      double pdf = (F_T*Norm[0] + F_L*Norm[1] +\
//      F_T*(-Norm[2]) + F_L*Norm[3] + P_1*F_T*Norm[4] +\
//      sqrt(F_L*F_T)*(P4p*(-Norm[5]) + P5p*(-Norm[6]))-\
//      sqrt(F_L*F_T)*(P6p*(-Norm[7]) - P8p*(-Norm[8]))+\
//      P_2*F_T*Norm[9] - P_3*F_T*(-Norm[10]));
//      
//        pdf = Norm[8]*pdfCost;
       
       
       double xmin = _x.min();
       double xdif = _x.max()-xmin;
              x    =(_x-xmin)/xdif;
       double ymin = _y.min();
       double ydif = _y.max()-ymin;
              y    =(_y-ymin)/ydif;
       double zmin = _z.min();
       double zdif = _z.max()-zmin;
              z    =(_z-zmin)/zdif;
       
       double sx[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   
       double sy[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
   
       double sz[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
       
       
       sx[0]=1.0;
       sy[0]=1.0;
       sz[0]=1.0;
       for( int i = 1; i <= _maxDegreeV ; ++i){
        sx[i]= sx[i-1]*(1.-x);
        sy[i]= sy[i-1]*(1.-y);
        sz[i]= sz[i-1]*(1.-z);
       }
       
//        fptype tx = 1.0;
//        fptype ty = 1.0;
//        fptype tz = 1.0;
//        for( unsigned int i = 0; i <=6 ; ++i){
//         int ii = 6-i;
//         
// //       for( int i = 1; i <= max(maxDegree3,max(maxDegree1,maxDegree2)) ; ++i){
//         if(ii<=_maxDegree1) {sx[ii]*= tx*device_coeffbinomial(_maxDegree1,ii);tx*=x;}
//         if(ii<=_maxDegree2) {sy[ii]*= ty*device_coeffbinomial(_maxDegree2,ii);ty*=y;}
//         if(ii<=_maxDegree3) {sz[ii]*= tz*device_coeffbinomial(_maxDegree3,ii);tz*=z;}
// 	
//        }
   
// 
       double bernknvalx = 0.;
       double bernknvaly = 0.;
       double bernknvalz = 0.;
       
       
       int ipar =0;
       double func =0.0;
       
       double tx = 1.;
       for(int i = 0; i <= _maxDegree1 ; ++i) {
         bernknvalx =  device_coeffbinomial(_maxDegree1,i)*tx*sx[_maxDegree1-i];
       
         double ty = 1.;
         for(int j = 0; j <= _maxDegree2 ; ++j) {
	  bernknvaly =  device_coeffbinomial(_maxDegree2,j)*ty*sy[_maxDegree2-j];
	 
          double tz = 1.;
          for(int k = 0; k <= _maxDegree3 ; ++k) {
	   bernknvalz =  device_coeffbinomial(_maxDegree3,k)*tz*sz[_maxDegree3-k];
//	   func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//           func += ((RooAbsReal&) _coefList[ipar]).getVal()*sx[_maxDegree1-i]*sy[_maxDegree2-j]*sz[_maxDegree3-k];
	   func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) _coefList[ipar]).getVal()<<std::endl;
	   ipar++;
	   tz *= z;
	  }
	   ty *= y;
         }
	   tx *= x;
       }
//        std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
//        std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
//        std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
//      if(func<1.E-30)  std::cout<<"evaluate = "<<func<<" x,y,z ="<<_x<<" "<<_y<<_z<<" "<<std::endl;
 
       if(func<1.E-30) func=1.E-30; 
//       func = func/(xdif*ydif*zdif); 
       func=func*pdf;
//        std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
//        std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
//        std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
//       exit(0);
//      func=func/(1.+maxDegree1)/(1.+maxDegree2)/(1.+maxDegree3);
//      std::cout<<"eval = "<<func/(xdif*ydif*zdif)<<std::endl;
      return  func;
//      return  func/(xdif*ydif*zdif);
} 
//==============================================================================================
// Double_t RooProdGenePdfBernsteinEffiNorm::evaluate() const 
// { 
//        double x    = _x; // x, y, z...
//        double y    = _y; // x, y, z...
//        double z    = _z; // x, y, z...
// 
//        double cosxq =	  x*x;
//        double cosyq =	  y*y;
//        double sinxq = (1.-cosxq);
//        double sinyq = (1.-cosyq);
//        double cos2x = (2.*cosxq-1.);
//        double sin2x =  2.*x*sqrt(sinxq);
//        double sin2y =  2.*y*sqrt(sinyq);
//        double sinx  =  sqrt(sinxq);
// //
// //
//        fptype   pdf  = 9.0/(32.0*TMath::Pi());
//       if     (_ilNorm==0){
//        pdf *= (3./4.*sinyq);
//       }else if(_ilNorm==1){ 
//        pdf *= cosyq;
//       }else if(_ilNorm==2){  
//        pdf *= -1./4.*sinyq*cos2x;
//       }else if(_ilNorm==3){  
//        pdf *= -cosyq*cos2x;
//        
//       }else if(_ilNorm==4){  
//        pdf *= 0.5*sinyq*sinxq*cos(2.*z);
//       }else if(_ilNorm==5){ 
//        pdf *= -cos(z)*(0.5*sin2y*sin2x);
//       }else if(_ilNorm==6){ 
//        pdf *= -cos(z)*sin2y*sinx;
//       }else if(_ilNorm==7){ 
//        pdf *= -sin(z)*sin2y*sinx;
//       }else if(_ilNorm==8){
//        pdf *= -0.5*sin(z)*sin2y*sin2x ;
//       }else if(_ilNorm==9){
//        pdf *= 2*sinyq*x;
//       }else if(_ilNorm==10){
//        pdf *= -sinyq*sinxq*sin(2.*z);
//       }else{
//         return -10000.;
//        }
// 
// //        fptype   pdf  = 9.0/(32.0*TMath::Pi());
// //        if     (_ilNorm==0){
// //         pdf *= 3./4.*(1.0-y)*(1.0+y);
// //        }else if(_ilNorm==1){
// //         pdf *= y*y;
// //        }else if(_ilNorm==2){
// //         pdf *= -1./4.*(1.0-y)*(1.0+y)*(2.*x*x-1.);
// //        }else if(_ilNorm==3){
// //         pdf *= -y*y*(2.*x*x-1.);
// //        }else if(_ilNorm==4){
// //         pdf *= 0.5*(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Cos(2.0*z);
// //        }else if(_ilNorm==5){
// //         pdf *= -TMath::Cos(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
// //        }else if(_ilNorm==6){
// //         pdf *= -TMath::Cos(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
// //        }else if(_ilNorm==7){
// //         pdf *= -TMath::Sin(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
// //        }else if(_ilNorm==8){
// //         pdf *= -TMath::Sin(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)) ;
// //        }else if(_ilNorm==9){
// //         pdf *= 2.0*(1.-y*y)*x;
// //        }else if(_ilNorm==10){
// //         pdf *= -(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Sin(2.0*z);
// //        }else{
// //         return -10000.;
// //        }
// //
// 
// 
// 
// 
// 
// 
// /*  
//       double F_L = 0.5;
//      double P4p = 0.0000e+00;
//      double P5p = 0.00000e+00;
//      double P6p = 0.00000e+00;
//      double P8p = 1.00000e+00;
//      double P_1 = 0.;
//      double P_2 =0.;
//      double P_3 = 0.00000e+00;
// 
// 
//   double F_T   = (1.-F_L);
// 
//   double cosxq =     x*x;
//   double cosyq =     y*y;
//   double sinxq = (1.-cosxq);
//   double sinyq = (1.-cosyq);
//   double cos2x = (2.*cosxq-1.);
//      double sin2x =  2.*x*sqrt(sinxq);
//      double sin2y =  2.*y*sqrt(sinyq);
// //double sin2x =  2.*x*TMath::Sin(TMath::ACos(x));
// //double sin2y =  2.*y*TMath::Sin(TMath::ACos(y));
// //  double sinx  =  TMath::Sin(TMath::ACos(x));
//   double sinx  =  sqrt(sinxq);
//   
// 
// //folded  
// //   double pdf = 9./(8.*TMath::Pi())*(2./3.*(F_S+A_S*y)*(1-x*x)+A5S*sqrt(1-y*y)*sqrt(1-x*x)*TMath::Cos(z))+ 
// //   (1.-F_S)*(2.*F_L*y*y*(1-x*x)+0.5*(1.-F_L)*(1-y*y)*(1.+x*x)+0.5*P_1*(1.-F_L)*(1-y*y)*(1-x*x)*TMath::Cos(2*z)+ 
// //    2*P5p*y*sqrt(F_L*(1.-F_L))*sqrt(1-y*y)*sqrt(1-x*x)*TMath::Cos(z));
//   
//   double pdf = 9./(32.*TMath::Pi())*(3./4.*F_T*sinyq+F_L*cosyq+\
//   (1./4.*F_T*sinyq-F_L*cosyq)*cos2x+0.5*P_1*F_T*sinyq*sinxq*TMath::Cos(2.*z)+\
//   sqrt(F_L*F_T)*TMath::Cos(z)*(0.5*P4p*sin2y*sin2x+P5p*sin2y*sinx)-\
//   sqrt(F_L*F_T)*TMath::Sin(z)*(P6p*sin2y*sinx-0.5*P8p*sin2y*sin2x)+\
//   2.*P_2*F_T*sinyq*x-P_3*F_T*sinyq*sinxq*TMath::Sin(2.*z));
//  */ 
// //        double Norm[11];
// //        
// //        fptype   pdfCost  = 9.0/(32.0*TMath::Pi());
// // 
//  //        Norm[0] = 3./4.*(1.0-y)*(1.0+y);
// //         Norm[1] = y*y;
// //         Norm[2] = -1./4.*(1.0-y)*(1.0+y)*(2.*x*x-1.);
// //         Norm[3] = -y*y*(2.*x*x-1.);
// //         Norm[4] = 0.5*(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Cos(2.0*z);
// //         Norm[5] = -TMath::Cos(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
// //         Norm[6] = -TMath::Cos(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
// //         Norm[7] = -TMath::Sin(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
// //         Norm[8] = -TMath::Sin(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)) ;
// //         Norm[9] = 2.0*(1.-y*y)*x;
// //         Norm[10] = -(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Sin(2.0*z);
// 
// // Norm[0]  =    (3./4.*sinyq);
// // Norm[1]  =   cosyq;
// // Norm[2]  =   -1./4.*sinyq*cos2x;
// // Norm[3]  =   -cosyq*cos2x;
// // Norm[4]  =   0.5*sinyq*sinxq*cos(2.*z);
// // Norm[5]  =   -cos(z)*(0.5*sin2y*sin2x);
// // Norm[6]  =   -cos(z)*sin2y*sinx;
// // Norm[7]  =   -sin(z)*sin2y*sinx;
// // Norm[8]  =   -0.5*sin(z)*sin2y*sin2x ;
// // Norm[9]  =   2*sinyq*x;
// // Norm[10] =   -sinyq*sinxq*sin(2.*z);
// //        
// //      double pdf = (F_T*Norm[0] + F_L*Norm[1] +\
// //      F_T*(-Norm[2]) + F_L*Norm[3] + P_1*F_T*Norm[4] +\
// //      sqrt(F_L*F_T)*(P4p*(-Norm[5]) + P5p*(-Norm[6]))-\
// //      sqrt(F_L*F_T)*(P6p*(-Norm[7]) - P8p*(-Norm[8]))+\
// //      P_2*F_T*Norm[9] - P_3*F_T*(-Norm[10]));
// //      
// //        pdf = Norm[8]*pdfCost;
//        
//        
// //        double xmin = _x.min();
// //        double xdif = _x.max()-xmin;
// //               x    =(_x-xmin)/xdif;
// //        double ymin = _y.min();
// //        double ydif = _y.max()-ymin;
// //               y    =(_y-ymin)/ydif;
// //        double zmin = _z.min();
// //        double zdif = _z.max()-zmin;
// //               z    =(_z-zmin)/zdif;
// //        
// //        double sx[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
// //    
// //        double sy[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
// //    
// //        double sz[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
// //        
// //        
// //        sx[0]=1.0;
// //        sy[0]=1.0;
// //        sz[0]=1.0;
// //        for( int i = 1; i <= _maxDegreeV ; ++i){
// //         sx[i]= sx[i-1]*(1.-x);
// //         sy[i]= sy[i-1]*(1.-y);
// //         sz[i]= sz[i-1]*(1.-z);
// //        }
// //        
// //        fptype tx = 1.0;
// //        fptype ty = 1.0;
// //        fptype tz = 1.0;
// //        for( unsigned int i = 0; i <=6 ; ++i){
// //         int ii = 6-i;
// //         
// // //       for( int i = 1; i <= max(maxDegree3,max(maxDegree1,maxDegree2)) ; ++i){
// //         if(ii<=_maxDegree1) {sx[ii]*= tx*device_coeffbinomial(_maxDegree1,ii);tx*=x;}
// //         if(ii<=_maxDegree2) {sy[ii]*= ty*device_coeffbinomial(_maxDegree2,ii);ty*=y;}
// //         if(ii<=_maxDegree3) {sz[ii]*= tz*device_coeffbinomial(_maxDegree3,ii);tz*=z;}
// // 	
// //        }
// //    
// // // 
// // //        double bernknvalx = 0.;
// // //        double bernknvaly = 0.;
// // //        double bernknvalz = 0.;
// //        
// //        
// //        int ipar =0;
// //        double func =0.0;
// //        
// // //       double tx = 1.;
// // //       long double sx = sx_ini;
// //        for(int i = 0; i <= _maxDegree1 ; ++i) {
// // //         bernknvalx =  device_coeffbinomial(_maxDegree1,i)*tx*sx[_maxDegree1-i];
// //        
// // //         double ty = 1.;
// // //         long double sy =sy_ini;
// //          for(int j = 0; j <= _maxDegree2 ; ++j) {
// // //	  bernknvaly =  device_coeffbinomial(_maxDegree2,j)*ty*sy[_maxDegree2-j];
// // 	 
// // //          double tz = 1.;
// // //          long double sz = sz_ini;
// //           for(int k = 0; k <= _maxDegree3 ; ++k) {
// // //	   bernknvalz =  device_coeffbinomial(_maxDegree3,k)*tz*sz[_maxDegree3-k];
// // //	   func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
// //            func += ((RooAbsReal&) _coefList[ipar]).getVal()*sx[_maxDegree1-i]*sy[_maxDegree2-j]*sz[_maxDegree3-k];
// // //	   func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
// // //	   std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) _coefList[ipar]).getVal()<<std::endl;
// // 	   ipar++;
// // //	   tz *= z;
// // 	  }
// // //	   ty *= y;
// //          }
// // //	   tx *= x;
// //        }
// // //        std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
// // //        std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
// // //        std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
// // //      if(func<1.E-30)  std::cout<<"evaluate = "<<func<<" x,y,z ="<<_x<<" "<<_y<<_z<<" "<<std::endl;
// //  
// //        if(func<1.E-30) func=1.E-30; 
// // //       func = func/(xdif*ydif*zdif); 
// //        func=func*pdf;
// // //        std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
// // //        std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
// // //        std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
// // //       exit(0);
// // //      func=func/(1.+maxDegree1)/(1.+maxDegree2)/(1.+maxDegree3);
// // //      std::cout<<"eval = "<<func/(xdif*ydif*zdif)<<std::endl;
// //       return  func/(xdif*ydif*zdif);
//          return pdf;
// } 

fptype RooProdGenePdfBernsteinEffiNorm::device_coeffbinomial  (fptype enne, fptype kappa) const {
//fptype RooProdGenePdfBernsteinEffiNorm::device_coeffbinomial  (int enne, int kappa) const {
        fptype factor=1.;
        for(fptype i = 1; i <=kappa; ++i) {
          factor *= (enne+1-i)/i; 
        }	 
        if (factor<=0 ){
	 printf("Error in RooProdGenePdfBernsteinEffiNorm coeffbinomial=> factor = %f enne=%f kappa=%f",factor,enne,kappa);
         return 0;
	} 
       return factor;
}
//==============================================================================================
fptype  RooProdGenePdfBernsteinEffiNorm::gridIntegral() {

       printf("Compute Integral %d on a grid \n",_ilNorm);

       unsigned int xNGrid = 100;
       unsigned int yNGrid = 100;
       unsigned int zNGrid = 100;
       double xmin = _x.min();
       double xdif = _x.max()-xmin;
       double ymin = _y.min();
       double ydif = _y.max()-ymin;
       double zmin = _z.min();
       double zdif = _z.max()-zmin;
       double xWidth = xdif/xNGrid;
       double yWidth = ydif/yNGrid;
       double zWidth = zdif/zNGrid;
       double vol3D = xWidth*yWidth*zWidth;
       printf("vol3d = %f\n",vol3D);
       double integ=0.;
       double pdf  =0.;
       double x = 0.;
       double y = 0.;
       double z = 0.;
       for( unsigned int igi = 0; igi <xNGrid ; ++igi){
        for( unsigned int igj = 0; igj <yNGrid ; ++igj){
         for( unsigned int igk = 0; igk <zNGrid ; ++igk){
	  x = xmin + xWidth*(0.5 + igi);
	  y = ymin + yWidth*(0.5 + igj);
	  z = zmin + zWidth*(0.5 + igk);
	 
 	 pdf  = 9.0/(32.0*TMath::Pi());
 	 if	(_ilNorm==0){
 	  pdf *= 3./4.*(1.0-y)*(1.0+y);
 	 }else if(_ilNorm==1){
 	  pdf *= y*y;
 	 }else if(_ilNorm==2){
 	  pdf *= -1./4.*(1.0-y)*(1.0+y)*(2.*x*x-1.);
 	 }else if(_ilNorm==3){
 	  pdf *= -y*y*(2.*x*x-1.);
 	 }else if(_ilNorm==4){
 	  pdf *= 0.5*(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Cos(2.0*z);
 	 }else if(_ilNorm==5){
 	  pdf *= -TMath::Cos(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
 	 }else if(_ilNorm==6){
 	  pdf *= -TMath::Cos(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
 	 }else if(_ilNorm==7){
 	  pdf *= -TMath::Sin(z)*2.*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x));
 	 }else if(_ilNorm==8){
 	  pdf *= -TMath::Sin(z)*2.*x*y*sqrt((1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)) ;
 	 }else if(_ilNorm==9){
 	  pdf *= 2.0*(1.-y*y)*x;
 	 }else if(_ilNorm==10){
 	  pdf *= -(1.0-y)*(1.0+y)*(1.0-x)*(1.0+x)*TMath::Sin(2.0*z);
 	 }else{
 	  return -10000.;
 	 }
 
 	 x    =( x-xmin)/xdif;
 	 y    =( y-ymin)/ydif;
 	 z    =( z-zmin)/zdif;
 
 	 double sx[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
 
 	 double sy[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
 
 	 double sz[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
 
 
 	 sx[0]=1.0;
 	 sy[0]=1.0;
 	 sz[0]=1.0;
 	 for( int i = 1; i <= _maxDegreeV ; ++i){
 	  sx[i]= sx[i-1]*(1.-x);
 	  sy[i]= sy[i-1]*(1.-y);
 	  sz[i]= sz[i-1]*(1.-z);
 	 }
 
 	 fptype tx = 1.0;
 	 fptype ty = 1.0;
 	 fptype tz = 1.0;
 	 for( unsigned int i = 0; i <=6 ; ++i){
 	  int ii = 6-i;
 
//	   for( int i = 1; i <= max(maxDegree3,max(maxDegree1,maxDegree2)) ; ++i){
 	  if(ii<=_maxDegree1) {sx[ii]*= tx*device_coeffbinomial(_maxDegree1,ii);tx*=x;}
 	  if(ii<=_maxDegree2) {sy[ii]*= ty*device_coeffbinomial(_maxDegree2,ii);ty*=y;}
 	  if(ii<=_maxDegree3) {sz[ii]*= tz*device_coeffbinomial(_maxDegree3,ii);tz*=z;}
	
 	 }
 
//
//	    double bernknvalx = 0.;
//	    double bernknvaly = 0.;
//	    double bernknvalz = 0.;
 
 
 	 int ipar =0;
 	 double func =0.0;
 
//	   double tx = 1.;
//	   long double sx = sx_ini;
 	 for(int i = 0; i <= _maxDegree1 ; ++i) {
//	     bernknvalx =  device_coeffbinomial(_maxDegree1,i)*tx*sx[_maxDegree1-i];
 
//	     double ty = 1.;
//	     long double sy =sy_ini;
 	   for(int j = 0; j <= _maxDegree2 ; ++j) {
//	    bernknvaly =  device_coeffbinomial(_maxDegree2,j)*ty*sy[_maxDegree2-j];
	
//	      double tz = 1.;
//	      long double sz = sz_ini;
 	    for(int k = 0; k <= _maxDegree3 ; ++k) {
//	     bernknvalz =  device_coeffbinomial(_maxDegree3,k)*tz*sz[_maxDegree3-k];
//	     func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
 	     func += ((RooAbsReal&) _coefList[ipar]).getVal()*sx[_maxDegree1-i]*sy[_maxDegree2-j]*sz[_maxDegree3-k];
//	     func += ((RooAbsReal&) _coefList[ipar]).getVal()*bernknvalx*bernknvaly*bernknvalz;
//	     std::cout<<"coeff p("<<ipar<<") = "<<((RooAbsReal&) _coefList[ipar]).getVal()<<std::endl;
	     ipar++;
//	     tz *= z;
	    }
//	     ty *= y;
 	   }
//	     tx *= x;
 	 }
//	    std::cout<<_x.min()<<" "<<_x.max()<<std::endl;
//	    std::cout<<_y.min()<<" "<<_y.max()<<std::endl;
//	    std::cout<<_z.min()<<" "<<_z.max()<<std::endl;
//	  if(func<1.E-30)  std::cout<<"evaluate = "<<func<<" x,y,z ="<<_x<<" "<<_y<<_z<<" "<<std::endl;
 
	   if(func<1.E-30) func=1.E-30;
//	   func = func/(xdif*ydif*zdif);

 	  integ += func*pdf;
         }
        }
       }
//       return integ;
       return integ*vol3D;
//       return integ*vol3D/(xdif*ydif*zdif);

}
//==============================================================================================
fptype  RooProdGenePdfBernsteinEffiNorm::device_bernsteinkn_func(fptype x, fptype enne, fptype kappa) const{
   return RooProdGenePdfBernsteinEffiNorm::device_coeffbinomial(enne,kappa)*pow(x,kappa)*pow(1.0-x,enne-kappa);
}
//==============================================================================================
