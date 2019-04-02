/***************************************************************************** 
 * Project: RooGenePdf                                                       * 
 *                                                                           * 
 * P.Dini fecit, Anno Domini MMXVIII                                         *
 * "Amicus Plato, sed magis amica veritas."                                  * 
 *                                                                           * 
 * Class to describe 3D angular PDF                                          * 
 *                                                                           * 
 *****************************************************************************/ 
#include "Riostream.h" 
#include "RooGenePdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h"
#include "RooRealProxy.h" 
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include <math.h> 
#include "TMath.h"  

ClassImp(RooGenePdf); 
   
  
   RooGenePdf::RooGenePdf(const char *name, const char *title, 
                        RooAbsReal& x, RooAbsReal& y, RooAbsReal& z, 
			const RooArgList& coefList
			) :
   RooAbsPdf(name,title), 
   _x("_x","_x",this,x),
   _y("_y","_y",this,y),
   _z("_z","_z",this,z), 
   _coefList("_coefList","_coefList",this)
{ 
    if(coefList.getSize()!=8){
    	 std::cout << "RooGenePdf ERROR: number of coefficients should be 6" << std::endl ;
    	 assert(0) ;
    }
    TIterator* coefIter = coefList.createIterator() ;
     RooAbsArg* coef ;
     while((coef = (RooAbsArg*)coefIter->Next())) {
       if (!dynamic_cast<RooAbsReal*>(coef)) {
    	 std::cout << "RooGenePdf::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
    	      << " is not of type RooAbsReal" << std::endl ;
    	 assert(0) ;
       }
       _coefList.add(*coef) ;
//   	 std::cout << "RooGenePdf::ctor(" << GetName() << ") add coefficient " << coef->GetName()<<std::endl ;
      }
//      assert(0);
     delete coefIter;
 } 

RooGenePdf::RooGenePdf(const RooGenePdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   _x("_x",this,other._x),
   _y("_y",this,other._y),
   _z("_z",this,other._z), 
   _coefList("_coefList",this,other._coefList)
{ 
} 
//
Double_t RooGenePdf::evaluate() const 
{ 
  
  double x= _x;
  double y= _y;
  double z= _z;
//   double F_S   = ((RooAbsPdf&) _coefList[0]).getVal();
//   double A_S   = ((RooAbsPdf&) _coefList[1]).getVal();
//   double A5S   = ((RooAbsPdf&) _coefList[2]).getVal();
  double F_L   = ((RooAbsPdf&) _coefList[0]).getVal();
  double P_1   = ((RooAbsPdf&) _coefList[1]).getVal();
  double P_2   = ((RooAbsPdf&) _coefList[2]).getVal();
  double P_3   = ((RooAbsPdf&) _coefList[3]).getVal();
  double P4p   = ((RooAbsPdf&) _coefList[4]).getVal();
  double P5p   = ((RooAbsPdf&) _coefList[5]).getVal();
  double P6p   = ((RooAbsPdf&) _coefList[6]).getVal();
  double P8p   = ((RooAbsPdf&) _coefList[7]).getVal();
  
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
  

//
  double pdf = 9./(32.*TMath::Pi())*(3./4.*F_T*sinyq+F_L*cosyq+\
  (1./4.*F_T*sinyq-F_L*cosyq)*cos2x+0.5*P_1*F_T*sinyq*sinxq*TMath::Cos(2.*z)+\
  sqrt(F_L*F_T)*TMath::Cos(z)*(P4p*sin2y*sin2x+P5p*sin2y*sinx)+\
  sqrt(F_L*F_T)*TMath::Sin(z)*(P6p*sin2y*sinx+P8p*sin2y*sin2x)+\
  2.*P_2*F_T*sinyq*x-P_3*F_T*sinyq*sinxq*TMath::Sin(2.*z));
//

  return pdf;
  
} 



 
