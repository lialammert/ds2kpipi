
// -- CLASS DESCRIPTION [PDF] --
// Cruijff function PDF


#include "RooFit.h"
#include <iostream>
#include <math.h>

#include "RooCruijff.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

//ClassImp(RooCruijff)

RooCruijff::RooCruijff() : RooAbsPdf() {} // Default constructor implementation

RooCruijff::RooCruijff(const char *name, const char *title,
		       RooAbsReal& _m, RooAbsReal& _m0, 
		       RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
		       RooAbsReal& _alphaL, RooAbsReal& _alphaR) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigmaL("sigmaL", "SigmaL", this, _sigmaL),
  sigmaR("sigmaR", "SigmaR", this, _sigmaR),
  alphaL("alphaL", "AlphaL", this, _alphaL),
  alphaR("alphaR", "AlphaR", this, _alphaR)
{
}

RooCruijff::RooCruijff(const RooCruijff& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigmaL("sigmaL", this, other.sigmaL),  sigmaR("sigmaR", this, other.sigmaR), 
  alphaL("alphaL", this, other.alphaL),alphaR("alphaR", this, other.alphaR)
{
}

Double_t RooCruijff::evaluate() const 
{
  // build the functional form
  double sigma = 0.0;
  double alpha = 0.0;
  double dx = (m - m0);
  if(dx<0){
    sigma = sigmaL;
    alpha = alphaL;
  } else {
    sigma = sigmaR;
    alpha = alphaR;
  }
  double f = 2*sigma*sigma + alpha*dx*dx ;
  return exp(-dx*dx/f) ;
}