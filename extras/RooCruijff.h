//--*-C++-*--
/*
Author: Joshua Thompson
Very similar to RooBifurCB, but not the same.

Based on HN post by Steve Sekula:
http://babar-hn.slac.stanford.edu:5090/HyperNews/get/rev-SemiLep-05-01/4/3/1/1.html
                                      (x - m0)^2
The form is basically f(x) = exp( - -------------------------- )
                                    2*sigma^2 + alpha(x-m0)^2

where sigma = sigmaL and alpha = alphaL when x<m0, and 
sigma = sigmaR and alpha = alphaR when x>m0, and where m0 is the
mean of the "gaussian" part of this function. So in total, excluding
the normalization, this function has 5 parameters.
*/
#ifndef ROO_CRUIJFF
#define ROO_CRUIJFF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooCruijff : public RooAbsPdf {
public:
  RooCruijff();

  RooCruijff(const char *name, const char *title, RooAbsReal& _m,
	     RooAbsReal& _m0, RooAbsReal& _sigmaL,  RooAbsReal& _sigmaR,
	     RooAbsReal& _alphaL, RooAbsReal& _alphaR);

  RooCruijff(const RooCruijff& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooCruijff(*this,newname); }

  inline virtual ~RooCruijff() { }

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigmaL;
  RooRealProxy sigmaR;
  RooRealProxy alphaL;
  RooRealProxy alphaR;

  Double_t evaluate() const;

private:

  ClassDef(RooCruijff, 1) // Cruijff function PDF
};

#endif