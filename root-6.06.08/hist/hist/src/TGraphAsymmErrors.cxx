// @(#)root/hist:$Id$
// Author: Rene Brun   03/03/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <string.h>

#include "Riostream.h"
#include "TEfficiency.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TArrow.h"
#include "TBox.h"
#include "TVirtualPad.h"
#include "TF1.h"
#include "TH1.h"
#include "TVector.h"
#include "TVectorD.h"
#include "TClass.h"
#include "Math/QuantFuncMathCore.h"

ClassImp(TGraphAsymmErrors)

/** \class TGraphAsymmErrors
    \ingroup Hist
TGraph with assymetric error bars.

The TGraphAsymmErrors painting is performed thanks to the TGraphPainter
class. All details about the various painting options are given in this class.
<p>
The picture below gives an example:

Begin_Macro(source)
{
   c1 = new TCanvas("c1","A Simple Graph with assymetric error bars",200,10,700,500);
   c1->SetFillColor(42);
   c1->SetGrid();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   const Int_t n = 10;
   Double_t x[n]   = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
   Double_t y[n]   = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
   Double_t exl[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
   Double_t eyl[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
   Double_t exh[n] = {.02,.08,.05,.05,.03,.03,.04,.05,.06,.03};
   Double_t eyh[n] = {.6,.5,.4,.3,.2,.2,.3,.4,.5,.6};
   gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
   gr->SetTitle("TGraphAsymmErrors Example");
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("ALP");
   return c1;
}
End_Macro
*/


////////////////////////////////////////////////////////////////////////////////
/// TGraphAsymmErrors default constructor.

TGraphAsymmErrors::TGraphAsymmErrors(): TGraph()
{
   fEXlow       = 0;
   fEYlow       = 0;
   fEXhigh      = 0;
   fEYhigh      = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphAsymmErrors copy constructor

TGraphAsymmErrors::TGraphAsymmErrors(const TGraphAsymmErrors &gr)
       : TGraph(gr)
{
   if (!CtorAllocate()) return;
   Int_t n = fNpoints*sizeof(Double_t);
   memcpy(fEXlow, gr.fEXlow, n);
   memcpy(fEYlow, gr.fEYlow, n);
   memcpy(fEXhigh, gr.fEXhigh, n);
   memcpy(fEYhigh, gr.fEYhigh, n);
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphAsymmErrors assignment operator

TGraphAsymmErrors& TGraphAsymmErrors::operator=(const TGraphAsymmErrors &gr)
{
   if(this!=&gr) {
      TGraph::operator=(gr);
      // delete arrays
      if (fEXlow) delete [] fEXlow;
      if (fEYlow) delete [] fEYlow;
      if (fEXhigh) delete [] fEXhigh;
      if (fEYhigh) delete [] fEYhigh;

      if (!CtorAllocate()) return *this;
      Int_t n = fNpoints*sizeof(Double_t);
      memcpy(fEXlow, gr.fEXlow, n);
      memcpy(fEYlow, gr.fEYlow, n);
      memcpy(fEXhigh, gr.fEXhigh, n);
      memcpy(fEYhigh, gr.fEYhigh, n);
   }
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphAsymmErrors normal constructor.
///
/// the arrays are preset to zero

TGraphAsymmErrors::TGraphAsymmErrors(Int_t n)
       : TGraph(n)
{
   if (!CtorAllocate()) return;
   FillZero(0, fNpoints);
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphAsymmErrors normal constructor.
///
/// if exl,h or eyl,h are null, the corresponding arrays are preset to zero

TGraphAsymmErrors::TGraphAsymmErrors(Int_t n, const Float_t *x, const Float_t *y, const Float_t *exl, const Float_t *exh, const Float_t *eyl, const Float_t *eyh)
       : TGraph(n,x,y)
{
   if (!CtorAllocate()) return;

   for (Int_t i=0;i<n;i++) {
      if (exl) fEXlow[i]  = exl[i];
      else     fEXlow[i]  = 0;
      if (exh) fEXhigh[i] = exh[i];
      else     fEXhigh[i] = 0;
      if (eyl) fEYlow[i]  = eyl[i];
      else     fEYlow[i]  = 0;
      if (eyh) fEYhigh[i] = eyh[i];
      else     fEYhigh[i] = 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphAsymmErrors normal constructor.
///
/// if exl,h or eyl,h are null, the corresponding arrays are preset to zero

TGraphAsymmErrors::TGraphAsymmErrors(Int_t n, const Double_t *x, const Double_t *y, const Double_t *exl, const Double_t *exh, const Double_t *eyl, const Double_t *eyh)
       : TGraph(n,x,y)
{
   if (!CtorAllocate()) return;

   n = fNpoints*sizeof(Double_t);
   if(exl) { memcpy(fEXlow, exl, n);
   } else { memset(fEXlow, 0, n); }
   if(exh) { memcpy(fEXhigh, exh, n);
   } else { memset(fEXhigh, 0, n); }
   if(eyl) { memcpy(fEYlow, eyl, n);
   } else { memset(fEYlow, 0, n); }
   if(eyh) { memcpy(fEYhigh, eyh, n);
   } else { memset(fEYhigh, 0, n); }
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor with six vectors of floats in input
/// A grapherrors is built with the X coordinates taken from vx and Y coord from vy
/// and the errors from vectors vexl/h and veyl/h.
/// The number of points in the graph is the minimum of number of points
/// in vx and vy.

TGraphAsymmErrors::TGraphAsymmErrors(const TVectorF  &vx, const TVectorF  &vy, const TVectorF  &vexl, const TVectorF  &vexh, const TVectorF  &veyl, const TVectorF  &veyh)
                  :TGraph()
{
   fNpoints = TMath::Min(vx.GetNrows(), vy.GetNrows());
   if (!TGraph::CtorAllocate()) return;
   if (!CtorAllocate()) return;
   Int_t ivxlow  = vx.GetLwb();
   Int_t ivylow  = vy.GetLwb();
   Int_t ivexllow = vexl.GetLwb();
   Int_t ivexhlow = vexh.GetLwb();
   Int_t iveyllow = veyl.GetLwb();
   Int_t iveyhlow = veyh.GetLwb();
      for (Int_t i=0;i<fNpoints;i++) {
      fX[i]      = vx(i+ivxlow);
      fY[i]      = vy(i+ivylow);
      fEXlow[i]  = vexl(i+ivexllow);
      fEYlow[i]  = veyl(i+iveyllow);
      fEXhigh[i] = vexh(i+ivexhlow);
      fEYhigh[i] = veyh(i+iveyhlow);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor with six vectors of doubles in input
/// A grapherrors is built with the X coordinates taken from vx and Y coord from vy
/// and the errors from vectors vexl/h and veyl/h.
/// The number of points in the graph is the minimum of number of points
/// in vx and vy.

TGraphAsymmErrors::TGraphAsymmErrors(const TVectorD &vx, const TVectorD &vy, const TVectorD &vexl, const TVectorD &vexh, const TVectorD &veyl, const TVectorD &veyh)
                  :TGraph()
{
   fNpoints = TMath::Min(vx.GetNrows(), vy.GetNrows());
   if (!TGraph::CtorAllocate()) return;
   if (!CtorAllocate()) return;
   Int_t ivxlow  = vx.GetLwb();
   Int_t ivylow  = vy.GetLwb();
   Int_t ivexllow = vexl.GetLwb();
   Int_t ivexhlow = vexh.GetLwb();
   Int_t iveyllow = veyl.GetLwb();
   Int_t iveyhlow = veyh.GetLwb();
      for (Int_t i=0;i<fNpoints;i++) {
      fX[i]      = vx(i+ivxlow);
      fY[i]      = vy(i+ivylow);
      fEXlow[i]  = vexl(i+ivexllow);
      fEYlow[i]  = veyl(i+iveyllow);
      fEXhigh[i] = vexh(i+ivexhlow);
      fEYhigh[i] = veyh(i+iveyhlow);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphAsymmErrors constructor importing its parameters from the TH1 object passed as argument
/// the low and high errors are set to the bin error of the histogram.

TGraphAsymmErrors::TGraphAsymmErrors(const TH1 *h)
       : TGraph(h)
{
   if (!CtorAllocate()) return;

   for (Int_t i=0;i<fNpoints;i++) {
      fEXlow[i]  = h->GetBinWidth(i+1)*gStyle->GetErrorX();
      fEXhigh[i] = fEXlow[i];
      fEYlow[i]  = h->GetBinError(i+1);
      fEYhigh[i] = fEYlow[i];
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Creates a TGraphAsymmErrors by dividing two input TH1 histograms:
/// pass/total. (see TGraphAsymmErrors::Divide)

TGraphAsymmErrors::TGraphAsymmErrors(const TH1* pass, const TH1* total, Option_t *option)
   : TGraph((pass)?pass->GetNbinsX():0)
{
   if (!pass || !total) {
      Error("TGraphAsymmErrors","Invalid histogram pointers");
      return;
   }
   if (!CtorAllocate()) return;

   std::string sname = "divide_" + std::string(pass->GetName()) + "_by_" +
      std::string(total->GetName());
   SetName(sname.c_str());
   SetTitle(pass->GetTitle());

   //copy style from pass
   pass->TAttLine::Copy(*this);
   pass->TAttFill::Copy(*this);
   pass->TAttMarker::Copy(*this);

   Divide(pass, total, option);
}


////////////////////////////////////////////////////////////////////////////////
/// TGraphAsymmErrors default destructor.

TGraphAsymmErrors::~TGraphAsymmErrors()
{
   if(fEXlow) delete [] fEXlow;
   if(fEXhigh) delete [] fEXhigh;
   if(fEYlow) delete [] fEYlow;
   if(fEYhigh) delete [] fEYhigh;
}


////////////////////////////////////////////////////////////////////////////////
/// apply a function to all data points
/// y = f(x,y)
///
/// Errors are calculated as eyh = f(x,y+eyh)-f(x,y) and
/// eyl = f(x,y)-f(x,y-eyl)
///
/// Special treatment has to be applied for the functions where the
/// role of "up" and "down" is reversed.
/// function suggested/implemented by Miroslav Helbich <helbich@mail.desy.de>

void TGraphAsymmErrors::Apply(TF1 *f)
{
   Double_t x,y,exl,exh,eyl,eyh,eyl_new,eyh_new,fxy;

   if (fHistogram) {
      delete fHistogram;
      fHistogram = 0;
   }
   for (Int_t i=0;i<GetN();i++) {
      GetPoint(i,x,y);
      exl=GetErrorXlow(i);
      exh=GetErrorXhigh(i);
      eyl=GetErrorYlow(i);
      eyh=GetErrorYhigh(i);

      fxy = f->Eval(x,y);
      SetPoint(i,x,fxy);

      // in the case of the functions like y-> -1*y the roles of the
      // upper and lower error bars is reversed
      if (f->Eval(x,y-eyl)<f->Eval(x,y+eyh)) {
         eyl_new = TMath::Abs(fxy - f->Eval(x,y-eyl));
         eyh_new = TMath::Abs(f->Eval(x,y+eyh) - fxy);
      }
      else {
         eyh_new = TMath::Abs(fxy - f->Eval(x,y-eyl));
         eyl_new = TMath::Abs(f->Eval(x,y+eyh) - fxy);
      }

      //error on x doesn't change
      SetPointError(i,exl,exh,eyl_new,eyh_new);
   }
   if (gPad) gPad->Modified();
}

////////////////////////////////////////////////////////////////////////////////
///This function is only kept for backward compatibility.
///You should rather use the Divide method.
///It calls Divide(pass,total,"cl=0.683 b(1,1) mode") which is equivalent to the
///former BayesDivide method.

void TGraphAsymmErrors::BayesDivide(const TH1* pass, const TH1* total, Option_t *)
{
   Divide(pass,total,"cl=0.683 b(1,1) mode");
}

////////////////////////////////////////////////////////////////////////////////
/// Fill this TGraphAsymmErrors by dividing two 1-dimensional histograms pass/total
///
/// This method serves two purposes:
///
/// 1) calculating efficiencies:
/// ----------------------------
///
/// The assumption is that the entries in "pass" are a subset of those in
/// "total". That is, we create an "efficiency" graph, where each entry is
/// between 0 and 1, inclusive.
///
/// If the histograms are not filled with unit weights, the number of effective
/// entries is used to normalise the bin contents which might lead to wrong results.
/// Begin_Latex effective entries = #frac{(#sum w_{i})^{2}}{#sum w_{i}^{2}}End_Latex
///
/// The points are assigned a x value at the center of each histogram bin.
/// The y values are Begin_Latex eff = #frac{pass}{total} End_Latex for all options except for the
/// bayesian methods where the result depends on the chosen option.
///
/// If the denominator becomes 0 or pass >  total, the corresponding bin is
/// skipped.
///
/// 2) calculating ratios of two Poisson means (option 'pois'):
/// --------------------------------------------------------------
///
/// The two histograms are interpreted as independent Poisson processes and the ratio
/// Begin_Latex #tau = #frac{n_{1}}{n_{2}} = #frac{#varepsilon}{1 - #varepsilon} with #varepsilon = #frac{n_{1}}{n_{1} + n_{2}} End_Latex
/// The histogram 'pass' is interpreted as n_{1} and the total histogram
/// is used for n_{2}
///
/// The (asymmetric) uncertainties of the Poisson ratio are linked to the uncertainties
/// of efficiency by a parameter transformation:
/// Begin_Latex #Delta #tau_{low/up} = #frac{1}{(1 - #varepsilon)^{2}} #Delta #varepsilon_{low/up} End_Latex
///
/// The x errors span each histogram bin (lowedge ... lowedge+width)
/// The y errors depend on the chosen statistic methode which can be determined
/// by the options given below. For a detailed description of the used statistic
/// calculations please have a look at the corresponding functions!
///
/// Options:
/// - v     : verbose mode: prints information about the number of used bins
///           and calculated efficiencies with their errors
/// - cl=x  : determine the used confidence level (0<x<1) (default is 0.683)
/// - cp    : Clopper-Pearson interval (see TEfficiency::ClopperPearson)
/// - w     : Wilson interval (see TEfficiency::Wilson)
/// - n     : normal approximation propagation (see TEfficiency::Normal)
/// - ac    : Agresti-Coull interval (see TEfficiency::AgrestiCoull)
/// - fc    : Feldman-Cousins interval (see TEfficiency::FeldmanCousinsInterval)
/// - midp  : Lancaster mid-P interval (see TEfficiency::MidPInterval)
/// - b(a,b): bayesian interval using a prior probability ~Beta(a,b); a,b > 0
///           (see TEfficiency::Bayesian)
/// - mode  : use mode of posterior for Bayesian interval (default is mean)
/// - shortest: use shortest interval (done by default if mode is set)
/// - central: use central interval (done by default if mode is NOT set)
/// - pois: interpret histograms as poisson ratio instead of efficiency
/// - e0    : plot (in Bayesian case) efficiency and interval for bins where total=0
///           (default is to skip them)
///
/// Note:
/// Unfortunately there is no straightforward approach for determining a confidence
/// interval for a given confidence level. The actual coverage probability of the
/// confidence interval oscillates significantly according to the total number of
/// events and the true efficiency. In order to decrease the impact of this
/// oscillation on the actual coverage probability a couple of approximations and
/// methodes has been developped. For a detailed discussion, please have a look at
/// this statistical paper:
/// <a href="http://www-stat.wharton.upenn.edu/~tcai/paper/Binomial-StatSci.pdf"
/// > http://www-stat.wharton.upenn.edu/~tcai/paper/Binomial-StatSci.pdf</a>

void TGraphAsymmErrors::Divide(const TH1* pass, const TH1* total, Option_t *opt)
{
   //check pointers
   if(!pass || !total) {
      Error("Divide","one of the passed pointers is zero");
      return;
   }

   //check dimension of histograms; only 1-dimensional ones are accepted
   if((pass->GetDimension() > 1) || (total->GetDimension() > 1)) {
      Error("Divide","passed histograms are not one-dimensional");
      return;
   }

   //check whether histograms are filled with weights -> use number of effective
   //entries
   Bool_t bEffective = false;
   //compare sum of weights with sum of squares of weights
   Double_t stats[TH1::kNstat];
   pass->GetStats(stats);
   if (TMath::Abs(stats[0] -stats[1]) > 1e-6)
      bEffective = true;
   total->GetStats(stats);
   if (TMath::Abs(stats[0] -stats[1]) > 1e-6)
      bEffective = true;

   // we do not want to ignore the weights
   // if (bEffective && (pass->GetSumw2()->fN == 0 || total->GetSumw2()->fN == 0) ) {
   //    Warning("Divide","histogram have been computed with weights but the sum of weight squares are not stored in the histogram. Error calculation is performed ignoring the weights");
   //    bEffective = false;
   // }

   //parse option
   TString option = opt;
   option.ToLower();

   Bool_t bVerbose = false;
   //pointer to function returning the boundaries of the confidence interval
   //(is only used in the frequentist cases.)
   Double_t (*pBound)(Double_t,Double_t,Double_t,Bool_t) = &TEfficiency::ClopperPearson; // default method
   //confidence level
   Double_t conf = 0.682689492137;
   //values for bayesian statistics
   Bool_t bIsBayesian = false;
   Double_t alpha = 1;
   Double_t beta = 1;

   //verbose mode
   if(option.Contains("v")) {
      option.ReplaceAll("v","");
      bVerbose = true;
   }

   //confidence level
   if(option.Contains("cl=")) {
      Double_t level = -1;
      // coverity [secure_coding : FALSE]
      sscanf(strstr(option.Data(),"cl="),"cl=%lf",&level);
      if((level > 0) && (level < 1))
         conf = level;
      else
         Warning("Divide","given confidence level %.3lf is invalid",level);
      option.ReplaceAll("cl=","");
   }

   //normal approximation
   if(option.Contains("n")) {
      option.ReplaceAll("n","");
      pBound = &TEfficiency::Normal;
   }

   //clopper pearson interval
   if(option.Contains("cp")) {
      option.ReplaceAll("cp","");
      pBound = &TEfficiency::ClopperPearson;
   }

   //wilson interval
   if(option.Contains("w")) {
      option.ReplaceAll("w","");
      pBound = &TEfficiency::Wilson;
   }

   //agresti coull interval
   if(option.Contains("ac")) {
      option.ReplaceAll("ac","");
      pBound = &TEfficiency::AgrestiCoull;
   }
   // Feldman-Cousins interval
   if(option.Contains("fc")) {
      option.ReplaceAll("fc","");
      pBound = &TEfficiency::FeldmanCousins;
   }
   // mid-P Lancaster interval
   if(option.Contains("midp")) {
      option.ReplaceAll("midp","");
      pBound = &TEfficiency::MidPInterval;
   }
   
   //bayesian with prior
   if(option.Contains("b(")) {
      Double_t a = 0;
      Double_t b = 0;
      sscanf(strstr(option.Data(),"b("),"b(%lf,%lf)",&a,&b);
      if(a > 0)
         alpha = a;
      else
         Warning("Divide","given shape parameter for alpha %.2lf is invalid",a);
      if(b > 0)
         beta = b;
      else
         Warning("Divide","given shape parameter for beta %.2lf is invalid",b);
      option.ReplaceAll("b(","");
      bIsBayesian = true;
   }

   // use posterior mode
   Bool_t usePosteriorMode = false;
   if(bIsBayesian && option.Contains("mode") ) {
      usePosteriorMode = true;
      option.ReplaceAll("mode","");
   }

   Bool_t plot0Bins = false;
   if(option.Contains("e0") ) {
      plot0Bins = true;
      option.ReplaceAll("e0","");
   }

   Bool_t useShortestInterval = false;
   if (bIsBayesian && ( option.Contains("sh") || (usePosteriorMode && !option.Contains("cen") ) ) ) {
      useShortestInterval = true;
   }

   // interpret as Poisson ratio
   Bool_t bPoissonRatio = false;
   if(option.Contains("pois") ) {
      bPoissonRatio = true;
      option.ReplaceAll("pois","");
   }

   // weights works only in case of Normal approximation or Bayesian for binomial interval
   // in case of Poisson ratio we can use weights by rescaling the obtained results using the effective entries
   if ( ( bEffective && !bPoissonRatio) && !bIsBayesian && pBound != &TEfficiency::Normal ) {
      Warning("Divide","Histograms have weights: only Normal or Bayesian error calculation is supported");
      Info("Divide","Using now the Normal approximation for weighted histograms");
   }

   if(bPoissonRatio)
   {
     if(pass->GetDimension() != total->GetDimension()) {
       Error("Divide","passed histograms are not of the same dimension");
       return;
     }

     if(!TEfficiency::CheckBinning(*pass,*total)) {
       Error("Divide","passed histograms are not consistent");
       return;
     }
   }
   else
   {
     //check consistency of histograms, allowing weights
     if(!TEfficiency::CheckConsistency(*pass,*total,"w")) {
       Error("Divide","passed histograms are not consistent");
       return;
     }
   }

   //Set the graph to have a number of points equal to the number of histogram
   //bins
   Int_t nbins = pass->GetNbinsX();
   Set(nbins);

   // Ok, now set the points for each bin
   // (Note: the TH1 bin content is shifted to the right by one:
   //  bin=0 is underflow, bin=nbins+1 is overflow.)

   //efficiency with lower and upper boundary of confidence interval
   double eff, low, upper;
   //this keeps track of the number of points added to the graph
   int npoint=0;
   //number of total and passed events
   Double_t t = 0 , p = 0;
   Double_t tw = 0, tw2 = 0, pw = 0, pw2 = 0, wratio = 1; // for the case of weights
   //loop over all bins and fill the graph
   for (Int_t b=1; b<=nbins; ++b) {

      // default value when total =0;
      eff = 0;
      low = 0;
      upper = 0;

      // special case in case of weights we have to consider the sum of weights and the sum of weight squares
       if(bEffective) {
          tw =  total->GetBinContent(b);
          tw2 = (total->GetSumw2()->fN > 0) ? total->GetSumw2()->At(b) : tw;
          pw =  pass->GetBinContent(b);
          pw2 = (pass->GetSumw2()->fN > 0) ? pass->GetSumw2()->At(b) : pw;

          if(bPoissonRatio)
          {
             // tw += pw;
             // tw2 += pw2;
             // compute effective entries
             // special case is (pw=0, pw2=0) in this case is like unweighted
             if (pw == 0 && pw2 == 0)
                p = pw;
             else
                p = (pw*pw)/pw2;

             if (tw == 0 && tw2 == 0)
                t = tw;
             else
                t = (tw*tw)/tw2;

             if (p > 0 && tw > 0)
                wratio = (pw*t)/(p * tw);
             else if (p == 0 && tw > 0)
                wratio = t/tw;
             else if (p > 0)
                wratio = pw/p;
             else {
                // case both pw and tw are zero - we skip these bins
                if (!plot0Bins) continue; // skip bins with total <= 0
             }

             t += p;
          }
          else
             if (tw <= 0 && !plot0Bins) continue; // skip bins with total <= 0

          // in the case of weights have the formula only for
          // the normal and  bayesian statistics (see below)

       }

       //use bin contents
       else {
          t = int( total->GetBinContent(b) + 0.5);
          p = int(pass->GetBinContent(b) + 0.5);

          if(bPoissonRatio) t += p;

          if (!t && !plot0Bins) continue; // skip bins with total = 0
       }


      //using bayesian statistics
      if(bIsBayesian) {
         double aa,bb;

         if ((bEffective && !bPoissonRatio) && tw2 <= 0) {
            // case of bins with zero errors
            eff = pw/tw;
            low = eff; upper = eff;
         }
         else {

            if (bEffective && !bPoissonRatio) {
               // tw/tw2 renormalize the weights
               double norm = tw/tw2;  // case of tw2 = 0 is treated above
               aa =  pw * norm + alpha;
               bb =  (tw - pw) * norm + beta;
            }
            else {
               aa = double(p) + alpha;
               bb = double(t-p) + beta;
            }
            if (usePosteriorMode)
               eff = TEfficiency::BetaMode(aa,bb);
            else
               eff = TEfficiency::BetaMean(aa,bb);

            if (useShortestInterval) {
               TEfficiency::BetaShortestInterval(conf,aa,bb,low,upper);
            }
            else {
               low = TEfficiency::BetaCentralInterval(conf,aa,bb,false);
               upper = TEfficiency::BetaCentralInterval(conf,aa,bb,true);
            }
         }
      }
      // case of non-bayesian statistics
      else {
         if (bEffective && !bPoissonRatio) {

            if (tw > 0) {

               eff = pw/tw;

               // use normal error calculation using variance of MLE with weights (F.James 8.5.2)
               // this is the same formula used in ROOT for TH1::Divide("B")

               double variance = ( pw2 * (1. - 2 * eff) + tw2 * eff *eff ) / ( tw * tw) ;
               double sigma = sqrt(variance);

               double prob = 0.5 * (1.-conf);
               double delta = ROOT::Math::normal_quantile_c(prob, sigma);
               low = eff - delta;
               upper = eff + delta;
               if (low < 0) low = 0;
               if (upper > 1) upper = 1.;
            }
         }
         else {
            // when not using weights (all cases) or in case of  Poisson ratio with weights
            if(t)
               eff = ((Double_t)p)/t;

            low = pBound(t,p,conf,false);
            upper = pBound(t,p,conf,true);
         }
      }
      // treat as Poisson ratio
      if(bPoissonRatio)
      {
        Double_t ratio = eff/(1 - eff);
        // take the intervals in eff as intervals in the Poisson ratio
        low = low/(1. - low);
        upper = upper/(1.-upper);
        eff = ratio;
        if (bEffective) {
           //scale result by the ratio of the weight
           eff *= wratio;
           low *= wratio;
           upper *= wratio; 
        }
      }
      //Set the point center and its errors
      SetPoint(npoint,pass->GetBinCenter(b),eff);
      SetPointError(npoint,
      pass->GetBinCenter(b)-pass->GetBinLowEdge(b),
      pass->GetBinLowEdge(b)-pass->GetBinCenter(b)+pass->GetBinWidth(b),
      eff-low,upper-eff);
      npoint++;//we have added a point to the graph
   }

   Set(npoint);//tell the graph how many points we've really added
   if (npoint < nbins)
      Warning("Divide","Number of graph points is different than histogram bins - %d points have been skipped",nbins-npoint);
   

   if (bVerbose) {
      Info("Divide","made a graph with %d points from %d bins",npoint,nbins);
      Info("Divide","used confidence level: %.2lf\n",conf);
      if(bIsBayesian)
         Info("Divide","used prior probability ~ beta(%.2lf,%.2lf)",alpha,beta);
      Print();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute Range

void TGraphAsymmErrors::ComputeRange(Double_t &xmin, Double_t &ymin, Double_t &xmax, Double_t &ymax) const
{
   TGraph::ComputeRange(xmin,ymin,xmax,ymax);

   for (Int_t i=0;i<fNpoints;i++) {
      if (fX[i] -fEXlow[i] < xmin) {
         if (gPad && gPad->GetLogx()) {
            if (fEXlow[i] < fX[i]) xmin = fX[i]-fEXlow[i];
            else                   xmin = TMath::Min(xmin,fX[i]/3);
         } else {
            xmin = fX[i]-fEXlow[i];
         }
      }
      if (fX[i] +fEXhigh[i] > xmax) xmax = fX[i]+fEXhigh[i];
      if (fY[i] -fEYlow[i] < ymin) {
         if (gPad && gPad->GetLogy()) {
            if (fEYlow[i] < fY[i]) ymin = fY[i]-fEYlow[i];
            else                   ymin = TMath::Min(ymin,fY[i]/3);
         } else {
            ymin = fY[i]-fEYlow[i];
         }
      }
      if (fY[i] +fEYhigh[i] > ymax) ymax = fY[i]+fEYhigh[i];
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Copy and release.

void TGraphAsymmErrors::CopyAndRelease(Double_t **newarrays,
                                       Int_t ibegin, Int_t iend, Int_t obegin)
{
   CopyPoints(newarrays, ibegin, iend, obegin);
   if (newarrays) {
      delete[] fEXlow;
      fEXlow = newarrays[0];
      delete[] fEXhigh;
      fEXhigh = newarrays[1];
      delete[] fEYlow;
      fEYlow = newarrays[2];
      delete[] fEYhigh;
      fEYhigh = newarrays[3];
      delete[] fX;
      fX = newarrays[4];
      delete[] fY;
      fY = newarrays[5];
      delete[] newarrays;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Copy errors from fE*** to arrays[***]
/// or to f*** Copy points.

Bool_t TGraphAsymmErrors::CopyPoints(Double_t **arrays,
                                     Int_t ibegin, Int_t iend, Int_t obegin)
{
   if (TGraph::CopyPoints(arrays ? arrays+4 : 0, ibegin, iend, obegin)) {
      Int_t n = (iend - ibegin)*sizeof(Double_t);
      if (arrays) {
         memmove(&arrays[0][obegin], &fEXlow[ibegin], n);
         memmove(&arrays[1][obegin], &fEXhigh[ibegin], n);
         memmove(&arrays[2][obegin], &fEYlow[ibegin], n);
         memmove(&arrays[3][obegin], &fEYhigh[ibegin], n);
      } else {
         memmove(&fEXlow[obegin], &fEXlow[ibegin], n);
         memmove(&fEXhigh[obegin], &fEXhigh[ibegin], n);
         memmove(&fEYlow[obegin], &fEYlow[ibegin], n);
         memmove(&fEYhigh[obegin], &fEYhigh[ibegin], n);
      }
      return kTRUE;
   } else {
      return kFALSE;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Should be called from ctors after fNpoints has been set
/// Note: This function should be called only from the constructor
/// since it does not delete previously existing arrays

Bool_t TGraphAsymmErrors::CtorAllocate(void)
{
   if (!fNpoints) {
      fEXlow = fEYlow = fEXhigh = fEYhigh = 0;
      return kFALSE;
   }
   fEXlow = new Double_t[fMaxSize];
   fEYlow = new Double_t[fMaxSize];
   fEXhigh = new Double_t[fMaxSize];
   fEYhigh = new Double_t[fMaxSize];
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///  protected function to perform the merge operation of a graph with asymmetric errors

Bool_t TGraphAsymmErrors::DoMerge(const TGraph *g)
{
   if (g->GetN() == 0) return kFALSE;

   Double_t * exl = g->GetEXlow();
   Double_t * exh = g->GetEXhigh();
   Double_t * eyl = g->GetEYlow();
   Double_t * eyh = g->GetEYhigh();
   if (exl == 0 || exh == 0 || eyl == 0 || eyh == 0) {
      if (g->IsA() != TGraph::Class() )
         Warning("DoMerge","Merging a %s is not compatible with a TGraphAsymmErrors - errors will be ignored",g->IsA()->GetName());
      return TGraph::DoMerge(g);
   }
   for (Int_t i = 0 ; i < g->GetN(); i++) {
      Int_t ipoint = GetN();
      Double_t x = g->GetX()[i];
      Double_t y = g->GetY()[i];
      SetPoint(ipoint, x, y);
      SetPointError(ipoint, exl[i], exh[i], eyl[i], eyh[i] );
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set zero values for point arrays in the range [begin, end)

void TGraphAsymmErrors::FillZero(Int_t begin, Int_t end,
                                 Bool_t from_ctor)
{
   if (!from_ctor) {
      TGraph::FillZero(begin, end, from_ctor);
   }
   Int_t n = (end - begin)*sizeof(Double_t);
   memset(fEXlow + begin, 0, n);
   memset(fEXhigh + begin, 0, n);
   memset(fEYlow + begin, 0, n);
   memset(fEYhigh + begin, 0, n);
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along X at point i.

Double_t TGraphAsymmErrors::GetErrorX(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (!fEXlow && !fEXhigh) return -1;
   Double_t elow=0, ehigh=0;
   if (fEXlow)  elow  = fEXlow[i];
   if (fEXhigh) ehigh = fEXhigh[i];
   return TMath::Sqrt(0.5*(elow*elow + ehigh*ehigh));
}


////////////////////////////////////////////////////////////////////////////////
/// This function is called by GraphFitChisquare.
/// It returns the error along Y at point i.

Double_t TGraphAsymmErrors::GetErrorY(Int_t i) const
{
   if (i < 0 || i >= fNpoints) return -1;
   if (!fEYlow && !fEYhigh) return -1;
   Double_t elow=0, ehigh=0;
   if (fEYlow)  elow  = fEYlow[i];
   if (fEYhigh) ehigh = fEYhigh[i];
   return TMath::Sqrt(0.5*(elow*elow + ehigh*ehigh));
}


////////////////////////////////////////////////////////////////////////////////
/// Get high error on X.

Double_t TGraphAsymmErrors::GetErrorXhigh(Int_t i) const
{
   if (i<0 || i>fNpoints) return -1;
   if (fEXhigh) return fEXhigh[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Get low error on X.

Double_t TGraphAsymmErrors::GetErrorXlow(Int_t i) const
{
   if (i<0 || i>fNpoints) return -1;
   if (fEXlow) return fEXlow[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Get high error on Y.

Double_t TGraphAsymmErrors::GetErrorYhigh(Int_t i) const
{
   if (i<0 || i>fNpoints) return -1;
   if (fEYhigh) return fEYhigh[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Get low error on Y.

Double_t TGraphAsymmErrors::GetErrorYlow(Int_t i) const
{
   if (i<0 || i>fNpoints) return -1;
   if (fEYlow) return fEYlow[i];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Adds all graphs with asymmetric errors from the collection to this graph.
/// Returns the total number of poins in the result or -1 in case of an error.

Int_t TGraphAsymmErrors::Merge(TCollection* li)
{
   TIter next(li);
   while (TObject* o = next()) {
      TGraph *g = dynamic_cast<TGraph*>(o);
      if (!g) {
         Error("Merge",
               "Cannot merge - an object which doesn't inherit from TGraph found in the list");
         return -1;
      }
      int n0 = GetN();
      int n1 = n0+g->GetN();
      Set(n1);
      Double_t * x = g->GetX();
      Double_t * y = g->GetY();
      Double_t * exlow  = g->GetEXlow();
      Double_t * exhigh = g->GetEXhigh();
      Double_t * eylow  = g->GetEYlow();
      Double_t * eyhigh = g->GetEYhigh();
      for (Int_t i = 0 ; i < g->GetN(); i++) {
         SetPoint(n0+i, x[i], y[i]);
         if (exlow)  fEXlow[n0+i]  = exlow[i];
         if (exhigh) fEXhigh[n0+i] = exhigh[i];
         if (eylow)  fEYlow[n0+i]  = eylow[i];
         if (eyhigh) fEYhigh[n0+i] = eyhigh[i];
      }
   }
   return GetN();
}

////////////////////////////////////////////////////////////////////////////////
/// Print graph and errors values.

void TGraphAsymmErrors::Print(Option_t *) const
{
   for (Int_t i=0;i<fNpoints;i++) {
      printf("x[%d]=%g, y[%d]=%g, exl[%d]=%g, exh[%d]=%g, eyl[%d]=%g, eyh[%d]=%g\n"
         ,i,fX[i],i,fY[i],i,fEXlow[i],i,fEXhigh[i],i,fEYlow[i],i,fEYhigh[i]);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Save primitive as a C++ statement(s) on output stream out

void TGraphAsymmErrors::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
   char quote = '"';
   out << "   " << std::endl;
   static Int_t frameNumber = 3000;
   frameNumber++;

   Int_t i;
   TString fXName   = TString(GetName()) + Form("_fx%d",frameNumber);
   TString fYName   = TString(GetName()) + Form("_fy%d",frameNumber);
   TString fElXName = TString(GetName()) + Form("_felx%d",frameNumber);
   TString fElYName = TString(GetName()) + Form("_fely%d",frameNumber);
   TString fEhXName = TString(GetName()) + Form("_fehx%d",frameNumber);
   TString fEhYName = TString(GetName()) + Form("_fehy%d",frameNumber);
   out << "   Double_t " << fXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fX[i] << "," << std::endl;
   out << "   " << fX[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fY[i] << "," << std::endl;
   out << "   " << fY[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fElXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEXlow[i] << "," << std::endl;
   out << "   " << fEXlow[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fElYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEYlow[i] << "," << std::endl;
   out << "   " << fEYlow[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEhXName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEXhigh[i] << "," << std::endl;
   out << "   " << fEXhigh[fNpoints-1] << "};" << std::endl;
   out << "   Double_t " << fEhYName << "[" << fNpoints << "] = {" << std::endl;
   for (i = 0; i < fNpoints-1; i++) out << "   " << fEYhigh[i] << "," << std::endl;
   out << "   " << fEYhigh[fNpoints-1] << "};" << std::endl;

   if (gROOT->ClassSaved(TGraphAsymmErrors::Class())) out<<"   ";
   else out << "   TGraphAsymmErrors *";
   out << "grae = new TGraphAsymmErrors("<< fNpoints << ","
                                    << fXName   << ","  << fYName  << ","
                                    << fElXName  << ","  << fEhXName << ","
                                    << fElYName  << ","  << fEhYName << ");"
                                    << std::endl;

   out << "   grae->SetName(" << quote << GetName() << quote << ");" << std::endl;
   out << "   grae->SetTitle(" << quote << GetTitle() << quote << ");" << std::endl;

   SaveFillAttributes(out, "grae", 0, 1001);
   SaveLineAttributes(out, "grae", 1, 1, 1);
   SaveMarkerAttributes(out, "grae", 1, 1, 1);

   if (fHistogram) {
      TString hname = fHistogram->GetName();
      hname += frameNumber;
      fHistogram->SetName(Form("Graph_%s",hname.Data()));
      fHistogram->SavePrimitive(out,"nodraw");
      out<<"   grae->SetHistogram("<<fHistogram->GetName()<<");"<<std::endl;
      out<<"   "<<std::endl;
   }

   // save list of functions
   TIter next(fFunctions);
   TObject *obj;
   while ((obj = next())) {
      obj->SavePrimitive(out, Form("nodraw #%d\n",++frameNumber));
      if (obj->InheritsFrom("TPaveStats")) {
         out << "   grae->GetListOfFunctions()->Add(ptstats);" << std::endl;
         out << "   ptstats->SetParent(grae->GetListOfFunctions());" << std::endl;
      } else {
         out << "   grae->GetListOfFunctions()->Add("
             << Form("%s%d",obj->GetName(),frameNumber) << ");" << std::endl;
      }
   }

   const char *l = strstr(option,"multigraph");
   if (l) {
      out<<"   multigraph->Add(grae,"<<quote<<l+10<<quote<<");"<<std::endl;
   } else {
      out<<"   grae->Draw("<<quote<<option<<quote<<");"<<std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set ex and ey values for point pointed by the mouse.

void TGraphAsymmErrors::SetPointError(Double_t exl, Double_t exh, Double_t eyl, Double_t eyh)
{
   Int_t px = gPad->GetEventX();
   Int_t py = gPad->GetEventY();

   //localize point to be deleted
   Int_t ipoint = -2;
   Int_t i;
   // start with a small window (in case the mouse is very close to one point)
   for (i=0;i<fNpoints;i++) {
      Int_t dpx = px - gPad->XtoAbsPixel(gPad->XtoPad(fX[i]));
      Int_t dpy = py - gPad->YtoAbsPixel(gPad->YtoPad(fY[i]));
      if (dpx*dpx+dpy*dpy < 25) {ipoint = i; break;}
   }
   if (ipoint == -2) return;

   fEXlow[ipoint]  = exl;
   fEYlow[ipoint]  = eyl;
   fEXhigh[ipoint] = exh;
   fEYhigh[ipoint] = eyh;
   gPad->Modified();
}


////////////////////////////////////////////////////////////////////////////////
/// Set ex and ey values for point number i.

void TGraphAsymmErrors::SetPointError(Int_t i, Double_t exl, Double_t exh, Double_t eyl, Double_t eyh)
{
   if (i < 0) return;
   if (i >= fNpoints) {
   // re-allocate the object
      TGraphAsymmErrors::SetPoint(i,0,0);
   }
   fEXlow[i]  = exl;
   fEYlow[i]  = eyl;
   fEXhigh[i] = exh;
   fEYhigh[i] = eyh;
}


////////////////////////////////////////////////////////////////////////////////
/// Set EXlow for point i

void TGraphAsymmErrors::SetPointEXlow(Int_t i, Double_t exl)
{
   if (i < 0) return;
   if (i >= fNpoints) {
   // re-allocate the object
      TGraphAsymmErrors::SetPoint(i,0,0);
   }
   fEXlow[i]  = exl;
}


////////////////////////////////////////////////////////////////////////////////
/// Set EXhigh for point i

void TGraphAsymmErrors::SetPointEXhigh(Int_t i, Double_t exh)
{
   if (i < 0) return;
   if (i >= fNpoints) {
   // re-allocate the object
      TGraphAsymmErrors::SetPoint(i,0,0);
   }
   fEXhigh[i]  = exh;
}


////////////////////////////////////////////////////////////////////////////////
/// Set EYlow for point i

void TGraphAsymmErrors::SetPointEYlow(Int_t i, Double_t eyl)
{
   if (i < 0) return;
   if (i >= fNpoints) {
   // re-allocate the object
      TGraphAsymmErrors::SetPoint(i,0,0);
   }
   fEYlow[i]  = eyl;
}


////////////////////////////////////////////////////////////////////////////////
/// Set EYhigh for point i

void TGraphAsymmErrors::SetPointEYhigh(Int_t i, Double_t eyh)
{
   if (i < 0) return;
   if (i >= fNpoints) {
   // re-allocate the object
      TGraphAsymmErrors::SetPoint(i,0,0);
   }
   fEYhigh[i]  = eyh;
}


////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TGraphAsymmErrors.

void TGraphAsymmErrors::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 2) {
         b.ReadClassBuffer(TGraphAsymmErrors::Class(), this, R__v, R__s, R__c);
         return;
      }
      //====process old versions before automatic schema evolution
      TGraph::Streamer(b);
      fEXlow  = new Double_t[fNpoints];
      fEYlow  = new Double_t[fNpoints];
      fEXhigh = new Double_t[fNpoints];
      fEYhigh = new Double_t[fNpoints];
      if (R__v < 2) {
         Float_t *exlow  = new Float_t[fNpoints];
         Float_t *eylow  = new Float_t[fNpoints];
         Float_t *exhigh = new Float_t[fNpoints];
         Float_t *eyhigh = new Float_t[fNpoints];
         b.ReadFastArray(exlow,fNpoints);
         b.ReadFastArray(eylow,fNpoints);
         b.ReadFastArray(exhigh,fNpoints);
         b.ReadFastArray(eyhigh,fNpoints);
         for (Int_t i=0;i<fNpoints;i++) {
            fEXlow[i]  = exlow[i];
            fEYlow[i]  = eylow[i];
            fEXhigh[i] = exhigh[i];
            fEYhigh[i] = eyhigh[i];
         }
         delete [] eylow;
         delete [] exlow;
         delete [] eyhigh;
         delete [] exhigh;
      } else {
         b.ReadFastArray(fEXlow,fNpoints);
         b.ReadFastArray(fEYlow,fNpoints);
         b.ReadFastArray(fEXhigh,fNpoints);
         b.ReadFastArray(fEYhigh,fNpoints);
      }
      b.CheckByteCount(R__s, R__c, TGraphAsymmErrors::IsA());
      //====end of old versions

   } else {
      b.WriteClassBuffer(TGraphAsymmErrors::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Swap points.

void TGraphAsymmErrors::SwapPoints(Int_t pos1, Int_t pos2)
{
   SwapValues(fEXlow,  pos1, pos2);
   SwapValues(fEXhigh, pos1, pos2);
   SwapValues(fEYlow,  pos1, pos2);
   SwapValues(fEYhigh, pos1, pos2);
   TGraph::SwapPoints(pos1, pos2);
}
