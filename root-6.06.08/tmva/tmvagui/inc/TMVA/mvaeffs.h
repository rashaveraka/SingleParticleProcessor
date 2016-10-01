#ifndef mvaeffs__HH
#define mvaeffs__HH
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;

#include "tmvaglob.h"

#include "RQ_OBJECT.h"

#include "TH1.h"
#include "TROOT.h"
#include "TList.h"
#include "TIterator.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH2.h"
#include "TFormula.h"
#include "TFile.h"
#include "TApplication.h"
#include "TKey.h"
#include "TClass.h"
#include "TGaxis.h"

#include "TGWindow.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"

namespace TMVA{

   void mvaeffs( TString fin = "TMVA.root", 
                 Bool_t useTMVAStyle = kTRUE, TString formula="S/sqrt(S+B)" );

   // this macro plots the signal and background efficiencies
   // as a function of the MVA cut.


   class MethodInfo : public TNamed {
   public:
   MethodInfo() :
      methodName(""),
         methodTitle(""),
         sig(0),
         bgd(0),
         origSigE(0),
         origBgdE(0),
         sigE(0),
         bgdE(0),
         purS(0),
         sSig(0),
         effpurS(0),
         canvas(0),
         line1(0),
         line2(0),
         rightAxis(0),
         maxSignificance(0),
         maxSignificanceErr(0)
         {}
      virtual ~MethodInfo();

      TString  methodName;
      TString  methodTitle;
      TH1*     sig;
      TH1*     bgd;
      TH1*     origSigE;
      TH1*     origBgdE;
      TH1*     sigE;
      TH1*     bgdE;
      TH1*     purS;
      TH1*     sSig;    
      TH1*     effpurS;
      TCanvas* canvas;
      TLatex*  line1;
      TLatex*  line2;
      TGaxis*  rightAxis;
      Double_t maxSignificance;
      Double_t maxSignificanceErr;

      void SetResultHists(); 

      ClassDef(MethodInfo,0);
   };

  class StatDialogMVAEffs {  

      RQ_OBJECT("StatDialogMVAEffs")
      
   public:

      StatDialogMVAEffs(const TGWindow* p, Float_t ns, Float_t nb);
      virtual ~StatDialogMVAEffs();
   
      void SetFormula(const TString& f) { fFormula = f; }
      TString GetFormula();
      TString GetFormulaString() { return fFormula; }
      TString GetLatexFormula();
   
      void ReadHistograms(TFile* file);
      void UpdateSignificanceHists();
      void DrawHistograms();

      void RaiseDialog() { if (fMain) { fMain->RaiseWindow(); fMain->Layout(); fMain->MapWindow(); } }

   private:

      TGMainFrame *fMain;
      Float_t fNSignal;
      Float_t fNBackground;  
      TString fFormula;
      TList * fInfoList;

      TGNumberEntry* fSigInput;
      TGNumberEntry* fBkgInput;

      TGHorizontalFrame* fButtons;
      TGTextButton* fDrawButton;
      TGTextButton* fCloseButton;

      Int_t maxLenTitle;

      void UpdateCanvases();

   public:

      // slots
      void SetNSignal(); //*SIGNAL*
      void SetNBackground(); //*SIGNAL*
      void Redraw(); //*SIGNAL*
      void Close(); //*SIGNAL*

      // result printing
      void PrintResults( const MethodInfo* info );
   };

}
#endif
