// @(#)root/roostats:$Id$
// Authors: Kevin Belasco        17/06/2009
// Authors: Kyle Cranmer         17/06/2009
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////////////////



#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#ifndef ROOSTATS_PdfProposal
#include "RooStats/PdfProposal.h"
#endif
#ifndef RooStats_RooStatsUtils
#include "RooStats/RooStatsUtils.h"
#endif
#ifndef ROO_ARG_SET
#include "RooArgSet.h"
#endif
#ifndef ROO_DATA_SET
#include "RooDataSet.h"
#endif
#ifndef ROO_ABS_PDF
#include "RooAbsPdf.h"
#endif
#ifndef ROO_MSG_SERVICE
#include "RooMsgService.h"
#endif
#ifndef ROO_REAL_VAR
#include "RooRealVar.h"
#endif
#ifndef ROOT_TIterator
#include "TIterator.h"
#endif

#include <map>

ClassImp(RooStats::PdfProposal);

using namespace RooFit;
using namespace RooStats;
using namespace std;

// By default, PdfProposal does NOT own the PDF that serves as the
// proposal density function
PdfProposal::PdfProposal() : ProposalFunction()
{
   fPdf = NULL;
   fOwnsPdf = kFALSE;
   fCacheSize = 1;
   fCachePosition = 0;
   fCache = NULL;
}

// By default, PdfProposal does NOT own the PDF that serves as the
// proposal density function
PdfProposal::PdfProposal(RooAbsPdf& pdf) : ProposalFunction()
{
   fPdf = &pdf;
   fOwnsPdf = kFALSE;
   fCacheSize = 1;
   fCachePosition = 0;
   fCache = NULL;
}

Bool_t PdfProposal::Equals(RooArgSet& x1, RooArgSet& x2)
{
   if (x1.equals(x2)) {
      TIterator* it = x1.createIterator();
      RooRealVar* r;
      while ((r = (RooRealVar*)it->Next()) != NULL)
         if (r->getVal() != x2.getRealValue(r->GetName())) {
            delete it;
            return kFALSE;
         }
      delete it;
      return kTRUE;
   }
   return kFALSE;
}

// Populate xPrime with a new proposed point
void PdfProposal::Propose(RooArgSet& xPrime, RooArgSet& x)
{
   if (fLastX.getSize() == 0) {
      // fLastX not yet initialized
      fLastX.addClone(x);
      // generate initial cache
      RooStats::SetParameters(&x, &fMaster);
      if (fMap.size() > 0) {
         for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
            fIt->first->setVal(fIt->second->getVal(&x));
      }
      fCache = fPdf->generate(xPrime, fCacheSize);
   }

   Bool_t moved = false;
   if (fMap.size() > 0) {
      moved = !Equals(fLastX, x);

      // if we've moved, set the values of the variables in the PDF to the
      // corresponding values of the variables in x, according to the
      // mappings (i.e. let the variables in x set the given values for the
      // PDF that will generate xPrime)
      if (moved) {
         // update the pdf parameters
         RooStats::SetParameters(&x, &fMaster);

         for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
            fIt->first->setVal(fIt->second->getVal(&x));

         // save the new x in fLastX
         RooStats::SetParameters(&x, &fLastX);
      }
   }

   // generate new cache if necessary
   if (moved || fCachePosition >= fCacheSize) {
      delete fCache;
      fCache = fPdf->generate(xPrime, fCacheSize);
      fCachePosition = 0;
   }

   const RooArgSet* proposal = fCache->get(fCachePosition);
   fCachePosition++;
   RooStats::SetParameters(proposal, &xPrime);
}

// Determine whether or not the proposal density is symmetric for
// points x1 and x2 - that is, whether the probabilty of reaching x2
// from x1 is equal to the probability of reaching x1 from x2
Bool_t PdfProposal::IsSymmetric(RooArgSet& /* x1 */, RooArgSet& /* x2 */)
{
   // kbelasco: is there a better way to do this?
   return false;
}

// Return the probability of proposing the point x1 given the starting
// point x2
Double_t PdfProposal::GetProposalDensity(RooArgSet& x1, RooArgSet& x2)
{
   RooStats::SetParameters(&x2, &fMaster);
   for (fIt = fMap.begin(); fIt != fMap.end(); fIt++)
      fIt->first->setVal(fIt->second->getVal(&x2));
   RooArgSet* temp = fPdf->getObservables(x1);
   RooStats::SetParameters(&x1, temp);
   delete temp;
   return fPdf->getVal(&x1); // could just as well use x2
}

void PdfProposal::AddMapping(RooRealVar& proposalParam, RooAbsReal& update)
{
   fMaster.add(*update.getParameters((RooAbsData*)NULL));
   if (update.getParameters((RooAbsData*)NULL)->getSize() == 0)
      fMaster.add(update);
   fMap.insert(pair<RooRealVar*, RooAbsReal*>(&proposalParam, &update));
}
