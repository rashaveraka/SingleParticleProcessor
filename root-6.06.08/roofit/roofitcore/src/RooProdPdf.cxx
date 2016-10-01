/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// RooProdPdf is an efficient implementation of a product of PDFs of the form
//
//  PDF_1 * PDF_2 * ... * PDF_N
//
// PDFs may share observables. If that is the case any irreducable subset
// of PDFS that share observables will be normalized with explicit numeric
// integration as any built-in normalization will no longer be valid.
//
// Alternatively, products using conditional PDFs can be defined, e.g.
//
//    F(x|y) * G(y)
//
// meaning a pdf F(x) _given_ y and a PDF G(y). In this contruction F is only
// normalized w.r.t x and G is normalized w.r.t y. The product in this construction
// is properly normalized.
//
// If exactly one of the component PDFs supports extended likelihood fits, the
// product will also be usable in extended mode, returning the number of expected
// events from the extendable component PDF. The extendable component does not
// have to appear in any specific place in the list.
//
// END_HTML
//

#include "RooProdPdf.h"
#include "RooRealProxy.h"
#include "RooProdGenContext.h"
#include "RooGenProdProj.h"
#include "RooProduct.h"
#include "RooNameReg.h"
#include "RooMsgService.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooAddition.h"
#include "RooGlobalFunc.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooRangeBoolean.h"
#include "RooCustomizer.h"
#include "RooRealIntegral.h"
#include "RooTrace.h"

#include <cstring>
#include <sstream>
#include <algorithm>

#ifndef _WIN32
#include <strings.h>
#else

static char *strtok_r(char *s1, const char *s2, char **lasts)
{
  char *ret;

  if (s1 == NULL)
    s1 = *lasts;
  while(*s1 && strchr(s2, *s1))
    ++s1;
  if(*s1 == '\0')
    return NULL;
  ret = s1;
  while(*s1 && !strchr(s2, *s1))
    ++s1;
  if(*s1)
    *s1++ = '\0';
  *lasts = s1;
  return ret;
}

#endif


#include "TSystem.h"

using namespace std;

ClassImp(RooProdPdf)
;



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooProdPdf::RooProdPdf() :
  _curNormSet(0),
  _cutOff(0),
  _extendedIndex(-1),
  _useDefaultGen(kFALSE),
  _refRangeName(0),
  _selfNorm(kTRUE)
{
  // Default constructor
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Dummy constructor

RooProdPdf::RooProdPdf(const char *name, const char *title, Double_t cutOff) :
  RooAbsPdf(name,title),
  _cacheMgr(this,10),
  _genCode(10),
  _cutOff(cutOff),
  _pdfList("!pdfs","List of PDFs",this),
  _extendedIndex(-1),
  _useDefaultGen(kFALSE),
  _refRangeName(0),
  _selfNorm(kTRUE)
{
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor with 2 PDFs (most frequent use case).
/// 
/// The optional cutOff parameter can be used as a speed optimization if
/// one or more of the PDF have sizable regions with very small values,
/// which would pull the entire product of PDFs to zero in those regions.
///
/// After each PDF multiplication, the running product is compared with
/// the cutOff parameter. If the running product is smaller than the
/// cutOff value, the product series is terminated and remaining PDFs
/// are not evaluated.
///
/// There is no magic value of the cutOff, the user should experiment
/// to find the appropriate balance between speed and precision.
/// If a cutoff is specified, the PDFs most likely to be small should
/// be put first in the product. The default cutOff value is zero.
///

RooProdPdf::RooProdPdf(const char *name, const char *title,
		       RooAbsPdf& pdf1, RooAbsPdf& pdf2, Double_t cutOff) :
  RooAbsPdf(name,title),
  _cacheMgr(this,10),
  _genCode(10),
  _cutOff(cutOff),
  _pdfList("!pdfs","List of PDFs",this),
  _extendedIndex(-1),
  _useDefaultGen(kFALSE),
  _refRangeName(0),
  _selfNorm(kTRUE)
{
  _pdfList.add(pdf1) ;
  RooArgSet* nset1 = new RooArgSet("nset") ;
  _pdfNSetList.Add(nset1) ;
  if (pdf1.canBeExtended()) {
    _extendedIndex = _pdfList.index(&pdf1) ;
  }

  _pdfList.add(pdf2) ;
  RooArgSet* nset2 = new RooArgSet("nset") ;
  _pdfNSetList.Add(nset2) ;

  if (pdf2.canBeExtended()) {
    if (_extendedIndex>=0) {
      // Protect against multiple extended terms
      coutW(InputArguments) << "RooProdPdf::RooProdPdf(" << GetName()
			    << ") multiple components with extended terms detected,"
			    << " product will not be extendible." << endl ;
      _extendedIndex=-1 ;
    } else {
      _extendedIndex=_pdfList.index(&pdf2) ;
    }
  }
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor from a list of PDFs
/// 
/// The optional cutOff parameter can be used as a speed optimization if
/// one or more of the PDF have sizable regions with very small values,
/// which would pull the entire product of PDFs to zero in those regions.
///
/// After each PDF multiplication, the running product is compared with
/// the cutOff parameter. If the running product is smaller than the
/// cutOff value, the product series is terminated and remaining PDFs
/// are not evaluated.
///
/// There is no magic value of the cutOff, the user should experiment
/// to find the appropriate balance between speed and precision.
/// If a cutoff is specified, the PDFs most likely to be small should
/// be put first in the product. The default cutOff value is zero.

RooProdPdf::RooProdPdf(const char* name, const char* title, const RooArgList& inPdfList, Double_t cutOff) :
  RooAbsPdf(name,title),
  _cacheMgr(this,10),
  _genCode(10),
  _cutOff(cutOff),
  _pdfList("!pdfs","List of PDFs",this),
  _extendedIndex(-1),
  _useDefaultGen(kFALSE),
  _refRangeName(0),
  _selfNorm(kTRUE)
{
  RooFIter iter = inPdfList.fwdIterator();
  RooAbsArg* arg ;
  Int_t numExtended(0) ;
  while((arg=(RooAbsArg*)iter.next())) {
    RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(arg) ;
    if (!pdf) {
      coutW(InputArguments) << "RooProdPdf::RooProdPdf(" << GetName() << ") list arg "
			    << arg->GetName() << " is not a PDF, ignored" << endl ;
      continue ;
    }
    _pdfList.add(*pdf) ;

    RooArgSet* nset = new RooArgSet("nset") ;
    _pdfNSetList.Add(nset) ;

    if (pdf->canBeExtended()) {
      _extendedIndex = _pdfList.index(pdf) ;
      numExtended++ ;
    }
  }

  // Protect against multiple extended terms
  if (numExtended>1) {
    coutW(InputArguments) << "RooProdPdf::RooProdPdf(" << GetName()
			  << ") WARNING: multiple components with extended terms detected,"
			  << " product will not be extendible." << endl ;
    _extendedIndex = -1 ;
  }

  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor from named argument list
///
/// fullPdf -- Set of 'regular' PDFS that are normalized over all their observables
/// Conditional(pdfSet,depSet) -- Add PDF to product with condition that it
///                               only be normalized over specified observables
///                               any remaining observables will be conditional
///                               observables
///                               
///
/// For example, given a PDF F(x,y) and G(y)
///
/// RooProdPdf("P","P",G,Conditional(F,x)) will construct a 2-dimensional PDF as follows:
/// 
///   P(x,y) = G(y)/Int[y]G(y) * F(x,y)/Int[x]G(x,y)
///
/// which is a well normalized and properly defined PDF, but different from the
///  
///  P'(x,y) = F(x,y)*G(y) / Int[x,y] F(x,y)*G(y)
///
/// In the former case the y distribution of P is identical to that of G, while
/// F only is used to determine the correlation between X and Y. In the latter
/// case the Y distribution is defined by the product of F and G.
///
/// This P(x,y) construction is analoguous to generating events from F(x,y) with
/// a prototype dataset sampled from G(y)

RooProdPdf::RooProdPdf(const char* name, const char* title, const RooArgSet& fullPdfSet,
		       const RooCmdArg& arg1, const RooCmdArg& arg2,
		       const RooCmdArg& arg3, const RooCmdArg& arg4,
		       const RooCmdArg& arg5, const RooCmdArg& arg6,
		       const RooCmdArg& arg7, const RooCmdArg& arg8) :
  RooAbsPdf(name,title),
  _cacheMgr(this,10),
  _genCode(10),
  _cutOff(0),
  _pdfList("!pdfs","List of PDFs",this),
  _extendedIndex(-1),
  _useDefaultGen(kFALSE),
  _refRangeName(0),
  _selfNorm(kTRUE)
{
  RooLinkedList l ;
  l.Add((TObject*)&arg1) ;  l.Add((TObject*)&arg2) ;
  l.Add((TObject*)&arg3) ;  l.Add((TObject*)&arg4) ;
  l.Add((TObject*)&arg5) ;  l.Add((TObject*)&arg6) ;
  l.Add((TObject*)&arg7) ;  l.Add((TObject*)&arg8) ;

  initializeFromCmdArgList(fullPdfSet,l) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor from named argument list
///
/// fullPdf -- Set of 'regular' PDFS that are normalized over all their observables
/// Conditional(pdfSet,depSet) -- Add PDF to product with condition that it
///                               only be normalized over specified observables
///                               any remaining observables will be conditional
///                               observables
///                               
///
/// For example, given a PDF F(x,y) and G(y)
///
/// RooProdPdf("P","P",G,Conditional(F,x)) will construct a 2-dimensional PDF as follows:
/// 
///   P(x,y) = G(y)/Int[y]G(y) * F(x,y)/Int[x]G(x,y)
///
/// which is a well normalized and properly defined PDF, but different from the
///  
///  P'(x,y) = F(x,y)*G(y) / Int[x,y] F(x,y)*G(y)
///
/// In the former case the y distribution of P is identical to that of G, while
/// F only is used to determine the correlation between X and Y. In the latter
/// case the Y distribution is defined by the product of F and G.
///
/// This P(x,y) construction is analoguous to generating events from F(x,y) with
/// a prototype dataset sampled from G(y)

RooProdPdf::RooProdPdf(const char* name, const char* title,
		       const RooCmdArg& arg1, const RooCmdArg& arg2,
		       const RooCmdArg& arg3, const RooCmdArg& arg4,
		       const RooCmdArg& arg5, const RooCmdArg& arg6,
		       const RooCmdArg& arg7, const RooCmdArg& arg8) :
  RooAbsPdf(name,title),
  _cacheMgr(this,10),
  _genCode(10),
  _cutOff(0),
  _pdfList("!pdfList","List of PDFs",this),
  _extendedIndex(-1),
  _useDefaultGen(kFALSE),
  _refRangeName(0),
  _selfNorm(kTRUE)
{
  RooLinkedList l ;
  l.Add((TObject*)&arg1) ;  l.Add((TObject*)&arg2) ;
  l.Add((TObject*)&arg3) ;  l.Add((TObject*)&arg4) ;
  l.Add((TObject*)&arg5) ;  l.Add((TObject*)&arg6) ;
  l.Add((TObject*)&arg7) ;  l.Add((TObject*)&arg8) ;

  initializeFromCmdArgList(RooArgSet(),l) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Internal constructor from list of named arguments  

RooProdPdf::RooProdPdf(const char* name, const char* title, const RooArgSet& fullPdfSet, const RooLinkedList& cmdArgList) :
  RooAbsPdf(name,title),
  _cacheMgr(this,10),
  _genCode(10),
  _cutOff(0),
  _pdfList("!pdfs","List of PDFs",this),
  _extendedIndex(-1),
  _useDefaultGen(kFALSE),
  _refRangeName(0),
  _selfNorm(kTRUE)
{
  initializeFromCmdArgList(fullPdfSet, cmdArgList) ;
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooProdPdf::RooProdPdf(const RooProdPdf& other, const char* name) :
  RooAbsPdf(other,name),
  _cacheMgr(other._cacheMgr,this),
  _genCode(other._genCode),
  _cutOff(other._cutOff),
  _pdfList("!pdfs",this,other._pdfList),
  _extendedIndex(other._extendedIndex),
  _useDefaultGen(other._useDefaultGen),
  _refRangeName(other._refRangeName),
  _selfNorm(other._selfNorm),
  _defNormSet(other._defNormSet)
{
  // Clone contents of normalizarion set list
  RooFIter iter = other._pdfNSetList.fwdIterator();
  RooArgSet* nset ;
  while((nset=(RooArgSet*)iter.next())) {
    RooArgSet* tmp = (RooArgSet*) nset->snapshot() ;
    tmp->setName(nset->GetName()) ;
    _pdfNSetList.Add(tmp) ;
  }
  TRACE_CREATE
}



////////////////////////////////////////////////////////////////////////////////
/// Initialize RooProdPdf configuration from given list of RooCmdArg configuration arguments
/// and set of 'regular' p.d.f.s in product

void RooProdPdf::initializeFromCmdArgList(const RooArgSet& fullPdfSet, const RooLinkedList& l)
{
  Int_t numExtended(0) ;

  // Process set of full PDFS
  RooFIter siter = fullPdfSet.fwdIterator() ;
  RooAbsPdf* pdf ;
  while((pdf=(RooAbsPdf*)siter.next())) {
    _pdfList.add(*pdf) ;
    RooArgSet* nset1 = new RooArgSet("nset") ;
    _pdfNSetList.Add(nset1) ;

    if (pdf->canBeExtended()) {
      _extendedIndex = _pdfList.index(pdf) ;
      numExtended++ ;
    }

  }

  // Process list of conditional PDFs
  RooFIter iter = l.fwdIterator();
  RooCmdArg* carg ;
  while((carg=(RooCmdArg*)iter.next())) {

    if (0 == strcmp(carg->GetName(), "Conditional")) {

      Int_t argType = carg->getInt(0) ;
      RooArgSet* pdfSet = (RooArgSet*) carg->getSet(0) ;
      RooArgSet* normSet = (RooArgSet*) carg->getSet(1) ;

      RooFIter siter2 = pdfSet->fwdIterator() ;
      RooAbsPdf* thePdf ;
      while ((thePdf=(RooAbsPdf*)siter2.next())) {
	_pdfList.add(*thePdf) ;

	RooArgSet* tmp = (RooArgSet*) normSet->snapshot() ;
	tmp->setName(0 == argType ? "nset" : "cset") ;
	_pdfNSetList.Add(tmp) ;       	

	if (thePdf->canBeExtended()) {
	  _extendedIndex = _pdfList.index(thePdf) ;
	  numExtended++ ;
	}

      }

    } else if (0 != strlen(carg->GetName())) {
      coutW(InputArguments) << "Unknown arg: " << carg->GetName() << endl ;
    }
  }

  // Protect against multiple extended terms
  if (numExtended>1) {
    coutW(InputArguments) << "RooProdPdf::RooProdPdf(" << GetName()
			  << ") WARNING: multiple components with extended terms detected,"
			  << " product will not be extendible." << endl ;
    _extendedIndex = -1 ;
  }


}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooProdPdf::~RooProdPdf()
{
  _pdfNSetList.Delete() ;
  TRACE_DESTROY
}



////////////////////////////////////////////////////////////////////////////////
/// Overload getVal() to intercept normalization set for use in evaluate()

Double_t RooProdPdf::getValV(const RooArgSet* set) const 
{
  _curNormSet = (RooArgSet*)set ;
  return RooAbsPdf::getValV(set) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate current value of object

Double_t RooProdPdf::evaluate() const 
{
  Int_t code ;
  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(_curNormSet,0,&code) ;

  // If cache doesn't have our configuration, recalculate here
  if (!cache) {
    RooArgList *plist(0) ;
    RooLinkedList *nlist(0) ;
    getPartIntList(_curNormSet,0,plist,nlist,code) ;
    cache = (CacheElem*) _cacheMgr.getObj(_curNormSet,0,&code) ;
  }


  return calculate(*cache) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate running product of pdfs terms, using the supplied
/// normalization set in 'normSetList' for each component

Double_t RooProdPdf::calculate(const RooArgList* partIntList, const RooLinkedList* normSetList) const
{
  RooAbsReal* partInt;
  RooArgSet* normSet;
  Double_t value = 1.0;
  RooFIter plIter = partIntList->fwdIterator(), nlIter = normSetList->fwdIterator();

    for (partInt = (RooAbsReal*) plIter.next(),
	normSet = (RooArgSet*) nlIter.next(); partInt && normSet;
      partInt = (RooAbsReal*) plIter.next(),
      normSet = (RooArgSet*) nlIter.next()) {
    const Double_t piVal = partInt->getVal(normSet->getSize() > 0 ? normSet : 0);
//  cout << "RooProdPdf::calculate(" << GetName() << ") partInt(" << partInt->GetName() << ") = " << piVal << " normSet = " << normSet << " " << (normSet->getSize()>0 ? *normSet : RooArgSet()) << endl ;
    value *= piVal;
    if (value <= _cutOff) { break; }
  }

  return value ;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate running product of pdfs terms, using the supplied
/// normalization set in 'normSetList' for each component

Double_t RooProdPdf::calculate(const RooProdPdf::CacheElem& cache, Bool_t /*verbose*/) const
{
  //cout << "RooProdPdf::calculate from cache" << endl ;

  if (cache._isRearranged) {
    if (dologD(Eval)) {
      cxcoutD(Eval) << "RooProdPdf::calculate(" << GetName() << ") rearranged product calculation"
                    << " calculate: num = " << cache._rearrangedNum->GetName() << " = " << cache._rearrangedNum->getVal() << endl ;
//       cache._rearrangedNum->printComponentTree("",0,5) ;
      cxcoutD(Eval) << "calculate: den = " << cache._rearrangedDen->GetName() << " = " << cache._rearrangedDen->getVal() << endl ;
//       cache._rearrangedDen->printComponentTree("",0,5) ;
    }

    return cache._rearrangedNum->getVal() / cache._rearrangedDen->getVal();
  } else {
    RooAbsReal* partInt;
    RooArgSet* normSet;
    RooFIter plIter = cache._partList.fwdIterator();
    RooFIter nlIter = cache._normList.fwdIterator();
    Double_t value = 1.0;
    for (partInt = (RooAbsReal*) plIter.next(),
	normSet = (RooArgSet*) nlIter.next(); partInt && normSet;
      partInt = (RooAbsReal*) plIter.next(),
      normSet = (RooArgSet*) nlIter.next()) {
      const Double_t piVal = partInt->getVal(normSet->getSize() > 0 ? normSet : 0);
      value *= piVal ;
      if (value <= _cutOff) break;
    }
//   cout << "return value = " << value << endl ;
    return value ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Factorize product in irreducible terms for given choice of integration/normalization

void RooProdPdf::factorizeProduct(const RooArgSet& normSet, const RooArgSet& intSet,
				  RooLinkedList& termList, RooLinkedList& normList,
				  RooLinkedList& impDepList, RooLinkedList& crossDepList,
				  RooLinkedList& intList) const
{
  // List of all term dependents: normalization and imported
  RooLinkedList depAllList;
  RooLinkedList depIntNoNormList;

  // Setup lists for factorization terms and their dependents
  RooArgSet* term(0);
  RooArgSet* termNormDeps(0);
  RooArgSet* termAllDeps(0);
  RooArgSet* termIntDeps(0);
  RooArgSet* termIntNoNormDeps(0);

  // Loop over the PDFs
  RooAbsPdf* pdf;
  RooArgSet* pdfNSetOrig;
  for (RooFIter pdfIter = _pdfList.fwdIterator(),
      nIter = _pdfNSetList.fwdIterator();
      (pdfNSetOrig = (RooArgSet*) nIter.next(),
       pdf = (RooAbsPdf*) pdfIter.next()); ) {
    RooArgSet* pdfNSet, *pdfCSet;

    // Reduce pdfNSet to actual dependents
    if (0 == strcmp("nset", pdfNSetOrig->GetName())) {
      pdfNSet = pdf->getObservables(*pdfNSetOrig);
      pdfCSet = new RooArgSet;
    } else if (0 == strcmp("cset", pdfNSetOrig->GetName())) {
      RooArgSet* tmp = pdf->getObservables(normSet);
      tmp->remove(*pdfNSetOrig, kTRUE, kTRUE);
      pdfCSet = pdfNSetOrig;
      pdfNSet = tmp;
    } else {
      // Legacy mode. Interpret at NSet for backward compatibility
      pdfNSet = pdf->getObservables(*pdfNSetOrig);
      pdfCSet = new RooArgSet;
    }


    RooArgSet pdfNormDeps; // Dependents to be normalized for the PDF
    RooArgSet pdfAllDeps; // All dependents of this PDF

    // Make list of all dependents of this PDF
    RooArgSet* tmp = pdf->getObservables(normSet);
    pdfAllDeps.add(*tmp);
    delete tmp;


//     cout << GetName() << ": pdf = " << pdf->GetName() << " pdfAllDeps = " << pdfAllDeps << " pdfNSet = " << *pdfNSet << " pdfCSet = " << *pdfCSet << endl;

    // Make list of normalization dependents for this PDF;
    if (pdfNSet->getSize() > 0) {
      // PDF is conditional
      RooArgSet* tmp2 = (RooArgSet*) pdfAllDeps.selectCommon(*pdfNSet);
      pdfNormDeps.add(*tmp2);
      delete tmp2;
    } else {
      // PDF is regular
      pdfNormDeps.add(pdfAllDeps);
    }

//     cout << GetName() << ": pdfNormDeps for " << pdf->GetName() << " = " << pdfNormDeps << endl;

    RooArgSet* pdfIntSet = pdf->getObservables(intSet) ;

    // WVE if we have no norm deps, conditional observables should be taken out of pdfIntSet
    if (0 == pdfNormDeps.getSize() && pdfCSet->getSize() > 0) {
      pdfIntSet->remove(*pdfCSet, kTRUE, kTRUE);
//       cout << GetName() << ": have no norm deps, removing conditional observables from intset" << endl;
    }

    RooArgSet pdfIntNoNormDeps(*pdfIntSet);
    pdfIntNoNormDeps.remove(pdfNormDeps, kTRUE, kTRUE);

//     cout << GetName() << ": pdf = " << pdf->GetName() << " intset = " << *pdfIntSet << " pdfIntNoNormDeps = " << pdfIntNoNormDeps << endl;

    // Check if this PDF has dependents overlapping with one of the existing terms
    Bool_t done(kFALSE);
    for (RooFIter lIter = termList.fwdIterator(),
	ldIter = normList.fwdIterator(),
	laIter = depAllList.fwdIterator();
      (termNormDeps = (RooArgSet*) ldIter.next(),
       termAllDeps = (RooArgSet*) laIter.next(),
       term = (RooArgSet*) lIter.next()); ) {
      // PDF should be added to existing term if
      // 1) It has overlapping normalization dependents with any other PDF in existing term
      // 2) It has overlapping dependents of any class for which integration is requested
      // 3) If normalization happens over multiple ranges, and those ranges are both defined
      //    in either observable

      Bool_t normOverlap = pdfNormDeps.overlaps(*termNormDeps);
      //Bool_t intOverlap =  pdfIntSet->overlaps(*termAllDeps);

      if (normOverlap) {
//  	cout << GetName() << ": this term overlaps with term " << (*term) << " in normalization observables" << endl;

	term->add(*pdf);
	termNormDeps->add(pdfNormDeps, kFALSE);
	termAllDeps->add(pdfAllDeps, kFALSE);
	if (!termIntDeps) {
	  termIntDeps = new RooArgSet("termIntDeps");
	}
	termIntDeps->add(*pdfIntSet, kFALSE);
	if (!termIntNoNormDeps) {
	  termIntNoNormDeps = new RooArgSet("termIntNoNormDeps");
	}
	termIntNoNormDeps->add(pdfIntNoNormDeps, kFALSE);
	done = kTRUE;
      }
    }

    // If not, create a new term
    if (!done) {
      if (!(0 == pdfNormDeps.getSize() && 0 == pdfAllDeps.getSize() &&
	    0 == pdfIntSet->getSize()) || 0 == normSet.getSize()) {
//   	cout << GetName() << ": creating new term" << endl;
	term = new RooArgSet("term");
	termNormDeps = new RooArgSet("termNormDeps");
	termAllDeps = new RooArgSet("termAllDeps");
	termIntDeps = new RooArgSet("termIntDeps");
	termIntNoNormDeps = new RooArgSet("termIntNoNormDeps");
	
	term->add(*pdf);
	termNormDeps->add(pdfNormDeps, kFALSE);
	termAllDeps->add(pdfAllDeps, kFALSE);
	termIntDeps->add(*pdfIntSet, kFALSE);
	termIntNoNormDeps->add(pdfIntNoNormDeps, kFALSE);
	
	termList.Add(term);
	normList.Add(termNormDeps);
	depAllList.Add(termAllDeps);
	intList.Add(termIntDeps);
	depIntNoNormList.Add(termIntNoNormDeps);
      }
    }

    // We own the reduced version of pdfNSet
    delete pdfNSet;
    delete pdfIntSet;
    if (pdfCSet != pdfNSetOrig) {
      delete pdfCSet;
    }
  }

  // Loop over list of terms again to determine 'imported' observables
  RooArgSet *normDeps, *allDeps, *intNoNormDeps;
  for (RooFIter lIter = termList.fwdIterator(),
      ldIter = normList.fwdIterator(),
      laIter = depAllList.fwdIterator(),
      innIter = depIntNoNormList.fwdIterator();
      (normDeps = (RooArgSet*) ldIter.next(),
       allDeps = (RooArgSet*) laIter.next(),
       intNoNormDeps = (RooArgSet*) innIter.next(),
       term=(RooArgSet*)lIter.next()); ) {
    // Make list of wholly imported dependents
    RooArgSet impDeps(*allDeps);
    impDeps.remove(*normDeps, kTRUE, kTRUE);
    impDepList.Add(impDeps.snapshot());
//     cout << GetName() << ": list of imported dependents for term " << (*term) << " set to " << impDeps << endl ;

    // Make list of cross dependents (term is self contained for these dependents,
    // but components import dependents from other components)
    RooArgSet* crossDeps = (RooArgSet*) intNoNormDeps->selectCommon(*normDeps);
    crossDepList.Add(crossDeps->snapshot());
//     cout << GetName() << ": list of cross dependents for term " << (*term) << " set to " << *crossDeps << endl ;
    delete crossDeps;
  }

  depAllList.Delete();
  depIntNoNormList.Delete();

  return;
}




////////////////////////////////////////////////////////////////////////////////
/// Return list of (partial) integrals of product terms for integration
/// of p.d.f over observables iset while normalization over observables nset.
/// Also return list of normalization sets to be used to evaluate 
/// each component in the list correctly.

void RooProdPdf::getPartIntList(const RooArgSet* nset, const RooArgSet* iset, 
				pRooArgList& partList, pRooLinkedList& nsetList, 
				Int_t& code, const char* isetRangeName) const 
{
//    cout << "   FOLKERT::RooProdPdf::getPartIntList(" << GetName() <<")  nset = " << (nset?*nset:RooArgSet()) << endl 
//         << "   _normRange = " << _normRange << endl 
//         << "   iset = " << (iset?*iset:RooArgSet()) << endl 
//         << "   isetRangeName = " << (isetRangeName?isetRangeName:"<null>") << endl ;

  // Check if this configuration was created before
  Int_t sterileIdx(-1);

  CacheElem* cache = (CacheElem*) _cacheMgr.getObj(nset,iset,&sterileIdx,isetRangeName);
  if (cache) {
    code = _cacheMgr.lastIndex();
    partList = &cache->_partList;
    nsetList = &cache->_normList;
    return;
  }

  // Create containers for partial integral components to be generated
  cache = new CacheElem;

  // Factorize the product in irreducible terms for this nset
  RooLinkedList terms, norms, imp, ints, cross;
//   cout << "RooProdPdf::getPIL -- now calling factorizeProduct()" << endl ;


  // Normalization set used for factorization
  RooArgSet factNset(nset ? (*nset) : _defNormSet);
//   cout << GetName() << "factNset = " << factNset << endl ;

  factorizeProduct(factNset, iset ? (*iset) : RooArgSet(), terms, norms, imp, cross, ints);

  RooArgSet *norm, *integ, *xdeps, *imps;

  // Group irriducible terms that need to be (partially) integrated together
  RooLinkedList groupedList;
  RooArgSet outerIntDeps;
//   cout << "RooProdPdf::getPIL -- now calling groupProductTerms()" << endl;
  groupProductTerms(groupedList, outerIntDeps, terms, norms, imp, ints, cross);
  RooFIter gIter = groupedList.fwdIterator();
  RooLinkedList* group;

  // Loop over groups
//   cout<<"FK: pdf("<<GetName()<<") Starting selecting F(x|y)!"<<endl;
  // Find groups of type F(x|y), i.e. termImpSet!=0, construct ratio object
  map<string, RooArgSet> ratioTerms;
  while ((group = (RooLinkedList*) gIter.next())) {
    if (1 == group->GetSize()) {
//       cout<<"FK: Starting Single Term"<<endl;

      RooArgSet* term = (RooArgSet*) group->At(0);
      	
      Int_t termIdx = terms.IndexOf(term);
      norm=(RooArgSet*) norms.At(termIdx);
      imps=(RooArgSet*)imp.At(termIdx);
      RooArgSet termNSet(*norm), termImpSet(*imps);

//       cout<<"FK: termImpSet.getSize()  = "<<termImpSet.getSize()<< " " << termImpSet << endl;
//       cout<<"FK: _refRangeName = "<<_refRangeName<<endl;

      if (termImpSet.getSize() > 0 && 0 != _refRangeName) {	

// 	cout << "WVE now here" << endl;

	// WVE we can skip this if the ref range is equal to the normalization range
	Bool_t rangeIdentical(kTRUE);
	RooFIter niter = termNSet.fwdIterator();
	RooRealVar* normObs;
// 	cout << "_normRange = " << _normRange << " _refRangeName = " << RooNameReg::str(_refRangeName) << endl ;
	while ((normObs = (RooRealVar*) niter.next())) {
	  //FK: Here the refRange should be compared to _normRange, if it's set, and to the normObs range if it's not set
	  if (_normRange.Length() > 0) {	
	    if (normObs->getMin(_normRange.Data()) != normObs->getMin(RooNameReg::str(_refRangeName))) rangeIdentical = kFALSE;
	    if (normObs->getMax(_normRange.Data()) != normObs->getMax(RooNameReg::str(_refRangeName))) rangeIdentical = kFALSE;	
	  }
	  else{
	    if (normObs->getMin() != normObs->getMin(RooNameReg::str(_refRangeName))) rangeIdentical = kFALSE;
	    if (normObs->getMax() != normObs->getMax(RooNameReg::str(_refRangeName))) rangeIdentical = kFALSE;
	  }
	}
// 	cout<<"FK: rangeIdentical Single = "<<(rangeIdentical ? 'T':'F')<<endl;
	// coverity[CONSTANT_EXPRESSION_RESULT]
	if (!rangeIdentical || 1) {
// 	  cout << "PREPARING RATIO HERE (SINGLE TERM)" << endl ;	
	  RooAbsReal* ratio = makeCondPdfRatioCorr(*(RooAbsReal*)term->first(), termNSet, termImpSet, normRange(), RooNameReg::str(_refRangeName));
	  ostringstream str; termImpSet.printValue(str);
// 	  cout << GetName() << "inserting ratio term" << endl;
	  ratioTerms[str.str()].add(*ratio);
	}	
      }
	
    } else {
//       cout<<"FK: Starting Composite Term"<<endl;

      RooArgSet compTermSet, compTermNorm;
      RooFIter tIter = group->fwdIterator();
      RooArgSet* term;
      while ((term = (RooArgSet*) tIter.next())) {
	
	Int_t termIdx = terms.IndexOf(term);
	norm=(RooArgSet*) norms.At(termIdx);
	imps=(RooArgSet*)imp.At(termIdx);	
	RooArgSet termNSet(*norm), termImpSet(*imps);
	
	if (termImpSet.getSize() > 0 && 0 != _refRangeName) {	

	  // WVE we can skip this if the ref range is equal to the normalization range
	  Bool_t rangeIdentical(kTRUE);
	  RooFIter niter = termNSet.fwdIterator();
	  RooRealVar* normObs;
	  //FK: Here the refRange should be compared to _normRange, if it's set, and to the normObs range if it's not set
	  if(_normRange.Length() > 0) {
	    while ((normObs = (RooRealVar*) niter.next())) {
	      if (normObs->getMin(_normRange.Data()) != normObs->getMin(RooNameReg::str(_refRangeName))) rangeIdentical = kFALSE;
	      if (normObs->getMax(_normRange.Data()) != normObs->getMax(RooNameReg::str(_refRangeName))) rangeIdentical = kFALSE;
	    }
	  } else {
	    while ((normObs = (RooRealVar*) niter.next())) {
	      if (normObs->getMin() != normObs->getMin(RooNameReg::str(_refRangeName))) rangeIdentical = kFALSE;
	      if (normObs->getMax() != normObs->getMax(RooNameReg::str(_refRangeName))) rangeIdentical = kFALSE;
	    }
	  }
// 	  cout<<"FK: rangeIdentical Composite = "<<(rangeIdentical ? 'T':'F') <<endl;
	  if (!rangeIdentical || 1) {	
// 	    cout << "PREPARING RATIO HERE (COMPOSITE TERM)" << endl ;	
	    RooAbsReal* ratio = makeCondPdfRatioCorr(*(RooAbsReal*)term->first(), termNSet, termImpSet, normRange(), RooNameReg::str(_refRangeName));
	    ostringstream str; termImpSet.printValue(str);
	    ratioTerms[str.str()].add(*ratio);
	  }
	}
      }
    }

  }

  // Find groups with y as termNSet
  // Replace G(y) with (G(y),ratio)
  gIter = groupedList.fwdIterator();
  while ((group = (RooLinkedList*) gIter.next())) {
    if (1 == group->GetSize()) {
      RooArgSet* term = (RooArgSet*) group->At(0);

      Int_t termIdx = terms.IndexOf(term);
      norm = (RooArgSet*) norms.At(termIdx);
      imps = (RooArgSet*) imp.At(termIdx);
      RooArgSet termNSet(*norm), termImpSet(*imps);

      // If termNset matches index of ratioTerms, insert ratio here
      ostringstream str; termNSet.printValue(str);
      if (ratioTerms[str.str()].getSize() > 0) {
//  	cout << "MUST INSERT RATIO OBJECT IN TERM (SINGLE) " << *term << endl;
	term->add(ratioTerms[str.str()]);
      }
    } else {
      RooArgSet compTermSet, compTermNorm;
      RooFIter tIter = group->fwdIterator();
      RooArgSet* term;
      while ((term = (RooArgSet*) tIter.next())) {
	Int_t termIdx = terms.IndexOf(term);
	norm = (RooArgSet*) norms.At(termIdx);
	imps = (RooArgSet*) imp.At(termIdx);	
	RooArgSet termNSet(*norm), termImpSet(*imps);
	
	// If termNset matches index of ratioTerms, insert ratio here
	ostringstream str; termNSet.printValue(str);
	if (ratioTerms[str.str()].getSize() > 0) {
//  	  cout << "MUST INSERT RATIO OBJECT IN TERM (COMPOSITE)" << *term << endl;
	  term->add(ratioTerms[str.str()]);
	}
      }
    }
  }

  gIter = groupedList.fwdIterator();
  while ((group = (RooLinkedList*) gIter.next())) {
//     cout << GetName() << ":now processing group" << endl;
//      group->Print("1");

    if (1 == group->GetSize()) {
//       cout << "processing atomic item" << endl;
      RooArgSet* term = (RooArgSet*) group->At(0);

      Int_t termIdx = terms.IndexOf(term);
      norm = (RooArgSet*) norms.At(termIdx);
      integ = (RooArgSet*) ints.At(termIdx);
      xdeps = (RooArgSet*) cross.At(termIdx);
      imps = (RooArgSet*) imp.At(termIdx);
	
      RooArgSet termNSet, termISet, termXSet, termImpSet;

      // Take list of normalization, integrated dependents from factorization algorithm
      termISet.add(*integ);
      termNSet.add(*norm);

      // Cross-imported integrated dependents
      termXSet.add(*xdeps);
      termImpSet.add(*imps);

//       cout << GetName() << ": termISet = " << termISet << endl;
//       cout << GetName() << ": termNSet = " << termNSet << endl;
//       cout << GetName() << ": termXSet = " << termXSet << endl;
//       cout << GetName() << ": termImpSet = " << termImpSet << endl;

      // Add prefab term to partIntList.
      Bool_t isOwned(kFALSE);
      vector<RooAbsReal*> func = processProductTerm(nset, iset, isetRangeName, term, termNSet, termISet, isOwned);
      if (func[0]) {
	cache->_partList.add(*func[0]);
	if (isOwned) cache->_ownedList.addOwned(*func[0]);
	
	cache->_normList.Add(norm->snapshot(kFALSE));

	cache->_numList.addOwned(*func[1]);
	cache->_denList.addOwned(*func[2]);
// 	cout << "func[0]=" << func[0]->IsA()->GetName() << "::" << func[0]->GetName() << endl;
// 	cout << "func[1]=" << func[1]->IsA()->GetName() << "::" << func[1]->GetName() << endl;
// 	cout << "func[2]=" << func[2]->IsA()->GetName() << "::" << func[2]->GetName() << endl;
      }
    } else {
//        cout << "processing composite item" << endl;
      RooArgSet compTermSet, compTermNorm, compTermNum, compTermDen;
      RooFIter tIter = group->fwdIterator();
      RooArgSet* term;
      while ((term = (RooArgSet*) tIter.next())) {
//   	cout << GetName() << ": processing term " << (*term) << " of composite item" << endl ;
	Int_t termIdx = terms.IndexOf(term);
	norm = (RooArgSet*) norms.At(termIdx);
	integ = (RooArgSet*) ints.At(termIdx);
	xdeps = (RooArgSet*) cross.At(termIdx);
	imps = (RooArgSet*) imp.At(termIdx);
	
	RooArgSet termNSet, termISet, termXSet, termImpSet;
	termISet.add(*integ);
	termNSet.add(*norm);
	termXSet.add(*xdeps);
	termImpSet.add(*imps);

	// Remove outer integration dependents from termISet
	termISet.remove(outerIntDeps, kTRUE, kTRUE);
//    	cout << "termISet = "; termISet.Print("1");

//  	cout << GetName() << ": termISet = " << termISet << endl;
//  	cout << GetName() << ": termNSet = " << termNSet << endl;
//   	cout << GetName() << ": termXSet = " << termXSet << endl;
//   	cout << GetName() << ": termImpSet = " << termImpSet << endl;

	Bool_t isOwned;
	vector<RooAbsReal*> func = processProductTerm(nset, iset, isetRangeName, term, termNSet, termISet, isOwned, kTRUE);
//    	cout << GetName() << ": created composite term component " << func[0]->GetName() << endl;
	if (func[0]) {
	  compTermSet.add(*func[0]);
	  if (isOwned) cache->_ownedList.addOwned(*func[0]);
	  compTermNorm.add(*norm, kFALSE);

	  compTermNum.add(*func[1]);
	  compTermDen.add(*func[2]);
	  //cache->_numList.add(*func[1]);
	  //cache->_denList.add(*func[2]);

// 	  cout << "func[0]=" << func[0]->IsA()->GetName() << "::" << func[0]->GetName() << endl;
// 	  cout << "func[1]=" << func[1]->IsA()->GetName() << "::" << func[1]->GetName() << endl;
// 	  cout << "func[2]=" << func[2]->IsA()->GetName() << "::" << func[2]->GetName() << endl;
	}
      }

//       cout << GetName() << ": constructing special composite product" << endl;
//       cout << GetName() << ": compTermSet = " ; compTermSet.Print("1");

      // WVE THIS NEEDS TO BE REARRANGED

      // compTermset is set van partial integrals to be multiplied
      // prodtmp = product (compTermSet)
      // inttmp = int ( prodtmp ) d (outerIntDeps) _range_isetRangeName

      const std::string prodname = makeRGPPName("SPECPROD", compTermSet, outerIntDeps, RooArgSet(), isetRangeName);
      RooProduct* prodtmp = new RooProduct(prodname.c_str(), prodname.c_str(), compTermSet);
      cache->_ownedList.addOwned(*prodtmp);

      const std::string intname = makeRGPPName("SPECINT", compTermSet, outerIntDeps, RooArgSet(), isetRangeName);
      RooRealIntegral* inttmp = new RooRealIntegral(intname.c_str(), intname.c_str(), *prodtmp, outerIntDeps, 0, 0, isetRangeName);
      inttmp->setStringAttribute("PROD_TERM_TYPE", "SPECINT");

      cache->_ownedList.addOwned(*inttmp);
      cache->_partList.add(*inttmp);

      // Product of numerator terms
      const string prodname_num = makeRGPPName("SPECPROD_NUM", compTermNum, RooArgSet(), RooArgSet(), 0);
      RooProduct* prodtmp_num = new RooProduct(prodname_num.c_str(), prodname_num.c_str(), compTermNum);
      prodtmp_num->addOwnedComponents(compTermNum);
      cache->_ownedList.addOwned(*prodtmp_num);

      // Product of denominator terms
      const string prodname_den = makeRGPPName("SPECPROD_DEN", compTermDen, RooArgSet(), RooArgSet(), 0);
      RooProduct* prodtmp_den = new RooProduct(prodname_den.c_str(), prodname_den.c_str(), compTermDen);
      prodtmp_den->addOwnedComponents(compTermDen);
      cache->_ownedList.addOwned(*prodtmp_den);

      // Ratio
      string name = Form("SPEC_RATIO(%s,%s)", prodname_num.c_str(), prodname_den.c_str());
      RooFormulaVar* ndr = new RooFormulaVar(name.c_str(), "@0/@1", RooArgList(*prodtmp_num, *prodtmp_den));

      // Integral of ratio
      RooAbsReal* numtmp = ndr->createIntegral(outerIntDeps,isetRangeName);
      numtmp->addOwnedComponents(*ndr);

      cache->_numList.addOwned(*numtmp);
      cache->_denList.addOwned(*(RooAbsArg*)RooFit::RooConst(1).clone("1"));
      cache->_normList.Add(compTermNorm.snapshot(kFALSE));
    }
  }

  // Store the partial integral list and return the assigned code
  code = _cacheMgr.setObj(nset, iset, (RooAbsCacheElement*)cache, RooNameReg::ptr(isetRangeName));

  // Fill references to be returned
  partList = &cache->_partList;
  nsetList = &cache->_normList;

  // WVE DEBUG PRINTING
//   cout << "RooProdPdf::getPartIntList(" << GetName() << ") made cache " << cache << " with the following nset pointers ";
//   TIterator* nliter = nsetList->MakeIterator();
//   RooArgSet* ns;
//   while((ns=(RooArgSet*)nliter->Next())) {
//     cout << ns << " ";
//   }
//   cout << endl;
//   delete nliter;

//   cout << "   FOLKERT::RooProdPdf::getPartIntList END(" << GetName() <<")  nset = " << (nset?*nset:RooArgSet()) << endl
//        << "   _normRange = " << _normRange << endl
//        << "   iset = " << (iset?*iset:RooArgSet()) << endl
//        << "   partList = ";
//   if(partListPointer) partListPointer->Print();
//   cout << "   nsetList = ";
//   if(nsetListPointer) nsetListPointer->Print("");
//   cout << "   code = " << code << endl
//        << "   isetRangeName = " << (isetRangeName?isetRangeName:"<null>") << endl;


  // Need to rearrange product in case of multiple ranges
  if (_normRange.Contains(",")) {
    rearrangeProduct(*cache);
  }

  // We own contents of all lists filled by factorizeProduct()
  groupedList.Delete();
  terms.Delete();
  ints.Delete();
  imp.Delete();
  norms.Delete();
  cross.Delete();
}




////////////////////////////////////////////////////////////////////////////////
/// For single normalization ranges

RooAbsReal* RooProdPdf::makeCondPdfRatioCorr(RooAbsReal& pdf, const RooArgSet& termNset, const RooArgSet& /*termImpSet*/, const char* normRangeTmp, const char* refRange) const
{
  RooAbsReal* ratio_num = pdf.createIntegral(termNset,normRangeTmp) ;
  RooAbsReal* ratio_den = pdf.createIntegral(termNset,refRange) ;
  RooFormulaVar* ratio = new RooFormulaVar(Form("ratio(%s,%s)",ratio_num->GetName(),ratio_den->GetName()),"@0/@1",
					   RooArgList(*ratio_num,*ratio_den)) ;

  ratio->addOwnedComponents(RooArgSet(*ratio_num,*ratio_den)) ;	
  ratio->setAttribute("RATIO_TERM") ;
  return ratio ;
}




////////////////////////////////////////////////////////////////////////////////

void RooProdPdf::rearrangeProduct(RooProdPdf::CacheElem& cache) const
{
  RooAbsReal* part, *num, *den ;
  RooArgSet nomList ;

  list<string> rangeComps ;
  {
    char* buf = new char[strlen(_normRange.Data()) + 1] ;
    strcpy(buf,_normRange.Data()) ;
    char* save(0) ;
    char* token = strtok_r(buf,",",&save) ;
    while(token) {
      rangeComps.push_back(token) ;
      token = strtok_r(0,",",&save) ;
    }
    delete[] buf;
  }


  map<string,RooArgSet> denListList ;
  RooArgSet specIntDeps ;
  string specIntRange ;

//   cout << "THIS IS REARRANGEPRODUCT" << endl ;

  RooFIter iterp = cache._partList.fwdIterator() ;
  RooFIter iter1 = cache._numList.fwdIterator() ;
  RooFIter iter2 = cache._denList.fwdIterator() ;
  RooFIter itern = cache._normList.fwdIterator() ;
  while((part=(RooAbsReal*)iterp.next())) {

    itern.next() ;
    num = (RooAbsReal*) iter1.next() ;
    den = (RooAbsReal*) iter2.next() ;

//     cout << "now processing part " << part->GetName() << " of type " << part->getStringAttribute("PROD_TERM_TYPE") << endl ;
//     cout << "corresponding numerator = " << num->GetName() << endl ;
//     cout << "corresponding denominator = " << den->GetName() << endl ;


    RooFormulaVar* ratio(0) ;
    RooArgSet origNumTerm ;

    if (string("SPECINT")==part->getStringAttribute("PROD_TERM_TYPE")) {

	RooRealIntegral* orig = (RooRealIntegral*) num;
	RooFormulaVar* specratio = (RooFormulaVar*) &orig->integrand() ;
	RooProduct* func = (RooProduct*) specratio->getParameter(0) ;

	RooArgSet* comps = orig->getComponents() ;
	RooFIter iter = comps->fwdIterator() ;
	RooAbsArg* carg ;
	while((carg=(RooAbsArg*)iter.next())) {
	  if (carg->getAttribute("RATIO_TERM")) {
	    ratio = (RooFormulaVar*)carg ;
	    break ;
	  }
	}
	delete comps ;

	if (ratio) {
	  RooCustomizer cust(*func,"blah") ;
	  cust.replaceArg(*ratio,RooFit::RooConst(1)) ;
	  RooAbsArg* funcCust = cust.build() ;
// 	  cout << "customized function = " << endl ;
// 	  funcCust->printComponentTree() ;
	  nomList.add(*funcCust) ;
	} else {	
	  nomList.add(*func) ;
	}


    } else {

      // Find the ratio term
      RooAbsReal* func = num;
      // If top level object is integral, navigate to integrand
      if (func->InheritsFrom(RooRealIntegral::Class())) {
	func = (RooAbsReal*) &((RooRealIntegral*)(func))->integrand();
      }
      if (func->InheritsFrom(RooProduct::Class())) {
// 	cout << "product term found: " ; func->Print() ;
	RooArgSet comps(((RooProduct*)(func))->components()) ;
	RooFIter iter = comps.fwdIterator() ;
	RooAbsArg* arg ;
	while((arg=(RooAbsArg*)iter.next())) {
	  if (arg->getAttribute("RATIO_TERM")) {
	    ratio = (RooFormulaVar*)(arg) ;
	  } else {
	    origNumTerm.add(*arg) ;
	  }
	}
      }

      if (ratio) {
// 	cout << "Found ratio term in numerator: " << ratio->GetName() << endl ;
// 	cout << "Adding only original term to numerator: " << origNumTerm << endl ;
	nomList.add(origNumTerm) ;
      } else {
	nomList.add(*num) ;
      }

    }

    for (list<string>::iterator iter = rangeComps.begin() ; iter != rangeComps.end() ; ++iter) {
      // If denominator is an integral, make a clone with the integration range adjusted to
      // the selected component of the normalization integral
//       cout << "NOW PROCESSING DENOMINATOR " << den->IsA()->GetName() << "::" << den->GetName() << endl ;

      if (string("SPECINT")==part->getStringAttribute("PROD_TERM_TYPE")) {

// 	cout << "create integral: SPECINT case" << endl ;
	RooRealIntegral* orig = (RooRealIntegral*) num;
	RooFormulaVar* specRatio = (RooFormulaVar*) &orig->integrand() ;
	specIntDeps.add(orig->intVars()) ;
	if (orig->intRange()) {
	  specIntRange = orig->intRange() ;
	}
	//RooProduct* numtmp = (RooProduct*) specRatio->getParameter(0) ;
	RooProduct* dentmp = (RooProduct*) specRatio->getParameter(1) ;

// 	cout << "numtmp = " << numtmp->IsA()->GetName() << "::" << numtmp->GetName() << endl ;
// 	cout << "dentmp = " << dentmp->IsA()->GetName() << "::" << dentmp->GetName() << endl ;

// 	cout << "denominator components are " << dentmp->components() << endl ;
	RooArgSet comps(dentmp->components()) ;
	RooFIter piter = comps.fwdIterator() ;
	RooAbsReal* parg ;
	while((parg=(RooAbsReal*)piter.next())) {
// 	  cout << "now processing denominator component " << parg->IsA()->GetName() << "::" << parg->GetName() << endl ;

	  if (ratio && parg->dependsOn(*ratio)) {
// 	    cout << "depends in value of ratio" << endl ;

	    // Make specialize ratio instance
	    RooAbsReal* specializedRatio = specializeRatio(*(RooFormulaVar*)ratio,iter->c_str()) ;
	
// 	    cout << "specRatio = " << endl ;
// 	    specializedRatio->printComponentTree() ;
	
	    // Replace generic ratio with specialized ratio
	    RooAbsArg *partCust(0) ;
	    if (parg->InheritsFrom(RooAddition::Class())) {
	


	      RooAddition* tmpadd = (RooAddition*)(parg) ;	

	      RooCustomizer cust(*tmpadd->list1().first(),Form("blah_%s",iter->c_str())) ;	
	      cust.replaceArg(*ratio,*specializedRatio) ;	
	      partCust = cust.build() ;
	
	    } else {
	      RooCustomizer cust(*parg,Form("blah_%s",iter->c_str())) ;	
	      cust.replaceArg(*ratio,*specializedRatio) ;	
	      partCust = cust.build() ;
	    }
	
	    // Print customized denominator
// 	    cout << "customized function = " << endl ;
// 	    partCust->printComponentTree() ;
	
	    RooAbsReal* specializedPartCust = specializeIntegral(*(RooAbsReal*)partCust,iter->c_str()) ;

	    // Finally divide again by ratio
	    string name = Form("%s_divided_by_ratio",specializedPartCust->GetName()) ;
	    RooFormulaVar* specIntFinal = new RooFormulaVar(name.c_str(),"@0/@1",RooArgList(*specializedPartCust,*specializedRatio)) ;
	
	    denListList[*iter].add(*specIntFinal) ;		
	  } else {

// 	    cout << "does NOT depend on value of ratio" << endl ;
// 	    parg->Print("t") ;

	    denListList[*iter].add(*specializeIntegral(*parg,iter->c_str())) ;		

	  }
	}
// 	cout << "end iteration over denominator components" << endl ;
      } else {

	if (ratio) {

	  RooAbsReal* specRatio = specializeRatio(*(RooFormulaVar*)ratio,iter->c_str()) ;

	  // If integral is 'Int r(y)*g(y) dy ' then divide a posteriori by r(y)
// 	  cout << "have ratio, orig den = " << den->GetName() << endl ;

	  RooArgSet tmp(origNumTerm) ;
	  tmp.add(*specRatio) ;
	  const string pname = makeRGPPName("PROD",tmp,RooArgSet(),RooArgSet(),0) ;
	  RooProduct* specDenProd = new RooProduct(pname.c_str(),pname.c_str(),tmp) ;
	  RooAbsReal* specInt(0) ;

	  if (den->InheritsFrom(RooRealIntegral::Class())) {
	    specInt = specDenProd->createIntegral(((RooRealIntegral*)den)->intVars(),iter->c_str()) ;
	  } else if (den->InheritsFrom(RooAddition::Class())) {	
	    RooAddition* orig = (RooAddition*)den ;
	    RooRealIntegral* origInt = (RooRealIntegral*) orig->list1().first() ;
	    specInt = specDenProd->createIntegral(origInt->intVars(),iter->c_str()) ;
	  } else {
	    throw string("this should not happen") ;
	  }
	
	  //RooAbsReal* specInt = specializeIntegral(*den,iter->c_str()) ;
	  string name = Form("%s_divided_by_ratio",specInt->GetName()) ;
	  RooFormulaVar* specIntFinal = new RooFormulaVar(name.c_str(),"@0/@1",RooArgList(*specInt,*specRatio)) ;
	  denListList[*iter].add(*specIntFinal) ;
	} else {
	  denListList[*iter].add(*specializeIntegral(*den,iter->c_str())) ;
	}
	
      }
    }

  }

  // Do not rearrage terms if numerator and denominator are effectively empty
  if (nomList.getSize()==0) {
    return ;
  }

  string name = Form("%s_numerator",GetName()) ;
  // WVE FIX THIS (2)

  RooAbsReal* numerator = new RooProduct(name.c_str(),name.c_str(),nomList) ;

  RooArgSet products ;
//   cout << "nomList = " << nomList << endl ;
  for (map<string,RooArgSet>::iterator iter = denListList.begin() ; iter != denListList.end() ; ++iter) {
//     cout << "denList[" << iter->first << "] = " << iter->second << endl ;
    name = Form("%s_denominator_comp_%s",GetName(),iter->first.c_str()) ;
    // WVE FIX THIS (2)
    RooProduct* prod_comp = new RooProduct(name.c_str(),name.c_str(),iter->second) ;
    products.add(*prod_comp) ;
  }
  name = Form("%s_denominator_sum",GetName()) ;
  RooAbsReal* norm = new RooAddition(name.c_str(),name.c_str(),products) ;
  norm->addOwnedComponents(products) ;

  if (specIntDeps.getSize()>0) {
    // Apply posterior integration required for SPECINT case

    string namesr = Form("SPEC_RATIO(%s,%s)",numerator->GetName(),norm->GetName()) ;
    RooFormulaVar* ndr = new RooFormulaVar(namesr.c_str(),"@0/@1",RooArgList(*numerator,*norm)) ;

    // Integral of ratio
    RooAbsReal* numtmp = ndr->createIntegral(specIntDeps,specIntRange.c_str()) ;

    numerator = numtmp ;
    norm = (RooAbsReal*) RooFit::RooConst(1).Clone() ;
  }


//   cout << "numerator" << endl ;
//   numerator->printComponentTree("",0,5) ;
//   cout << "denominator" << endl ;
//   norm->printComponentTree("",0,5) ;


  // WVE DEBUG
  //RooMsgService::instance().debugWorkspace()->import(RooArgSet(*numerator,*norm)) ;

  cache._rearrangedNum = numerator ;
  cache._rearrangedDen = norm ;
  cache._isRearranged = kTRUE ;

}


////////////////////////////////////////////////////////////////////////////////

RooAbsReal* RooProdPdf::specializeRatio(RooFormulaVar& input, const char* targetRangeName) const
{
  RooRealIntegral* numint = (RooRealIntegral*) input.getParameter(0) ;
  RooRealIntegral* denint = (RooRealIntegral*) input.getParameter(1) ;

  RooAbsReal* numint_spec = specializeIntegral(*numint,targetRangeName) ;

  RooAbsReal* ret =  new RooFormulaVar(Form("ratio(%s,%s)",numint_spec->GetName(),denint->GetName()),"@0/@1",RooArgList(*numint_spec,*denint)) ;
  ret->addOwnedComponents(*numint_spec) ;

  return ret ;
}



////////////////////////////////////////////////////////////////////////////////

RooAbsReal* RooProdPdf::specializeIntegral(RooAbsReal& input, const char* targetRangeName) const
{
  if (input.InheritsFrom(RooRealIntegral::Class())) {

    // If input is integral, recreate integral but override integration range to be targetRangeName
    RooRealIntegral* orig = (RooRealIntegral*)&input ;
//     cout << "creating integral: integrand =  " << orig->integrand().GetName() << " vars = " << orig->intVars() << " range = " << targetRangeName << endl ;
    return orig->integrand().createIntegral(orig->intVars(),targetRangeName) ;

  } else if (input.InheritsFrom(RooAddition::Class())) {

    // If input is sum of integrals, recreate integral from first component of set, but override integration range to be targetRangeName
    RooAddition* orig = (RooAddition*)&input ;
    RooRealIntegral* origInt = (RooRealIntegral*) orig->list1().first() ;
//     cout << "creating integral from addition: integrand =  " << origInt->integrand().GetName() << " vars = " << origInt->intVars() << " range = " << targetRangeName << endl ;
    return origInt->integrand().createIntegral(origInt->intVars(),targetRangeName) ;

  } else {

//     cout << "specializeIntegral: unknown input type " << input.IsA()->GetName() << "::" << input.GetName() << endl ;
  }

  return &input ;
}


////////////////////////////////////////////////////////////////////////////////
/// Group product into terms that can be calculated independently

void RooProdPdf::groupProductTerms(RooLinkedList& groupedTerms, RooArgSet& outerIntDeps, 
				   const RooLinkedList& terms, const RooLinkedList& norms, 
				   const RooLinkedList& imps, const RooLinkedList& ints, const RooLinkedList& /*cross*/) const
{
  // Start out with each term in its own group
  RooFIter tIter = terms.fwdIterator() ;
  RooArgSet* term ;
  while((term=(RooArgSet*)tIter.next())) {
    RooLinkedList* group = new RooLinkedList ;
    group->Add(term) ;
    groupedTerms.Add(group) ;
  }

  // Make list of imported dependents that occur in any term
  RooArgSet allImpDeps ;
  RooFIter iIter = imps.fwdIterator() ;
  RooArgSet *impDeps ;
  while((impDeps=(RooArgSet*)iIter.next())) {
    allImpDeps.add(*impDeps,kFALSE) ;
  }

  // Make list of integrated dependents that occur in any term
  RooArgSet allIntDeps ;
  iIter = ints.fwdIterator() ;
  RooArgSet *intDeps ;
  while((intDeps=(RooArgSet*)iIter.next())) {
    allIntDeps.add(*intDeps,kFALSE) ;
  }

  RooArgSet* tmp = (RooArgSet*) allIntDeps.selectCommon(allImpDeps) ;
  outerIntDeps.removeAll() ;
  outerIntDeps.add(*tmp) ;
  delete tmp ;

  // Now iteratively merge groups that should be (partially) integrated together
  RooFIter oidIter = outerIntDeps.fwdIterator() ;
  RooAbsArg* outerIntDep ;
  while ((outerIntDep =(RooAbsArg*)oidIter.next())) {

    // Collect groups that feature this dependent
    RooLinkedList* newGroup = 0 ;

    // Loop over groups
    RooLinkedList* group ;
    RooFIter glIter = groupedTerms.fwdIterator() ;
    Bool_t needMerge = kFALSE ;
    while((group=(RooLinkedList*)glIter.next())) {

      // See if any term in this group depends in any ay on outerDepInt
      RooArgSet* term2 ;
      RooFIter tIter2 = group->fwdIterator() ;
      while((term2=(RooArgSet*)tIter2.next())) {

	Int_t termIdx = terms.IndexOf(term2) ;
	RooArgSet* termNormDeps = (RooArgSet*) norms.At(termIdx) ;
	RooArgSet* termIntDeps = (RooArgSet*) ints.At(termIdx) ;
	RooArgSet* termImpDeps = (RooArgSet*) imps.At(termIdx) ;

	if (termNormDeps->contains(*outerIntDep) ||
	    termIntDeps->contains(*outerIntDep) ||
	    termImpDeps->contains(*outerIntDep)) {
	  needMerge = kTRUE ;
	}
	
      }

      if (needMerge) {
	// Create composite group if not yet existing
	if (newGroup==0) {
	  newGroup = new RooLinkedList ;
	}
	
	// Add terms of this group to new term
	tIter2 = group->fwdIterator() ;
	while((term2=(RooArgSet*)tIter2.next())) {
	  newGroup->Add(term2) ;	
	}

	// Remove this group from list and delete it (but not its contents)
	groupedTerms.Remove(group) ;
	delete group ;
      }
    }
    // If a new group has been created to merge terms dependent on current outerIntDep, add it to group list
    if (newGroup) {
      groupedTerms.Add(newGroup) ;
    }

  }
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate integrals of factorized product terms over observables iset while normalized
/// to observables in nset.

std::vector<RooAbsReal*> RooProdPdf::processProductTerm(const RooArgSet* nset, const RooArgSet* iset, const char* isetRangeName,
							const RooArgSet* term,const RooArgSet& termNSet, const RooArgSet& termISet,
							Bool_t& isOwned, Bool_t forceWrap) const
{
//    cout << "   FOLKERT::RooProdPdf(" << GetName() <<") processProductTerm nset = " << (nset?*nset:RooArgSet()) << endl 
//          << "   _normRange = " << _normRange << endl 
//          << "   iset = " << (iset?*iset:RooArgSet()) << endl 
//          << "   isetRangeName = " << (isetRangeName?isetRangeName:"<null>") << endl 
//          << "   term = " << (term?*term:RooArgSet()) << endl 
//          << "   termNSet = " << termNSet << endl 
//          << "   termISet = " << termISet << endl 
//          << "   isOwned = " << isOwned << endl 
//          << "   forceWrap = " << forceWrap << endl ;

  vector<RooAbsReal*> ret(3) ; ret[0] = 0 ; ret[1] = 0 ; ret[2] = 0 ;

  // CASE I: factorizing term: term is integrated over all normalizing observables
  // -----------------------------------------------------------------------------
  // Check if all observbales of this term are integrated. If so the term cancels
  if (termNSet.getSize()>0 && termNSet.getSize()==termISet.getSize() && isetRangeName==0) {


    //cout << "processProductTerm(" << GetName() << ") case I " << endl ;

    // Term factorizes
    return ret ;
  }

  // CASE II: Dropped terms: if term is entirely unnormalized, it should be dropped
  // ------------------------------------------------------------------------------
  if (nset && termNSet.getSize()==0) {

    //cout << "processProductTerm(" << GetName() << ") case II " << endl ;

    // Drop terms that are not asked to be normalized
    return ret ;
  }

  if (iset && termISet.getSize()>0) {
    if (term->getSize()==1) {

      // CASE IIIa: Normalized and partially integrated single PDF term
      //---------------------------------------------------------------

      RooAbsPdf* pdf = (RooAbsPdf*) term->first() ;

      RooAbsReal* partInt = pdf->createIntegral(termISet,termNSet,isetRangeName) ;
      //partInt->setOperMode(operMode()) ;
      partInt->setStringAttribute("PROD_TERM_TYPE","IIIa") ;

      isOwned=kTRUE ;

      //cout << "processProductTerm(" << GetName() << ") case IIIa func = " << partInt->GetName() << endl ;

      ret[0] = partInt ;

      // Split mode results
      ret[1] = pdf->createIntegral(termISet,isetRangeName) ;
      ret[2] = pdf->createIntegral(termNSet,normRange()) ;

      return ret ;


    } else {

      // CASE IIIb: Normalized and partially integrated composite PDF term
      //---------------------------------------------------------------

      // Use auxiliary class RooGenProdProj to calculate this term
      const std::string name = makeRGPPName("GENPROJ_",*term,termISet,termNSet,isetRangeName) ;
      RooAbsReal* partInt = new RooGenProdProj(name.c_str(),name.c_str(),*term,termISet,termNSet,isetRangeName) ;
      partInt->setStringAttribute("PROD_TERM_TYPE","IIIb") ;
      //partInt->setOperMode(operMode()) ;

      //cout << "processProductTerm(" << GetName() << ") case IIIb func = " << partInt->GetName() << endl ;

      isOwned=kTRUE ;
      ret[0] = partInt ;

      const std::string name1 = makeRGPPName("PROD",*term,RooArgSet(),RooArgSet(),0) ;

      // WVE FIX THIS
      RooProduct* tmp_prod = new RooProduct(name1.c_str(),name1.c_str(),*term) ;

      ret[1] = tmp_prod->createIntegral(termISet,isetRangeName) ;
      ret[2] = tmp_prod->createIntegral(termNSet,normRange()) ;

      return ret ;
    }
  }

  // CASE IVa: Normalized non-integrated composite PDF term
  // -------------------------------------------------------
  if (nset && nset->getSize()>0 && term->getSize()>1) {
    // Composite term needs normalized integration

    const std::string name = makeRGPPName("GENPROJ_",*term,termISet,termNSet,isetRangeName) ;
    RooAbsReal* partInt = new RooGenProdProj(name.c_str(),name.c_str(),*term,termISet,termNSet,isetRangeName,normRange()) ;
    partInt->setExpensiveObjectCache(expensiveObjectCache()) ;

    partInt->setStringAttribute("PROD_TERM_TYPE","IVa") ;
    //partInt->setOperMode(operMode()) ;

    //cout << "processProductTerm(" << GetName() << ") case IVa func = " << partInt->GetName() << endl ;

    isOwned=kTRUE ;
    ret[0] = partInt ;

    const std::string name1 = makeRGPPName("PROD",*term,RooArgSet(),RooArgSet(),0) ;

    // WVE FIX THIS
    RooProduct* tmp_prod = new RooProduct(name1.c_str(),name1.c_str(),*term) ;

    ret[1] = tmp_prod->createIntegral(termISet,isetRangeName) ;
    ret[2] = tmp_prod->createIntegral(termNSet,normRange()) ;

    return ret ;
  }

  // CASE IVb: Normalized, non-integrated single PDF term
  // -----------------------------------------------------
  RooFIter pIter = term->fwdIterator() ;
  RooAbsPdf* pdf ;
  while((pdf=(RooAbsPdf*)pIter.next())) {

    if (forceWrap) {

      // Construct representative name of normalization wrapper
      TString name(pdf->GetName()) ;
      name.Append("_NORM[") ;
      RooFIter nIter = termNSet.fwdIterator() ;
      RooAbsArg* arg ;
      Bool_t first(kTRUE) ;
      while((arg=(RooAbsArg*)nIter.next())) {
	if (!first) {
	  name.Append(",") ;
	} else {
	  first=kFALSE ;
	}		
	name.Append(arg->GetName()) ;
      }
      if (normRange()) {
	name.Append("|") ;
	name.Append(normRange()) ;
      }
      name.Append("]") ;

      RooAbsReal* partInt = new RooRealIntegral(name.Data(),name.Data(),*pdf,RooArgSet(),&termNSet) ;
      partInt->setStringAttribute("PROD_TERM_TYPE","IVb") ;
      isOwned=kTRUE ;

      //cout << "processProductTerm(" << GetName() << ") case IVb func = " << partInt->GetName() << endl ;

      ret[0] = partInt ;

      ret[1] = pdf->createIntegral(RooArgSet()) ;
      ret[2] = pdf->createIntegral(termNSet,normRange()) ;

      return ret ;


    } else {
      isOwned=kFALSE ;

      //cout << "processProductTerm(" << GetName() << ") case IVb func = " << pdf->GetName() << endl ;


      pdf->setStringAttribute("PROD_TERM_TYPE","IVb") ;
      ret[0] = pdf ;

      ret[1] = pdf->createIntegral(RooArgSet()) ;
      ret[2] = termNSet.getSize()>0 ? pdf->createIntegral(termNSet,normRange()) : ((RooAbsReal*)RooFit::RooConst(1).clone("1")) ;
      return ret  ;
    }
  }

  coutE(Eval) << "RooProdPdf::processProductTerm(" << GetName() << ") unidentified term!!!" << endl ;
  return ret ;
}




////////////////////////////////////////////////////////////////////////////////
/// Make an appropriate automatic name for a RooGenProdProj object in getPartIntList() 

std::string RooProdPdf::makeRGPPName(const char* pfx, const RooArgSet& term, const RooArgSet& iset,
				     const RooArgSet& nset, const char* isetRangeName) const
{
  // Make an appropriate automatic name for a RooGenProdProj object in getPartIntList()

  std::ostringstream os(pfx);
  os << "[";

  RooFIter pIter = term.fwdIterator() ;
  // Encode component names
  Bool_t first(kTRUE) ;
  RooAbsPdf* pdf ;
  while ((pdf=(RooAbsPdf*)pIter.next())) {
    if (!first) os << "_X_";
    first = kFALSE;
    os << pdf->GetName();
  }
  os << "]" << integralNameSuffix(iset,&nset,isetRangeName,kTRUE);

  return os.str();
}



////////////////////////////////////////////////////////////////////////////////
/// Force RooRealIntegral to offer all observables for internal integration

Bool_t RooProdPdf::forceAnalyticalInt(const RooAbsArg& /*dep*/) const 
{
  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Determine which part (if any) of given integral can be performed analytically.
/// If any analytical integration is possible, return integration scenario code.
///
/// RooProdPdf implements two strategies in implementing analytical integrals
///
/// First, PDF components whose entire set of dependents are requested to be integrated
/// can be dropped from the product, as they will integrate out to 1 by construction
///
/// Second, RooProdPdf queries each remaining component PDF for its analytical integration 
/// capability of the requested set ('allVars'). It finds the largest common set of variables 
/// that can be integrated by all remaining components. If such a set exists, it reconfirms that 
/// each component is capable of analytically integrating the common set, and combines the components 
/// individual integration codes into a single integration code valid for RooProdPdf.

Int_t RooProdPdf::getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, 
					  const RooArgSet* normSet, const char* rangeName) const 
{
  if (_forceNumInt) return 0 ;

  // Declare that we can analytically integrate all requested observables
  analVars.add(allVars) ;

  // Retrieve (or create) the required partial integral list
  Int_t code ;
  RooArgList *plist(0) ;
  RooLinkedList *nlist(0) ;
  getPartIntList(normSet,&allVars,plist,nlist,code,rangeName) ;

  return code+1 ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return analytical integral defined by given scenario code

Double_t RooProdPdf::analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName) const 
{
  // No integration scenario
  if (code==0) {
    return getVal(normSet) ;
  }


  // WVE needs adaptation for rangename feature

  // Partial integration scenarios
  CacheElem* cache = (CacheElem*) _cacheMgr.getObjByIndex(code-1) ;

  RooArgList* partIntList(0) ;
  RooLinkedList* normList(0) ;

  // If cache has been sterilized, revive this slot
  if (cache==0) {
    RooArgSet* vars = getParameters(RooArgSet()) ;
    RooArgSet* nset = _cacheMgr.nameSet1ByIndex(code-1)->select(*vars) ;
    RooArgSet* iset = _cacheMgr.nameSet2ByIndex(code-1)->select(*vars) ;

    Int_t code2(-1) ;
    getPartIntList(nset,iset,partIntList,normList,code2,rangeName) ;

    delete vars ;

    // preceding call to getPartIntList guarantees non-null return
    // coverity[NULL_RETURNS]
    cache = (CacheElem*) _cacheMgr.getObj(nset,iset,&code2,rangeName) ;

    delete nset ;
    delete iset ;

  } else {

    partIntList = &cache->_partList ;
    normList = &cache->_normList ;

  }

//   Double_t val = calculate(partIntList,normList) ;

  Double_t val = calculate(*cache,kTRUE) ;
//   cout << "RPP::aIWN(" << GetName() << ") ,code = " << code << ", value = " << val << endl ;

  return val ;
}



////////////////////////////////////////////////////////////////////////////////
/// Obsolete

Bool_t RooProdPdf::checkObservables(const RooArgSet* /*nset*/) const 
{ return kFALSE ; }




////////////////////////////////////////////////////////////////////////////////
/// If this product contains exactly one extendable p.d.f return the extension abilities of
/// that p.d.f, otherwise return CanNotBeExtended

RooAbsPdf::ExtendMode RooProdPdf::extendMode() const
{
  return (_extendedIndex>=0) ? ((RooAbsPdf*)_pdfList.at(_extendedIndex))->extendMode() : CanNotBeExtended ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return the expected number of events associated with the extentable input p.d.f
/// in the product. If there is no extendable term, return zero and issue and error

Double_t RooProdPdf::expectedEvents(const RooArgSet* nset) const 
{
  if (_extendedIndex<0) {
    coutE(Generation) << "ERROR: Requesting expected number of events from a RooProdPdf that does not contain an extended p.d.f" << endl ;
  }
  assert(_extendedIndex>=0) ;
  return ((RooAbsPdf*)_pdfList.at(_extendedIndex))->expectedEvents(nset) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return generator context optimized for generating events from product p.d.f.s

RooAbsGenContext* RooProdPdf::genContext(const RooArgSet &vars, const RooDataSet *prototype, 
					 const RooArgSet* auxProto, Bool_t verbose) const 
{
  if (_useDefaultGen) return RooAbsPdf::genContext(vars,prototype,auxProto,verbose) ;
  return new RooProdGenContext(*this,vars,prototype,auxProto,verbose) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Query internal generation capabilities of component p.d.f.s and aggregate capabilities
/// into master configuration passed to the generator context

Int_t RooProdPdf::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
{
  if (!_useDefaultGen) return 0 ;

  // Find the subset directVars that only depend on a single PDF in the product
  RooArgSet directSafe ;
  RooFIter dIter = directVars.fwdIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)dIter.next())) {
    if (isDirectGenSafe(*arg)) directSafe.add(*arg) ;
  }


  // Now find direct integrator for relevant components ;
  RooAbsPdf* pdf ;
  std::vector<Int_t> code;
  code.reserve(64);
  RooFIter pdfIter = _pdfList.fwdIterator();
  while((pdf=(RooAbsPdf*)pdfIter.next())) {
    RooArgSet pdfDirect ;
    Int_t pdfCode = pdf->getGenerator(directSafe,pdfDirect,staticInitOK);
    code.push_back(pdfCode);
    if (pdfCode != 0) {
      generateVars.add(pdfDirect) ;
    }
  }


  if (generateVars.getSize()>0) {
    Int_t masterCode = _genCode.store(code) ;
    return masterCode+1 ;
  } else {
    return 0 ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Forward one-time initialization call to component generation initialization
/// methods.

void RooProdPdf::initGenerator(Int_t code)
{
  if (!_useDefaultGen) return ;

  const std::vector<Int_t>& codeList = _genCode.retrieve(code-1) ;
  RooAbsPdf* pdf ;
  Int_t i(0) ;
  RooFIter pdfIter = _pdfList.fwdIterator();
  while((pdf=(RooAbsPdf*)pdfIter.next())) {
    if (codeList[i]!=0) {
      pdf->initGenerator(codeList[i]) ;
    }
    i++ ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Generate a single event with configuration specified by 'code'
/// Defer internal generation to components as encoded in the _genCode
/// registry for given generator code.

void RooProdPdf::generateEvent(Int_t code)
{
  if (!_useDefaultGen) return ;

  const std::vector<Int_t>& codeList = _genCode.retrieve(code-1) ;
  RooAbsPdf* pdf ;
  Int_t i(0) ;
  RooFIter pdfIter = _pdfList.fwdIterator();
  while((pdf=(RooAbsPdf*)pdfIter.next())) {
    if (codeList[i]!=0) {
      pdf->generateEvent(codeList[i]) ;
    }
    i++ ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooProdPdf::CacheElem::~CacheElem() 
{
  _normList.Delete() ; //WVE THIS IS AN INTENTIAL LEAK -- MUST FIX LATER
  if (_rearrangedNum) delete _rearrangedNum ;
  if (_rearrangedDen) delete _rearrangedDen ;
//   cout << "RooProdPdf::CacheElem dtor, this = " << this << endl ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return RooAbsArg components contained in the cache

RooArgList RooProdPdf::CacheElem::containedArgs(Action) 
{
  RooArgList ret ;
  ret.add(_partList) ;
  ret.add(_numList) ;
  ret.add(_denList) ;
  if (_rearrangedNum) ret.add(*_rearrangedNum) ;
  if (_rearrangedDen) ret.add(*_rearrangedDen) ;
  return ret ;

}



////////////////////////////////////////////////////////////////////////////////
/// Hook function to print cache contents in tree printing of RooProdPdf

void RooProdPdf::CacheElem::printCompactTreeHook(ostream& os, const char* indent, Int_t curElem, Int_t maxElem) 
{
   if (curElem==0) {
     os << indent << "RooProdPdf begin partial integral cache" << endl ;
   }

   RooFIter iter = _partList.fwdIterator() ;
   RooAbsArg* arg ;
   TString indent2(indent) ;
   indent2 += Form("[%d] ",curElem) ;
   while((arg=(RooAbsArg*)iter.next())) {
     arg->printCompactTree(os,indent2) ;
   }

   if (curElem==maxElem) {
     os << indent << "RooProdPdf end partial integral cache" << endl ;
   }
}



////////////////////////////////////////////////////////////////////////////////
/// Forward determination of safety of internal generator code to
/// component p.d.f that would generate the given observable

Bool_t RooProdPdf::isDirectGenSafe(const RooAbsArg& arg) const 
{
  // Only override base class behaviour if default generator method is enabled
  if (!_useDefaultGen) return RooAbsPdf::isDirectGenSafe(arg) ;

  // Argument may appear in only one PDF component
  RooAbsPdf* pdf, *thePdf(0) ;
  RooFIter pdfIter = _pdfList.fwdIterator();
  while((pdf=(RooAbsPdf*)pdfIter.next())) {

    if (pdf->dependsOn(arg)) {
      // Found PDF depending on arg

      // If multiple PDFs depend on arg directGen is not safe
      if (thePdf) return kFALSE ;

      thePdf = pdf ;
    }
  }
  // Forward call to relevant component PDF
  return thePdf?(thePdf->isDirectGenSafe(arg)):kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Look up user specified normalization set for given input PDF component

RooArgSet* RooProdPdf::findPdfNSet(RooAbsPdf& pdf) const 
{
  Int_t idx = _pdfList.index(&pdf) ;
  if (idx<0) return 0 ;
  return (RooArgSet*) _pdfNSetList.At(idx) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return all parameter constraint p.d.f.s on parameters listed in constrainedParams
/// The observables set is required to distinguish unambiguously p.d.f in terms 
/// of observables and parameters, which are not constraints, and p.d.fs in terms
/// of parameters only, which can serve as constraints p.d.f.s

RooArgSet* RooProdPdf::getConstraints(const RooArgSet& observables, RooArgSet& constrainedParams, Bool_t stripDisconnected) const
{
  RooArgSet constraints ;
  RooArgSet pdfParams, conParams ;

  // Loop over p.d.f. components
  RooFIter piter = _pdfList.fwdIterator() ;
  RooAbsPdf* pdf ;
  while((pdf=(RooAbsPdf*)piter.next())) {
    // A constraint term is a p.d.f that does not depend on any of the listed observables
    // but does depends on any of the parameters that should be constrained
    if (!pdf->dependsOnValue(observables) && pdf->dependsOnValue(constrainedParams)) {
      constraints.add(*pdf) ;
      RooArgSet* tmp = pdf->getParameters(observables) ;
      conParams.add(*tmp,kTRUE) ;
      delete tmp ;
    } else {
      RooArgSet* tmp = pdf->getParameters(observables) ;
      pdfParams.add(*tmp,kTRUE) ;
      delete tmp ;
    }
  }

  // Strip any constraints that are completely decoupled from the other product terms
  RooArgSet* finalConstraints = new RooArgSet("constraints") ;
  RooFIter citer = constraints.fwdIterator() ;
  while((pdf=(RooAbsPdf*)citer.next())) {
    if (pdf->dependsOnValue(pdfParams) || !stripDisconnected) {
      finalConstraints->add(*pdf) ;
    } else {
      coutI(Minimization) << "RooProdPdf::getConstraints(" << GetName() << ") omitting term " << pdf->GetName()
			  << " as constraint term as it does not share any parameters with the other pdfs in product. "
			  << "To force inclusion in likelihood, add an explicit Constrain() argument for the target parameter" << endl ;
    }
  }

  // Now remove from constrainedParams all parameters that occur exclusively in constraint term and not in regular pdf term

  RooArgSet* cexl = (RooArgSet*) conParams.selectCommon(constrainedParams) ;
  cexl->remove(pdfParams,kTRUE,kTRUE) ;
  constrainedParams.remove(*cexl,kTRUE,kTRUE) ;
  delete cexl ;

  return finalConstraints ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return all parameter constraint p.d.f.s on parameters listed in constrainedParams
/// The observables set is required to distinguish unambiguously p.d.f in terms 
/// of observables and parameters, which are not constraints, and p.d.fs in terms
/// of parameters only, which can serve as constraints p.d.f.s

RooArgSet* RooProdPdf::getConnectedParameters(const RooArgSet& observables) const
{
  RooFIter iter = _pdfList.fwdIterator() ;
  RooAbsArg* arg ;
  RooArgSet* connectedPars  = new RooArgSet("connectedPars") ;
  while((arg=iter.next())) {
    // Check if term is relevant
    if (arg->dependsOn(observables)) {
      RooArgSet* tmp = arg->getParameters(observables) ;
      connectedPars->add(*tmp) ;
      delete tmp ;
    } else {
    }
  }
  return connectedPars ;
}




////////////////////////////////////////////////////////////////////////////////

void RooProdPdf::getParametersHook(const RooArgSet* nset, RooArgSet* params, Bool_t stripDisconnected) const 
{
  if (!stripDisconnected) return ;
  if (!nset || nset->getSize()==0) return ;

  // Get/create appropriate term list for this normalization set
  RooArgList *plist(0) ;
  RooLinkedList *nlist(0) ;
  Int_t code ;
  getPartIntList(nset,0,plist,nlist,code) ;

  // Strip any terms from params that do not depend on any term
  RooFIter piter = params->fwdIterator() ;
  RooAbsReal* term, *param ;
  RooArgSet tostrip ;
  while((param=(RooAbsReal*)piter.next())) {
    Bool_t anyDep(kFALSE) ;
    RooFIter titer = plist->fwdIterator() ;
    while((term=(RooAbsReal*)titer.next())) {
      if (term->dependsOnValue(*param)) {
	anyDep=kTRUE ;
      }
    }
    if (!anyDep) {
      tostrip.add(*param) ;
    }
  }

  if (tostrip.getSize()>0) {
    params->remove(tostrip,kTRUE,kTRUE);
  }

}



////////////////////////////////////////////////////////////////////////////////
/// Interface function used by test statistics to freeze choice of range
/// for interpretation of conditional product terms

void RooProdPdf::selectNormalizationRange(const char* rangeName, Bool_t force) 
{
  if (!force && _refRangeName) {
    return ;
  }

  fixRefRange(rangeName) ;
}




////////////////////////////////////////////////////////////////////////////////

void RooProdPdf::fixRefRange(const char* rangeName)
{
  _refRangeName = (TNamed*)RooNameReg::ptr(rangeName) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Forward the plot sampling hint from the p.d.f. that defines the observable obs  

std::list<Double_t>* RooProdPdf::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  RooAbsPdf* pdf ;
  RooFIter pdfIter = _pdfList.fwdIterator();
  while((pdf=(RooAbsPdf*)pdfIter.next())) {
    list<Double_t>* hint = pdf->plotSamplingHint(obs,xlo,xhi) ;
    if (hint) {
      return hint ;
    }
  }

  return 0 ;
}



////////////////////////////////////////////////////////////////////////////////
/// If all components that depend on obs are binned that so is the product

Bool_t RooProdPdf::isBinnedDistribution(const RooArgSet& obs) const 
{
  RooAbsPdf* pdf ;
  RooFIter pdfIter = _pdfList.fwdIterator();
  while((pdf=(RooAbsPdf*)pdfIter.next())) {
    if (pdf->dependsOn(obs) && !pdf->isBinnedDistribution(obs)) {
      return kFALSE ;
    }
  }

  return kTRUE  ;
}






////////////////////////////////////////////////////////////////////////////////
/// Forward the plot sampling hint from the p.d.f. that defines the observable obs  

std::list<Double_t>* RooProdPdf::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  RooAbsPdf* pdf ;
  RooFIter pdfIter = _pdfList.fwdIterator();
  while((pdf=(RooAbsPdf*)pdfIter.next())) {
    list<Double_t>* hint = pdf->binBoundaries(obs,xlo,xhi) ;
    if (hint) {
      return hint ;
    }
  }

  return 0 ;
}


////////////////////////////////////////////////////////////////////////////////
/// Label OK'ed components of a RooProdPdf with cache-and-track, _and_ label all RooProdPdf
/// descendants with extra informationa about (conditional) normalization, needed to be able
/// to Cache-And-Track them outside the RooprodPdf context

void RooProdPdf::setCacheAndTrackHints(RooArgSet& trackNodes) 
{
  RooFIter piter = pdfList().fwdIterator() ;
  RooAbsArg* parg ;
  while ((parg=piter.next())) {

    if (parg->canNodeBeCached()==Always) {
      trackNodes.add(*parg) ;
      //cout << "tracking node RooProdPdf component " << parg->IsA()->GetName() << "::" << parg->GetName() << endl ;	

      // Additional processing to fix normalization sets in case product defines conditional observables
      RooArgSet* pdf_nset = findPdfNSet((RooAbsPdf&)(*parg)) ;
      if (pdf_nset) {
	// Check if conditional normalization is specified		
	if (string("nset")==pdf_nset->GetName() && pdf_nset->getSize()>0) {
	  RooNameSet n(*pdf_nset) ;
	  parg->setStringAttribute("CATNormSet",n.content()) ;
	}
	if (string("cset")==pdf_nset->GetName()) {
	  RooNameSet c(*pdf_nset) ;
	  parg->setStringAttribute("CATCondSet",c.content()) ;
	}
      } else {
	coutW(Optimization) << "RooProdPdf::setCacheAndTrackHints(" << GetName() << ") WARNING product pdf does not specify a normalization set for component " << parg->GetName() << endl ;
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Customized printing of arguments of a RooProdPdf to more intuitively reflect the contents of the
/// product operator construction

void RooProdPdf::printMetaArgs(ostream& os) const 
{
  RooFIter niter = _pdfNSetList.fwdIterator() ;
  for (int i=0 ; i<_pdfList.getSize() ; i++) {
    if (i>0) os << " * " ;
    RooArgSet* ncset = (RooArgSet*) niter.next() ;
    os << _pdfList.at(i)->GetName() ;
    if (ncset->getSize()>0) {
      if (string("nset")==ncset->GetName()) {
	os << *ncset  ;
      } else {
	os << "|" ;
	RooFIter nciter = ncset->fwdIterator() ;
	RooAbsArg* arg ;
	Bool_t first(kTRUE) ;
	while((arg=(RooAbsArg*)nciter.next())) {
	  if (!first) {
	    os << "," ;
	  } else {
	    first = kFALSE ;
	  }	
	  os << arg->GetName() ;	
	}
      }
    }
  }
  os << " " ;
}



////////////////////////////////////////////////////////////////////////////////
/// Implement support for node removal

Bool_t RooProdPdf::redirectServersHook(const RooAbsCollection& /*newServerList*/, Bool_t /*mustReplaceAll*/, Bool_t nameChange, Bool_t /*isRecursive*/) 
{
  if (nameChange && _pdfList.find("REMOVAL_DUMMY")) {

    cxcoutD(LinkStateMgmt) << "RooProdPdf::redirectServersHook(" << GetName() << "): removing REMOVAL_DUMMY" << endl ;

    // Remove node from _pdfList proxy and remove corresponding entry from normset list
    RooAbsArg* pdfDel = _pdfList.find("REMOVAL_DUMMY") ;

    TObject* setDel = _pdfNSetList.At(_pdfList.index("REMOVAL_DUMMY")) ;
    _pdfList.remove(*pdfDel) ;
    _pdfNSetList.Remove(setDel) ;

    // Clear caches
    _cacheMgr.reset() ;
  }
  return kFALSE ;
}
