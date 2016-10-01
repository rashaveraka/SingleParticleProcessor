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
// RooGenericPdf is a concrete implementation of a probability density function,
// which takes a RooArgList of servers and a C++ expression string defining how
// its value should be calculated from the given list of servers.
// A fully numerical integration is automatically performed to normalize the given
// expression. RooGenericPdf uses a RooFormula object to perform the expression evaluation
//
// The string expression can be any valid TFormula expression referring to the
// listed servers either by name or by their ordinal list position:
//
//   RooGenericPdf("gen","x*y",RooArgList(x,y))  or
//   RooGenericPdf("gen","@0*@1",RooArgList(x,y)) 
//
// The latter form, while slightly less readable, is more versatile because it
// doesn't hardcode any of the variable names it expects
// END_HTML
//

#include "RooFit.h"
#include "Riostream.h"

#include "RooGenericPdf.h"
#include "RooGenericPdf.h"
#include "RooStreamParser.h"
#include "RooMsgService.h"
#include "RooArgList.h"



using namespace std;

ClassImp(RooGenericPdf)



////////////////////////////////////////////////////////////////////////////////
/// Constructor with formula expression and list of input variables

RooGenericPdf::RooGenericPdf(const char *name, const char *title, const RooArgList& dependents) : 
  RooAbsPdf(name,title), 
  _actualVars("actualVars","Variables used by PDF expression",this),
  _formula(0),
  _formExpr(title)
{  
  _actualVars.add(dependents) ; 

  if (_actualVars.getSize()==0) _value = traceEval(0) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor with a name, title, formula expression and a list of variables

RooGenericPdf::RooGenericPdf(const char *name, const char *title, 
			     const char* inFormula, const RooArgList& dependents) : 
  RooAbsPdf(name,title), 
  _actualVars("actualVars","Variables used by PDF expression",this),
  _formula(0),
  _formExpr(inFormula)
{  
  _actualVars.add(dependents) ; 

  if (_actualVars.getSize()==0) _value = traceEval(0) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooGenericPdf::RooGenericPdf(const RooGenericPdf& other, const char* name) : 
  RooAbsPdf(other, name), 
  _actualVars("actualVars",this,other._actualVars),
  _formula(0),
  _formExpr(other._formExpr)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooGenericPdf::~RooGenericPdf() 
{
  if (_formula) delete _formula ;
}



////////////////////////////////////////////////////////////////////////////////

RooFormula& RooGenericPdf::formula() const
{
  if (!_formula) {
    _formula = new RooFormula(GetName(),_formExpr.Data(),_actualVars) ;
  } 
  return *_formula ;
}



////////////////////////////////////////////////////////////////////////////////
/// Calculate current value of this object

Double_t RooGenericPdf::evaluate() const
{
  return formula().eval(_normSet) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Change formula expression to given expression

Bool_t RooGenericPdf::setFormula(const char* inFormula) 
{
  if (formula().reCompile(inFormula)) return kTRUE ;

  _formExpr = inFormula ;
  setValueDirty() ;
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Check if given value is valid

Bool_t RooGenericPdf::isValidReal(Double_t /*value*/, Bool_t /*printError*/) const 
{
  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Propagate server changes to embedded formula object

Bool_t RooGenericPdf::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t /*isRecursive*/)
{
  if (_formula) {
     return _formula->changeDependents(newServerList,mustReplaceAll,nameChange) ;
  } else {
    return kTRUE ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Print info about this object to the specified stream. 

void RooGenericPdf::printMultiline(ostream& os, Int_t content, Bool_t verbose, TString indent) const
{
  RooAbsPdf::printMultiline(os,content,verbose,indent);
  if (verbose) {
    os << " --- RooGenericPdf --- " << endl ;
    indent.Append("  ");
    os << indent ;
    formula().printMultiline(os,content,verbose,indent);
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Add formula expression as meta argument in printing interface

void RooGenericPdf::printMetaArgs(ostream& os) const 
{
  os << "formula=\"" << _formExpr << "\" " ;
}



////////////////////////////////////////////////////////////////////////////////
/// Read object contents from given stream

Bool_t RooGenericPdf::readFromStream(istream& is, Bool_t compact, Bool_t /*verbose*/)
{
  if (compact) {
    coutE(InputArguments) << "RooGenericPdf::readFromStream(" << GetName() << "): can't read in compact mode" << endl ;
    return kTRUE ;
  } else {
    RooStreamParser parser(is) ;
    return setFormula(parser.readLine()) ;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Write object contents to given stream

void RooGenericPdf::writeToStream(ostream& os, Bool_t compact) const
{
  if (compact) {
    os << getVal() << endl ;
  } else {
    os << GetTitle() ;
  }
}



