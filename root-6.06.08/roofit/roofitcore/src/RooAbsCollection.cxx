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
// RooAbsCollection is an abstract container object that can hold
// multiple RooAbsArg objects.  Collections are ordered and can
// contain multiple objects of the same name, (but a derived
// implementation can enforce unique names). The storage of objects in
// implement through class RooLinkedList, a doubly linked list with an
// an optional hash-table lookup mechanism for fast indexing of large
// collections. 
// END_HTML
//
//

#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include "Riostream.h"
#include "TClass.h"
#include "TStopwatch.h"
#include "TRegexp.h"
#include "RooAbsCollection.h"
#include "RooStreamParser.h"
#include "RooFormula.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
#include "RooStringVar.h"
#include "RooTrace.h"
#include "RooArgList.h"
#include "RooLinkedListIter.h"
#include "RooCmdConfig.h"
#include "RooRealVar.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
using namespace std ;

#if (__GNUC__==3&&__GNUC_MINOR__==2&&__GNUC_PATCHLEVEL__==3)
char* operator+( streampos&, char* );
#endif

ClassImp(RooAbsCollection)
  ;

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooAbsCollection::RooAbsCollection() :
  _list(0),
  _ownCont(kFALSE), 
  _name(),
  _allRRV(kTRUE)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Empty collection constructor

RooAbsCollection::RooAbsCollection(const char *name) :
  _list(0),
  _ownCont(kFALSE), 
  _name(name),
  _allRRV(kTRUE)
{
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. Note that a copy of a collection is always non-owning,
/// even the source collection is owning. To create an owning copy of
/// a collection (owning or not), use the snaphot() method.

RooAbsCollection::RooAbsCollection(const RooAbsCollection& other, const char *name) :
  TObject(other),
  RooPrintable(other),
  _list(other._list.getHashTableSize()) , 
  _ownCont(kFALSE), 
  _name(name),
  _allRRV(other._allRRV)
{
  RooTrace::create(this) ;
  if (!name) setName(other.GetName()) ;
  
  // Transfer contents (not owned)
  RooFIter iterat= other.fwdIterator();
  RooAbsArg *arg = 0;
  while((arg= iterat.next())) {
    add(*arg);
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Destructor

RooAbsCollection::~RooAbsCollection() 
{
  // Delete all variables in our list if we own them
  if(_ownCont){ 
    safeDeleteList() ;
    //_list.Delete();
  }
}



////////////////////////////////////////////////////////////////////////////////

RooLinkedListIter RooAbsCollection::iterator(Bool_t dir) const 
{
    return _list.iterator(dir) ;
}

////////////////////////////////////////////////////////////////////////////////
/// Examine client server dependencies in list and
/// delete contents in safe order: any client
/// is deleted before a server is deleted

void RooAbsCollection::safeDeleteList() 
{
  // Handle trivial case here
  if (_list.GetSize() > 1) {
    std::vector<RooAbsArg*> tmp;
    tmp.reserve(_list.GetSize());
    do {
      tmp.clear();
      {
        RooFIter it = _list.fwdIterator();
        RooAbsArg* arg;
        while ((arg = it.next())) {
    	  // Check if arg depends on remainder of list      
	  if (!arg->dependsOn(*this, arg)) tmp.push_back(arg);
        }
      }
      // sort and uniquify, in case some elements occur more than once
      std::sort(tmp.begin(), tmp.end());
      // okay, can remove and delete what's in tmp
      for (std::vector<RooAbsArg*>::iterator it = tmp.begin(),
	  end = std::unique(tmp.begin(), tmp.end()); end != it; ++it) {
	while (_list.Remove(*it)) { };
	delete *it;
      }
    } while (!tmp.empty() && _list.GetSize() > 1);

    // Check if there are any remaining elements
    if (_list.GetSize() > 1) {
      coutW(ObjectHandling) << "RooAbsCollection::safeDeleteList(" << GetName() 
	<< ") WARNING: unable to delete following elements in client-server order " ;
      Print("1") ;
    }
  }

  // Built-in delete remaining elements
  _list.Delete() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Take a snap shot of current collection contents:
/// An owning collection is returned containing clones of 
/// 
///     - Elements in this collection 
///     - External dependents of all elements
///       and recursively any dependents of those dependents
///       (if deepCopy flag is set)
///
/// If deepCopy is specified, the client-server links between the cloned
/// list elements and the cloned external dependents are reconnected to
/// each other, making the snapshot a completely self-contained entity.
///
///

RooAbsCollection* RooAbsCollection::snapshot(Bool_t deepCopy) const
{
  // First create empty list
  TString snapName ;
  if (TString(GetName()).Length()>0) {
    snapName.Append("Snapshot of ") ;
    snapName.Append(GetName()) ;
  }
  RooAbsCollection* output = (RooAbsCollection*) create(snapName.Data()) ;
  if (deepCopy || getSize()>1000) {
    output->setHashTableSize(1000) ;
  }
  Bool_t error = snapshot(*output,deepCopy) ;
  if (error) {
    delete output ;
    return 0 ;
  }
  output->setHashTableSize(0) ;
  return output ;
}



////////////////////////////////////////////////////////////////////////////////
/// Take a snap shot of current collection contents:
/// An owning collection is returned containing clones of 
/// 
///     - Elements in this collection 
///     - External dependents of all elements
///       and recursively any dependents of those dependents
///       (if deepCopy flag is set)
///
/// If deepCopy is specified, the client-server links between the cloned
/// list elements and the cloned external dependents are reconnected to
/// each other, making the snapshot a completely self-contained entity.
///
///

Bool_t RooAbsCollection::snapshot(RooAbsCollection& output, Bool_t deepCopy) const 
{
  // Copy contents
  RooFIter iterat= fwdIterator();
  RooAbsArg *orig = 0;
  while((0 != (orig= iterat.next()))) {
    RooAbsArg *copy= (RooAbsArg*)orig->Clone();
    output.add(*copy);
  }

  RooFIter vIter = output.fwdIterator() ;
  RooAbsArg* var ;

  // Add external dependents
  Bool_t error(kFALSE) ;
  if (deepCopy) {
    // Recursively add clones of all servers
    while ((var=vIter.next())) {
      error |= output.addServerClonesToList(*var) ;
    }
  }

  // Handle eventual error conditions
  if (error) {
    coutE(ObjectHandling) << "RooAbsCollection::snapshot(): Errors occurred in deep clone process, snapshot not created" << endl ;
    output._ownCont = kTRUE ;    
    return kTRUE ;
  }

   // Redirect all server connections to internal list members
  vIter = output.fwdIterator() ;
  while ((var=vIter.next())) {
    var->redirectServers(output,deepCopy) ;
  }


  // Transfer ownership of contents to list
  output._ownCont = kTRUE ;
  return kFALSE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Add clones of servers of given argument to list

Bool_t RooAbsCollection::addServerClonesToList(const RooAbsArg& var)
{
  Bool_t ret(kFALSE) ;

  RooFIter sIter = var.serverMIterator() ;
  RooAbsArg* server ;
  while ((server=sIter.next())) {
    RooAbsArg* tmp = find(*server) ;
    if (!tmp) {
      RooAbsArg* serverClone = (RooAbsArg*)server->Clone() ;      
      serverClone->setAttribute("SnapShot_ExtRefClone") ;
      _list.Add(serverClone) ;      
      if (_allRRV && dynamic_cast<RooRealVar*>(serverClone)==0) {
	_allRRV=kFALSE ;
      }
      ret |= addServerClonesToList(*server) ;
    } else {
    }
  }
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// The assignment operator sets the value of any argument in our set
/// that also appears in the other set.

RooAbsCollection &RooAbsCollection::operator=(const RooAbsCollection& other) 
{
  if (&other==this) return *this ;

  RooAbsArg *elem, *theirs ;
  RooFIter iter = _list.fwdIterator() ;
  while((elem=iter.next())) {
    theirs= other.find(*elem);
    if(!theirs) continue;
    theirs->syncCache() ;
    elem->copyCache(theirs) ;
   elem->setAttribute("Constant",theirs->isConstant()) ;
  }
  return *this;
}



////////////////////////////////////////////////////////////////////////////////
/// The assignment operator sets the value of any argument in our set
/// that also appears in the other set.

RooAbsCollection &RooAbsCollection::assignValueOnly(const RooAbsCollection& other, Bool_t oneSafe) 
{
  if (&other==this) return *this ;
  
  // Short cut for 1 element assignment
  if (getSize()==1 && getSize()==other.getSize() && oneSafe) {
    other.first()->syncCache() ;
    first()->copyCache(other.first(),kTRUE) ;
    return *this ;
  }

  RooAbsArg *elem, *theirs ;
  RooFIter iter = _list.fwdIterator() ;
  while((elem=iter.next())) {
    theirs= other.find(*elem);
    if(!theirs) continue;
    theirs->syncCache() ;
    elem->copyCache(theirs,kTRUE) ;
  }
  return *this;
}



////////////////////////////////////////////////////////////////////////////////
/// Functional equivalent of operator=() but assumes this and other collection
/// have same layout. Also no attributes are copied

void RooAbsCollection::assignFast(const RooAbsCollection& other, Bool_t setValDirty) 
{
  if (&other==this) return ;

  RooFIter iter = _list.fwdIterator(), iter2 = other._list.fwdIterator() ;
  
  if (_allRRV) {
    
    RooRealVar *elem, *theirs ;
    // All contents are know to be RooRealVars - fast version of assignment
    while((elem=(RooRealVar*)iter.next())) {      
      // Identical size of iterators is documented assumption of method
      // coverity[NULL_RETURNS]
      theirs= (RooRealVar*)iter2.next() ;      
      elem->copyCacheFast(*theirs,setValDirty) ;
    }


  } else {
    
    RooAbsArg *elem, *theirs ;
    while((elem=iter.next())) {
      
      // Identical size of iterators is documented assumption of method
      // coverity[NULL_RETURNS]
      theirs= iter2.next() ;
      
      theirs->syncCache() ;
      elem->copyCache(theirs,kTRUE,setValDirty) ;
    }

  }
}



////////////////////////////////////////////////////////////////////////////////
/// Add the specified argument to list. Returns kTRUE if successful, or
/// else kFALSE if a variable of the same name is already in the list.
/// This method can only be called on a list that is flagged as owning
/// all of its contents, or else on an empty list (which will force the
/// list into that mode).

Bool_t RooAbsCollection::addOwned(RooAbsArg& var, Bool_t silent) 
{
  // check that we own our variables or else are empty
  if(!_ownCont && (getSize() > 0) && !silent) {
    coutE(ObjectHandling) << ClassName() << "::" << GetName() << "::addOwned: can only add to an owned list" << endl;
    return kFALSE;
  }
  _ownCont= kTRUE;

  _list.Add((RooAbsArg*)&var);
  if (_allRRV && dynamic_cast<RooRealVar*>(&var)==0) {
    _allRRV=kFALSE ;
  }

  return kTRUE;
}



////////////////////////////////////////////////////////////////////////////////
/// Add a clone of the specified argument to list. Returns a pointer to
/// the clone if successful, or else zero if a variable of the same name
/// is already in the list or the list does *not* own its variables (in
/// this case, try add() instead.) Calling addClone() on an empty list
/// forces it to take ownership of all its subsequent variables.

RooAbsArg *RooAbsCollection::addClone(const RooAbsArg& var, Bool_t silent) 
{
  // check that we own our variables or else are empty
  if(!_ownCont && (getSize() > 0) && !silent) {
    coutE(ObjectHandling) << ClassName() << "::" << GetName() << "::addClone: can only add to an owned list" << endl;
    return 0;
  }
  _ownCont= kTRUE;

  // add a pointer to a clone of this variable to our list (we now own it!)
  RooAbsArg *clone2= (RooAbsArg*)var.Clone();
  if(0 != clone2) _list.Add((RooAbsArg*)clone2);
  if (_allRRV && dynamic_cast<const RooRealVar*>(&var)==0) {
    _allRRV=kFALSE ;
  }

  return clone2;
}



////////////////////////////////////////////////////////////////////////////////
/// Add the specified argument to list. Returns kTRUE if successful, or
/// else kFALSE if a variable of the same name is already in the list
/// or the list owns its variables (in this case, try addClone() or addOwned() instead).

Bool_t RooAbsCollection::add(const RooAbsArg& var, Bool_t silent) 
{
  // check that this isn't a copy of a list
  if(_ownCont && !silent) {
    coutE(ObjectHandling) << ClassName() << "::" << GetName() << "::add: cannot add to an owned list" << endl;
    return kFALSE;
  }

  // add a pointer to this variable to our list (we don't own it!)
  _list.Add((RooAbsArg*)&var);
  if (_allRRV && dynamic_cast<const RooRealVar*>(&var)==0) {
    _allRRV=kFALSE ;
  }
  return kTRUE;
}



////////////////////////////////////////////////////////////////////////////////
/// Add a collection of arguments to this collection by calling add()
/// for each element in the source collection

Bool_t RooAbsCollection::add(const RooAbsCollection& list, Bool_t silent)
{
  Bool_t result(false) ;

  Int_t n= list.getSize() ;
  for(Int_t index= 0; index < n; index++) {
    result |= add((RooAbsArg&)*list._list.At(index),silent) ;
  }

  return result;  
}



////////////////////////////////////////////////////////////////////////////////
/// Add a collection of arguments to this collection by calling addOwned()
/// for each element in the source collection

Bool_t RooAbsCollection::addOwned(const RooAbsCollection& list, Bool_t silent)
{
  Bool_t result(false) ;

  Int_t n= list.getSize() ;
  for(Int_t index= 0; index < n; index++) {
    result |= addOwned((RooAbsArg&)*list._list.At(index),silent) ;
  }

  return result;  
}



////////////////////////////////////////////////////////////////////////////////
/// Add a collection of arguments to this collection by calling addOwned()
/// for each element in the source collection

void RooAbsCollection::addClone(const RooAbsCollection& list, Bool_t silent)
{
  Int_t n= list.getSize() ;
  for(Int_t index= 0; index < n; index++) {
    addClone((RooAbsArg&)*list._list.At(index),silent) ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Replace any args in our set with args of the same name from the other set
/// and return kTRUE for success. Fails if this list is a copy of another.

Bool_t RooAbsCollection::replace(const RooAbsCollection &other) 
{
  // check that this isn't a copy of a list
  if(_ownCont) {
    coutE(ObjectHandling) << "RooAbsCollection: cannot replace variables in a copied list" << endl;
    return kFALSE;
  }

  // loop over elements in the other list
  RooFIter otherArgs= other.fwdIterator();
  const RooAbsArg *arg = 0;
  while((arg= (const RooAbsArg*)otherArgs.next())) {

    // do we have an arg of the same name in our set?
    RooAbsArg *found= find(*arg);
    if(found) replace(*found,*arg);
  }
  return kTRUE;
}



////////////////////////////////////////////////////////////////////////////////
/// Replace var1 with var2 and return kTRUE for success. Fails if
/// this list is a copy of another, if var1 is not already in this set,
/// or if var2 is already in this set. var1 and var2 do not need to have
/// the same name.

Bool_t RooAbsCollection::replace(const RooAbsArg& var1, const RooAbsArg& var2) 
{
  // check that this isn't a copy of a list
  if(_ownCont) {
    coutE(ObjectHandling) << "RooAbsCollection: cannot replace variables in a copied list" << endl;
    return kFALSE;
  }

  // is var1 already in this list?
  const char *name= var1.GetName();

  if (!_list.FindObject(&var1)) {
    coutE(ObjectHandling) << "RooAbsCollection: variable \"" << name << "\" is not in the list"
	 << " and cannot be replaced" << endl;
    return kFALSE;
  }


  // is var2's name already in this list?
  if (dynamic_cast<RooArgSet*>(this)) {
    RooAbsArg *other = find(var2);
    if(other != 0 && other != &var1) {
      coutE(ObjectHandling) << "RooAbsCollection: cannot replace \"" << name
	   << "\" with already existing \"" << var2.GetName() << "\"" << endl;
      return kFALSE;
    }
  }

  // replace var1 with var2
  _list.Replace(&var1,&var2) ;
//   _list.AddBefore((RooAbsArg*)&var1,(RooAbsArg*)&var2);
//   _list.Remove((RooAbsArg*)&var1);

  if (_allRRV && dynamic_cast<const RooRealVar*>(&var2)==0) {
    _allRRV=kFALSE ;
  }

  return kTRUE;
}



////////////////////////////////////////////////////////////////////////////////
/// Remove the specified argument from our list. Return kFALSE if
/// the specified argument is not found in our list. An exact pointer
/// match is required, not just a match by name. A variable can be
/// removed from a copied list and will be deleted at the same time.

Bool_t RooAbsCollection::remove(const RooAbsArg& var, Bool_t , Bool_t matchByNameOnly) 
{
  // is var already in this list?
  TString name(var.GetName()) ;
  Bool_t anyFound(kFALSE) ;

  RooAbsArg* arg;
  while ((arg = (RooAbsArg*) _list.FindObject(&var))) {
    anyFound = kTRUE;
    _list.Remove(arg);
  }
  if (matchByNameOnly) {
      while ((arg = (RooAbsArg*) _list.FindObject(name.Data()))) {
	  anyFound = kTRUE;
	  _list.Remove(arg);
	  if (_ownCont) delete arg;
      }
  }
  
  return anyFound ;
}



////////////////////////////////////////////////////////////////////////////////
/// Remove each argument in the input list from our list using remove(const RooAbsArg&).
/// Return kFALSE in case of problems.

Bool_t RooAbsCollection::remove(const RooAbsCollection& list, Bool_t silent, Bool_t matchByNameOnly) 
{
  Bool_t result(false) ;

  Int_t n= list.getSize() ;
  for(Int_t index= 0; index < n; index++) {
    result |= remove((RooAbsArg&)*list._list.At(index),silent,matchByNameOnly) ;
  }

  return result;
}



////////////////////////////////////////////////////////////////////////////////
/// Remove all arguments from our set, deleting them if we own them.
/// This effectively restores our object to the state it would have
/// just after calling the RooAbsCollection(const char*) constructor.

void RooAbsCollection::removeAll() 
{
  if(_ownCont) {
    safeDeleteList() ;
    _ownCont= kFALSE;
  }
  else {
    _list.Clear();
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Set given attribute in each element of the collection by
/// calling each elements setAttribute() function.

void RooAbsCollection::setAttribAll(const Text_t* name, Bool_t value) 
{
  RooFIter iter= fwdIterator() ;
  RooAbsArg* arg ;
  while ((arg=iter.next())) {
    arg->setAttribute(name,value) ;
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Create a subset of the current collection, consisting only of those
/// elements with the specified attribute set. The caller is responsibe
/// for deleting the returned collection

RooAbsCollection* RooAbsCollection::selectByAttrib(const char* name, Bool_t value) const
{
  TString selName(GetName()) ;
  selName.Append("_selection") ;
  RooAbsCollection *sel = (RooAbsCollection*) create(selName.Data()) ;
  
  // Scan set contents for matching attribute
  RooFIter iter= fwdIterator() ;
  RooAbsArg* arg ;
  while ((arg=iter.next())) {
    if (arg->getAttribute(name)==value)
      sel->add(*arg) ;
  }

  return sel ;
}




////////////////////////////////////////////////////////////////////////////////
/// Create a subset of the current collection, consisting only of those
/// elements that are contained as well in the given reference collection.
/// The caller is responsible for deleting the returned collection

RooAbsCollection* RooAbsCollection::selectCommon(const RooAbsCollection& refColl) const 
{
  // Create output set
  TString selName(GetName()) ;
  selName.Append("_selection") ;
  RooAbsCollection *sel = (RooAbsCollection*) create(selName.Data()) ; 

  // Scan set contents for matching attribute
  RooFIter iter= fwdIterator() ;
  RooAbsArg* arg ;
  while ((arg=iter.next())) {
    if (refColl.find(*arg))
      sel->add(*arg) ;
  }

  return sel ;
}



////////////////////////////////////////////////////////////////////////////////
/// Create a subset of the current collection, consisting only of those
/// elements with names matching the wildcard expressions in nameList,
/// supplied as a comma separated list

RooAbsCollection* RooAbsCollection::selectByName(const char* nameList, Bool_t verbose) const 
{
  // Create output set
  TString selName(GetName()) ;
  selName.Append("_selection") ;
  RooAbsCollection *sel = (RooAbsCollection*) create(selName.Data()) ; 

  const size_t bufSize = strlen(nameList) + 1;
  char* buf = new char[bufSize] ;
  strlcpy(buf,nameList,bufSize) ;
  char* wcExpr = strtok(buf,",") ;
  while(wcExpr) {
    TRegexp rexp(wcExpr,kTRUE) ;
    if (verbose) {
      cxcoutD(ObjectHandling) << "RooAbsCollection::selectByName(" << GetName() << ") processing expression '" << wcExpr << "'" << endl ;
    }

    RooFIter iter = fwdIterator() ;
    RooAbsArg* arg ;
    while((arg=iter.next())) {
      if (TString(arg->GetName()).Index(rexp)>=0) {
	if (verbose) {
	  cxcoutD(ObjectHandling) << "RooAbsCollection::selectByName(" << GetName() << ") selected element " << arg->GetName() << endl ;
	}
	sel->add(*arg) ;
      }
    }
    wcExpr = strtok(0,",") ;
  }
  delete[] buf ;

  return sel ;
}




////////////////////////////////////////////////////////////////////////////////
/// Check if this and other collection have identically named contents

Bool_t RooAbsCollection::equals(const RooAbsCollection& otherColl) const
{
  // First check equal length 
  if (getSize() != otherColl.getSize()) return kFALSE ;

  // Then check that each element of our list also occurs in the other list
  RooFIter iter = fwdIterator() ;
  RooAbsArg* arg ;
  while((arg=iter.next())) {
    if (!otherColl.find(*arg)) {
      return kFALSE ;
    }
  }
  return kTRUE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Check if this and other collection have common entries

Bool_t RooAbsCollection::overlaps(const RooAbsCollection& otherColl) const 
{
  RooFIter iter = fwdIterator() ;
  RooAbsArg* arg ;
  while((arg=iter.next())) {
    if (otherColl.find(*arg)) {
      return kTRUE ;
    }
  }
  return kFALSE ;
}




////////////////////////////////////////////////////////////////////////////////
/// Find object with given name in list. A null pointer 
/// is returned if no object with the given name is found

RooAbsArg *RooAbsCollection::find(const char *name) const 
{
  return (RooAbsArg*) _list.find(name);
}



////////////////////////////////////////////////////////////////////////////////
/// Find object with given name in list. A null pointer 
/// is returned if no object with the given name is found

RooAbsArg *RooAbsCollection::find(const RooAbsArg& arg) const 
{
  return (RooAbsArg*) _list.findArg(&arg);
}



////////////////////////////////////////////////////////////////////////////////
/// Return comma separated list of contained object names as STL string

string RooAbsCollection::contentsString() const 
{
  string retVal ;
  RooFIter iter = fwdIterator() ;
  RooAbsArg* arg ;
  Bool_t isFirst(kTRUE) ;
  while((arg=iter.next())) {
    if (isFirst) {
      isFirst=kFALSE ;
    } else {
      retVal += "," ;
    }
    retVal += arg->GetName() ;
  }
  return retVal ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return collection name

void RooAbsCollection::printName(ostream& os) const 
{
  os << GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return collection title

void RooAbsCollection::printTitle(ostream& os) const 
{
  os << GetTitle() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return collection class name

void RooAbsCollection::printClassName(ostream& os) const 
{
  os << IsA()->GetName() ;
}



////////////////////////////////////////////////////////////////////////////////
/// Define default RooPrinable print options for given Print() flag string
/// For inline printing only show value of objects, for default print show
/// name,class name value and extras of each object. In verbose mode
/// also add object adress, argument and title

Int_t RooAbsCollection::defaultPrintContents(Option_t* opt) const 
{
  if (opt && TString(opt)=="I") {
    return kValue ;
  }
  if (opt && TString(opt).Contains("v")) {
    return kAddress|kName|kArgs|kClassName|kValue|kTitle|kExtras ;
  }
  return kName|kClassName|kValue ;
}





////////////////////////////////////////////////////////////////////////////////
/// Print value of collection, i.e. a comma separated list of contained
/// object names

void RooAbsCollection::printValue(ostream& os) const
{
  Bool_t first2(kTRUE) ;
  os << "(" ;
  RooFIter iter = fwdIterator() ;
  RooAbsArg* arg ;
  while((arg=iter.next())) {
    if (!first2) {
      os << "," ;
    } else {
      first2 = kFALSE ;
    }
    os << arg->GetName() ;
    
  }
  os << ")" ;  
}



////////////////////////////////////////////////////////////////////////////////
/// Implement multiline printin of collection, one line for each ontained object showing
/// the requested content

void RooAbsCollection::printMultiline(ostream&os, Int_t contents, Bool_t /*verbose*/, TString indent) const
{
  if (TString(GetName()).Length()>0 && (contents&kCollectionHeader)) {
    os << indent << ClassName() << "::" << GetName() << ":" << (_ownCont?" (Owning contents)":"") << endl;
  }

  RooFIter iterat= fwdIterator();
  int index= 0;
  RooAbsArg *next = 0;
  TString deeper(indent);
  deeper.Append("     ");
  
  // Adjust the with of the name field to fit the largest name, if requesed
  Int_t maxNameLen(1) ;
  Int_t nameFieldLengthSaved = RooPrintable::_nameLength ;
  if (nameFieldLengthSaved==0) {
    while((next=iterat.next())) {
      Int_t len = strlen(next->GetName()) ;
      if (len>maxNameLen) maxNameLen = len ;
    }
    iterat = fwdIterator() ;
    RooPrintable::nameFieldLength(maxNameLen+1) ;
  }
  
  while((0 != (next= iterat.next()))) {
    os << indent << setw(3) << ++index << ") ";
    next->printStream(os,contents,kSingleLine,"");
  }
  
  // Reset name field length, if modified
  RooPrintable::nameFieldLength(nameFieldLengthSaved) ;
}



////////////////////////////////////////////////////////////////////////////////
/// Base contents dumper for debugging purposes

void RooAbsCollection::dump() const 
{
  RooFIter iter = fwdIterator() ;
  RooAbsArg* arg ;
  while((arg=iter.next())) {
    cout << arg << " " << arg->IsA()->GetName() << "::" << arg->GetName() << " (" << arg->GetTitle() << ")" << endl ;
  }
}



////////////////////////////////////////////////////////////////////////////////
/// Output content of collection as LaTex table. By default a table with two columns is created: the left
/// column contains the name of each variable, the right column the value.
///
/// The following optional named arguments can be used to modify the default behavior
///
///   Columns(Int_t ncol)                    -- Fold table into multiple columns, i.e. ncol=3 will result in 3 x 2 = 6 total columns
///   Sibling(const RooAbsCollection& other) -- Define sibling list. The sibling list is assumed to have objects with the same
///                                             name in the same order. If this is not the case warnings will be printed. If a single
///                                             sibling list is specified, 3 columns will be output: the (common) name, the value of this
///                                             list and the value in the sibling list. Multiple sibling lists can be specified by 
///                                             repeating the Sibling() command. 
///   Format(const char* str)                -- Classic format string, provided for backward compatibility
///   Format(...)                            -- Formatting arguments, details are given below
///   OutputFile(const char* fname)          -- Send output to file with given name rather than standard output
///
/// The Format(const char* what,...) has the following structure
///
///   const char* what          -- Controls what is shown. "N" adds name, "E" adds error, 
///                                "A" shows asymmetric error, "U" shows unit, "H" hides the value
///   FixedPrecision(int n)     -- Controls precision, set fixed number of digits
///   AutoPrecision(int n)      -- Controls precision. Number of shown digits is calculated from error 
///                                + n specified additional digits (1 is sensible default)
///   VerbatimName(Bool_t flag) -- Put variable name in a \verb+   + clause.
///
/// Example use: list.printLatex(Columns(2), Format("NEU",AutoPrecision(1),VerbatimName()) ) ;

void RooAbsCollection::printLatex(const RooCmdArg& arg1, const RooCmdArg& arg2,
				  const RooCmdArg& arg3, const RooCmdArg& arg4,	
				  const RooCmdArg& arg5, const RooCmdArg& arg6,	
				  const RooCmdArg& arg7, const RooCmdArg& arg8) const
{

  
  // Define configuration for this method
  RooCmdConfig pc("RooAbsCollection::printLatex()") ;
  pc.defineInt("ncol","Columns",0,1) ;
  pc.defineString("outputFile","OutputFile",0,"") ;
  pc.defineString("format","Format",0,"NEYVU") ;
  pc.defineInt("sigDigit","Format",0,1) ;
  pc.defineObject("siblings","Sibling",0,0,kTRUE) ;
  pc.defineInt("dummy","FormatArgs",0,0) ;
  pc.defineMutex("Format","FormatArgs") ;
 
  // Stuff all arguments in a list
  RooLinkedList cmdList;
  cmdList.Add(const_cast<RooCmdArg*>(&arg1)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg2)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg3)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg4)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg5)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg6)) ;
  cmdList.Add(const_cast<RooCmdArg*>(&arg7)) ;  cmdList.Add(const_cast<RooCmdArg*>(&arg8)) ;

  // Process & check varargs 
  pc.process(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8) ;
  if (!pc.ok(kTRUE)) {
    return ;
  }

  const char* outFile = pc.getString("outputFile") ;
  if (outFile && strlen(outFile)) {
    ofstream ofs(outFile) ;
    if (pc.hasProcessed("FormatArgs")) {
      RooCmdArg* formatCmd = static_cast<RooCmdArg*>(cmdList.FindObject("FormatArgs")) ;
      formatCmd->addArg(RooFit::LatexTableStyle()) ;
      printLatex(ofs,pc.getInt("ncol"),0,0,pc.getObjectList("siblings"),formatCmd) ;    
    } else {
      printLatex(ofs,pc.getInt("ncol"),pc.getString("format"),pc.getInt("sigDigit"),pc.getObjectList("siblings")) ;
    }
  } else {
    if (pc.hasProcessed("FormatArgs")) {
      RooCmdArg* formatCmd = static_cast<RooCmdArg*>(cmdList.FindObject("FormatArgs")) ;
      formatCmd->addArg(RooFit::LatexTableStyle()) ;
      printLatex(cout,pc.getInt("ncol"),0,0,pc.getObjectList("siblings"),formatCmd) ;    
    } else {
      printLatex(cout,pc.getInt("ncol"),pc.getString("format"),pc.getInt("sigDigit"),pc.getObjectList("siblings")) ;
    }
  }
}




////////////////////////////////////////////////////////////////////////////////
/// Internal implementation function of printLatex

void RooAbsCollection::printLatex(ostream& ofs, Int_t ncol, const char* option, Int_t sigDigit, const RooLinkedList& siblingList, const RooCmdArg* formatCmd) const 
{
  // Count number of rows to print
  Int_t nrow = (Int_t) (getSize() / ncol + 0.99) ;
  Int_t i,j,k ;

  // Sibling list do not need to print their name as it is supposed to be the same  
  TString sibOption ;
  RooCmdArg sibFormatCmd ;
  if (option) {
    sibOption = option ;
    sibOption.ReplaceAll("N","") ;
    sibOption.ReplaceAll("n","") ;
  } else {
    sibFormatCmd = *formatCmd ;
    TString tmp = formatCmd->_s[0] ;
    tmp.ReplaceAll("N","") ;    
    tmp.ReplaceAll("n","") ;    
    static char buf[100] ;
    strlcpy(buf,tmp.Data(),100) ;
    sibFormatCmd._s[0] = buf ;
  }


  // Make list of lists ;
  RooLinkedList listList ;
  listList.Add((RooAbsArg*)this) ;
  RooFIter sIter = siblingList.fwdIterator() ;
  RooAbsCollection* col ;
  while((col=(RooAbsCollection*)sIter.next())) {
    listList.Add(col) ;
  }

  RooLinkedList listListRRV ;

  // Make list of RRV-only components
  RooFIter lIter = listList.fwdIterator() ;
  RooArgList* prevList = 0 ;
  while((col=(RooAbsCollection*)lIter.next())) {
    RooArgList* list = new RooArgList ;
    RooFIter iter = col->fwdIterator() ;
    RooAbsArg* arg ;
    while((arg=iter.next())) {    
      
      RooRealVar* rrv = dynamic_cast<RooRealVar*>(arg) ;
      if (rrv) {
	list->add(*rrv) ;
      } else {
	coutW(InputArguments) << "RooAbsCollection::printLatex: can only print RooRealVar in LateX, skipping non-RooRealVar object named "
	     << arg->GetName() << endl ;      
      }
      if (prevList && TString(rrv->GetName()).CompareTo(prevList->at(list->getSize()-1)->GetName())) {
	coutW(InputArguments) << "RooAbsCollection::printLatex: WARNING: naming and/or ordering of sibling list is different" << endl ;
      }
    }
    listListRRV.Add(list) ;
    if (prevList && list->getSize() != prevList->getSize()) {
      coutW(InputArguments) << "RooAbsCollection::printLatex: ERROR: sibling list(s) must have same length as self" << endl ;
      delete list ;
      listListRRV.Delete() ;
      return ;
    }
    prevList = list ;
  }

  // Construct table header
  Int_t nlist = listListRRV.GetSize() ;
  TString subheader = "l" ;
  for (k=0 ; k<nlist ; k++) subheader += "c" ;

  TString header = "\\begin{tabular}{" ;
  for (j=0 ; j<ncol ; j++) {
    if (j>0) header += "|" ;
    header += subheader ;
  }
  header += "}" ;
  ofs << header << endl ;


  // Print contents, delegating actual printing to RooRealVar::format()
  for (i=0 ; i<nrow ; i++) {
    for (j=0 ; j<ncol ; j++) {
      for (k=0 ; k<nlist ; k++) {
	RooRealVar* par = (RooRealVar*) ((RooArgList*)listListRRV.At(k))->at(i+j*nrow) ;
	if (par) {
	  if (option) {
	    TString* tmp = par->format(sigDigit,(k==0)?option:sibOption.Data()) ;
	    ofs << *tmp ;
	    delete tmp ;
	  } else {
	    TString* tmp = par->format((k==0)?*formatCmd:sibFormatCmd) ;
	    ofs << *tmp ;
	    delete tmp ;
	  }
	}
	if (!(j==ncol-1 && k==nlist-1)) {
	  ofs << " & " ;
	}
      }
    }
    ofs << "\\\\" << endl ;
  }
  
  ofs << "\\end{tabular}" << endl ;
  listListRRV.Delete() ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return true if all contained object report to have their
/// value inside the specified range

Bool_t RooAbsCollection::allInRange(const char* rangeSpec) const
{
  if (!rangeSpec) return kTRUE ;

  // Parse rangeSpec specification
  vector<string> cutVec ;
  if (rangeSpec && strlen(rangeSpec)>0) {
    if (strchr(rangeSpec,',')==0) {
      cutVec.push_back(rangeSpec) ;
    } else {
      const size_t bufSize = strlen(rangeSpec)+1;
      char* buf = new char[bufSize] ;
      strlcpy(buf,rangeSpec,bufSize) ;
      const char* oneRange = strtok(buf,",") ;
      while(oneRange) {
	cutVec.push_back(oneRange) ;
	oneRange = strtok(0,",") ;
      }
      delete[] buf ;
    }
  }


  RooFIter iter = _list.fwdIterator() ;

  // Apply range based selection criteria
  Bool_t selectByRange = kTRUE ;
  RooAbsArg* arg ;
  while((arg=iter.next())) {
    Bool_t selectThisArg = kFALSE ;
    UInt_t icut ;
    for (icut=0 ; icut<cutVec.size() ; icut++) {
      if (arg->inRange(cutVec[icut].c_str())) {
	selectThisArg = kTRUE ;
	break ;
      }
    }
    if (!selectThisArg) {
      selectByRange = kFALSE ;
      break ;
    }
  }

  return selectByRange ;
}



////////////////////////////////////////////////////////////////////////////////

void RooAbsCollection::makeStructureTag() 
{
}


////////////////////////////////////////////////////////////////////////////////

void RooAbsCollection::makeTypedStructureTag() 
{
}





