// @(#)root/tree:$Id$
// Author: Rene Brun   05/07/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TVirtualIndex
Abstract interface for Tree Index
*/

#include "TVirtualIndex.h"
#include "TTree.h"

ClassImp(TVirtualIndex)

////////////////////////////////////////////////////////////////////////////////
/// Default constructor for TVirtualIndex

TVirtualIndex::TVirtualIndex(): TNamed()
{
   fTree         = 0;
}

////////////////////////////////////////////////////////////////////////////////

TVirtualIndex::~TVirtualIndex()
{
}
