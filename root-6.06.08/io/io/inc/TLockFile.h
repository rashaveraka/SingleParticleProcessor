// @(#)root/io:$Id$
// Author: Jan Fiete Grosse-Oetringhaus, 04.06.07

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TLockFile
#define ROOT_TLockFile

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif


class TLockFile : public TObject {

private:
   TLockFile(const TLockFile&);             // not implemented
   TLockFile& operator=(const TLockFile&);  // not implemented

protected:
   TString fPath;         ///< Path to file holding the lock

   Bool_t Lock(const char *path, Int_t timeLimit);

public:
   TLockFile(const char *path, Int_t timeLimit = 0);
   virtual ~TLockFile();

   ClassDef(TLockFile, 0) //Lock an object using a file
};

#endif
