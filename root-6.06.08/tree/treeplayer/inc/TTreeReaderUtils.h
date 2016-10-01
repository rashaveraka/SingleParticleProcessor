// @(#)root/tree:$Id$
// Author: Axel Naumann, 2010-10-12

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeReaderUtils
#define ROOT_TTreeReaderUtils


////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TTreeReaderUtils                                                       //
//                                                                        //
// TTreeReader's helpers.                                                 //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TBranchProxyDirector
#include "TBranchProxyDirector.h"
#endif
#ifndef ROOT_TBranchProxy
#include "TBranchProxy.h"
#endif
#include "TTreeReaderValue.h"

class TDictionary;
class TTree;

namespace ROOT {
   namespace Detail {
      class TBranchProxy;
   }

namespace Internal {
   class TBranchProxyDirector;
   class TTreeReaderArrayBase;

   class TNamedBranchProxy: public TObject {
   public:
      TNamedBranchProxy(): fDict(0), fContentDict(0) {}
      TNamedBranchProxy(TBranchProxyDirector* boss, TBranch* branch, const char* membername):
         fProxy(boss, branch, membername), fDict(0), fContentDict(0) {}

      const char* GetName() const { return fProxy.GetBranchName(); }
      const Detail::TBranchProxy* GetProxy() const { return &fProxy; }
      Detail::TBranchProxy* GetProxy() { return &fProxy; }
      TDictionary* GetDict() const { return fDict; }
      void SetDict(TDictionary* dict) { fDict = dict; }
      TDictionary* GetContentDict() const { return fContentDict; }
      void SetContentDict(TDictionary* dict) { fContentDict = dict; }

   private:
      Detail::TBranchProxy fProxy;
      TDictionary*       fDict;
      TDictionary*       fContentDict; // type of content, if a collection
      ClassDef(TNamedBranchProxy, 0); // branch proxy with a name
   };

   // Used by TTreeReaderArray
   class TVirtualCollectionReader {
   public:
      TTreeReaderValueBase::EReadStatus fReadStatus;

      TVirtualCollectionReader() : fReadStatus(TTreeReaderValueBase::kReadNothingYet) {}

      virtual ~TVirtualCollectionReader();
      virtual size_t GetSize(Detail::TBranchProxy*) = 0;
      virtual void* At(Detail::TBranchProxy*, size_t /*idx*/) = 0;
   };

}
}

#endif // defined TTreeReaderUtils
