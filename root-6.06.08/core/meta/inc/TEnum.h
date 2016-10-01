// @(#)root/meta:$Id$
// Author: Bianca-Cristina Cristescu   09/07/13

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEnum
#define ROOT_TEnum

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEnum                                                                //
//                                                                      //
// TEnum class defines enum type.                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif
#ifndef ROOT_THashList
#include "THashList.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif
#ifndef ROOT_TDictionary
#include "TDictionary.h"
#endif
#ifndef ROOT_Rtypeinfo
#include "Rtypeinfo.h"
#endif

class TClass;
class TEnumConstant;

class TEnum : public TDictionary {

private:
   THashList fConstantList;     //list of constants the enum type
   void     *fInfo;             //!interpreter implementation provided declaration
   TClass   *fClass;            //!owning class
   std::string fQualName;       // fully qualified type name

public:

   enum ESearchAction {kNone                 = 0,
                       kAutoload             = 1,
                       kInterpLookup         = 2,
                       kALoadAndInterpLookup = 3
                      };

   TEnum(): fInfo(0), fClass(0) {}
   TEnum(const char *name, void *info, TClass *cls);
   virtual ~TEnum();

   void                  AddConstant(TEnumConstant *constant);
   TClass               *GetClass() const {
      return fClass;
   }
   const TSeqCollection *GetConstants() const {
      return &fConstantList;
   }
   const TEnumConstant  *GetConstant(const char *name) const {
      return (TEnumConstant *) fConstantList.FindObject(name);
   }
   DeclId_t              GetDeclId() const {
      return (DeclId_t)fInfo;
   }
   Bool_t                IsValid();
   Long_t                Property() const;
   void                  SetClass(TClass *cl) {
      fClass = cl;
   }
   void                  Update(DeclId_t id);
   const char*           GetQualifiedName() const { return fQualName.c_str(); }
   static TEnum         *GetEnum(const std::type_info &ti, ESearchAction sa = kALoadAndInterpLookup);
   static TEnum         *GetEnum(const char *enumName, ESearchAction sa = kALoadAndInterpLookup);

   ClassDef(TEnum, 2) //Enum type class
};

#endif
