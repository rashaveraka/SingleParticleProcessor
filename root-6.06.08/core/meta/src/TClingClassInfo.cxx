// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TClingClassInfo

Emulation of the CINT ClassInfo class.

The CINT C++ interpreter provides an interface to metadata about
a class through the ClassInfo class.  This class provides the same
functionality, using an interface as close as possible to ClassInfo
but the class metadata comes from the Clang C++ compiler, not CINT.
*/

#include "TClingClassInfo.h"

#include "TClassEdit.h"
#include "TClingBaseClassInfo.h"
#include "TClingCallFunc.h"
#include "TClingMethodInfo.h"
#include "TDictionary.h"
#include "TClingTypeInfo.h"
#include "TError.h"
#include "TMetaUtils.h"
#include "ThreadLocalStorage.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/LookupHelper.h"
#include "cling/Utils/AST.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/DeclTemplate.h"
#include "clang/AST/GlobalDecl.h"
#include "clang/AST/PrettyPrinter.h"
#include "clang/AST/RecordLayout.h"
#include "clang/AST/Type.h"
#include "clang/Basic/Specifiers.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Sema/Sema.h"

#include "llvm/ExecutionEngine/GenericValue.h"
#include "llvm/Support/Casting.h"
#include "llvm/Support/raw_ostream.h"

#include <sstream>
#include <string>

using namespace clang;
using namespace ROOT;

static std::string FullyQualifiedName(const Decl *decl) {
   // Return the fully qualified name without worrying about normalizing it.
   std::string buf;
   if (const NamedDecl* ND = llvm::dyn_cast<NamedDecl>(decl)) {
      PrintingPolicy Policy(decl->getASTContext().getPrintingPolicy());
      llvm::raw_string_ostream stream(buf);
      ND->getNameForDiagnostic(stream, Policy, /*Qualified=*/true);
   }
   return buf;
}

TClingClassInfo::TClingClassInfo(cling::Interpreter *interp, Bool_t all)
   : fInterp(interp), fFirstTime(true), fDescend(false), fIterAll(all),
     fDecl(0), fType(0), fOffsetCache(0)
{
   TranslationUnitDecl *TU =
      interp->getCI()->getASTContext().getTranslationUnitDecl();
   // Could trigger deserialization of decls.
   cling::Interpreter::PushTransactionRAII RAII(interp);
   if (fIterAll)
      fIter = TU->decls_begin();
   else
      fIter = TU->noload_decls_begin();

   InternalNext();
   fFirstTime = true;
   // CINT had this odd behavior where a ClassInfo created without any
   // argument/input was set as an iterator that was ready to be iterated
   // on but was set an not IsValid *BUT* a few routine where using this
   // state as representing the global namespace (These routines include the
   // GetMethod routines and CallFunc::SetFunc, but do not include many others
   // (such as Property etc).  To be somewhat backward compatible, let's make
   // this state actually valid (i.e., representing both the ready-for-first-
   // iteration iterator *and* the global namespace) so that code that was
   // working with CINT (grabbing the default initialized ClassInfo
   // to look at the global namespace) is working again (and, yes, things that
   // used to not work like 'asking' the filename on this will go 'further'
   // but oh well).
   fDecl = TU;
   fType = 0;
}

TClingClassInfo::TClingClassInfo(cling::Interpreter *interp, const char *name)
   : fInterp(interp), fFirstTime(true), fDescend(false), fIterAll(kTRUE), fDecl(0),
     fType(0), fTitle(""), fOffsetCache(0)
{
   const cling::LookupHelper& lh = fInterp->getLookupHelper();
   const Type *type = 0;
   const Decl *decl = lh.findScope(name,
                                   gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                   : cling::LookupHelper::NoDiagnostics,
                                   &type, /* intantiateTemplate= */ true );
   if (!decl) {
      std::string buf = TClassEdit::InsertStd(name);
      if (buf != name) {
         decl = lh.findScope(buf,
                             gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                             : cling::LookupHelper::NoDiagnostics,
                             &type, /* intantiateTemplate= */ true );
      }
   }
   if (!decl && type) {
      const TagType *tagtype =type->getAs<TagType>();
      if (tagtype) {
         decl = tagtype->getDecl();
      }
   }
   fDecl = decl;
   fType = type;
   if (decl && decl->isInvalidDecl()) {
      Error("TClingClassInfo", "Found an invalid decl for %s.",name);
      fDecl = nullptr;
      fType = nullptr;
   }
}

TClingClassInfo::TClingClassInfo(cling::Interpreter *interp,
                                 const Type &tag)
   : fInterp(interp), fFirstTime(true), fDescend(false), fIterAll(kTRUE),
     fDecl(0), fType(0), fTitle(""), fOffsetCache(0)
{
   Init(tag);
}

void TClingClassInfo::AddBaseOffsetValue(const clang::Decl* decl, ptrdiff_t offset)
{
   // Add the offset value from this class to the non-virtual base class
   // determined by the parameter decl.

   OffsetPtrFunc_t executableFunc = 0;
   fOffsetCache[decl] = std::make_pair(offset, executableFunc);
}

long TClingClassInfo::ClassProperty() const
{
   if (!IsValid()) {
      return 0L;
   }
   long property = 0L;
   const RecordDecl *RD = llvm::dyn_cast<RecordDecl>(fDecl);
   if (!RD) {
      // We are an enum or namespace.
      // The cint interface always returns 0L for these guys.
      return property;
   }
   if (RD->isUnion()) {
      // The cint interface always returns 0L for these guys.
      return property;
   }
   // We now have a class or a struct.
   const CXXRecordDecl *CRD =
      llvm::dyn_cast<CXXRecordDecl>(fDecl);
   property |= kClassIsValid;
   if (CRD->isAbstract()) {
      property |= kClassIsAbstract;
   }
   if (CRD->hasUserDeclaredConstructor()) {
      property |= kClassHasExplicitCtor;
   }
   if (
      !CRD->hasUserDeclaredConstructor() &&
      !CRD->hasTrivialDefaultConstructor()
   ) {
      property |= kClassHasImplicitCtor;
   }
   if (
      CRD->hasUserProvidedDefaultConstructor() ||
      !CRD->hasTrivialDefaultConstructor()
   ) {
      property |= kClassHasDefaultCtor;
   }
   if (CRD->hasUserDeclaredDestructor()) {
      property |= kClassHasExplicitDtor;
   }
   else if (!CRD->hasTrivialDestructor()) {
      property |= kClassHasImplicitDtor;
   }
   if (CRD->hasUserDeclaredCopyAssignment()) {
      property |= kClassHasAssignOpr;
   }
   if (CRD->isPolymorphic()) {
      property |= kClassHasVirtual;
   }
   return property;
}

void TClingClassInfo::Delete(void *arena, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   // Invoke operator delete on a pointer to an object
   // of this class type.
   if (!IsValid()) {
      Error("TClingClassInfo::Delete()", "Called while invalid!");
      return;
   }
   if (!IsLoaded()) {
      Error("TClingClassInfo::Delete()", "Class is not loaded: %s",
            FullyQualifiedName(fDecl).c_str());
      return;
   }
   TClingCallFunc cf(fInterp,normCtxt);
   cf.ExecDestructor(this, arena, /*nary=*/0, /*withFree=*/true);
}

void TClingClassInfo::DeleteArray(void *arena, bool dtorOnly, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   // Invoke operator delete[] on a pointer to an array object
   // of this class type.
   if (!IsLoaded()) {
      return;
   }
   if (dtorOnly) {
      // There is no syntax in C++ for invoking the placement delete array
      // operator, so we have to placement destroy each element by hand.
      // Unfortunately we do not know how many elements to delete.
      //TClingCallFunc cf(fInterp);
      //cf.ExecDestructor(this, arena, nary, /*withFree=*/false);
      Error("DeleteArray", "Placement delete of an array is unsupported!\n");
      return;
   }
   TClingCallFunc cf(fInterp,normCtxt);
   cf.ExecDestructor(this, arena, /*nary=*/1, /*withFree=*/true);
}

void TClingClassInfo::Destruct(void *arena, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   // Invoke placement operator delete on a pointer to an array object
   // of this class type.
   if (!IsLoaded()) {
      return;
   }
   TClingCallFunc cf(fInterp,normCtxt);
   cf.ExecDestructor(this, arena, /*nary=*/0, /*withFree=*/false);
}

const FunctionTemplateDecl *TClingClassInfo::GetFunctionTemplate(const char *fname) const
{
   // Return any method or function in this scope with the name 'fname'.

   if (!IsLoaded()) {
      return 0;
   }

   if (fType) {
      const TypedefType *TT = llvm::dyn_cast<TypedefType>(fType);
      if (TT) {
         llvm::StringRef tname(TT->getDecl()->getName());
         if (tname.equals(fname)) {
            const NamedDecl *ndecl = llvm::dyn_cast<NamedDecl>(fDecl);
            if (ndecl && !ndecl->getName().equals(fname)) {
               // Constructor name matching the typedef type, use the decl name instead.
               return GetFunctionTemplate(ndecl->getName().str().c_str());
            }
         }
      }
   }
   const cling::LookupHelper &lh = fInterp->getLookupHelper();
   const FunctionTemplateDecl *fd
      = lh.findFunctionTemplate(fDecl, fname,
                                gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                : cling::LookupHelper::NoDiagnostics, false);
   if (fd) return fd->getCanonicalDecl();
   return 0;
}

const clang::ValueDecl *TClingClassInfo::GetDataMember(const char *name) const
{
   // Return the value decl (if any) corresponding to a data member which
   // the given name declared in this scope.

   const cling::LookupHelper &lh = fInterp->getLookupHelper();
   const ValueDecl *vd
      = lh.findDataMember(fDecl, name,
                          gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                          : cling::LookupHelper::NoDiagnostics);
   if (vd) return llvm::dyn_cast<ValueDecl>(vd->getCanonicalDecl());
   else return 0;
}

TClingMethodInfo TClingClassInfo::GetMethod(const char *fname) const
{
   // Return any method or function in this scope with the name 'fname'.

   if (!IsLoaded()) {
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }

   R__LOCKGUARD(gInterpreterMutex);

   if (fType) {
      const TypedefType *TT = llvm::dyn_cast<TypedefType>(fType);
      if (TT) {
         llvm::StringRef tname(TT->getDecl()->getName());
         if (tname.equals(fname)) {
            const NamedDecl *ndecl = llvm::dyn_cast<NamedDecl>(fDecl);
            if (ndecl && !ndecl->getName().equals(fname)) {
               // Constructor name matching the typedef type, use the decl name instead.
               return GetMethod(ndecl->getName().str().c_str());
            }
         }
      }
   }
   const cling::LookupHelper &lh = fInterp->getLookupHelper();
   const FunctionDecl *fd
      = lh.findAnyFunction(fDecl, fname,
                           gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                           : cling::LookupHelper::NoDiagnostics,
                           false);
   if (!fd) {
      // Function not found.
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }
   TClingMethodInfo tmi(fInterp);
   tmi.Init(fd);
   return tmi;
}

TClingMethodInfo TClingClassInfo::GetMethod(const char *fname,
      const char *proto, long *poffset, EFunctionMatchMode mode /*= kConversionMatch*/,
      EInheritanceMode imode /*= kWithInheritance*/) const
{
   return GetMethod(fname,proto,false,poffset,mode,imode);
}

TClingMethodInfo TClingClassInfo::GetMethod(const char *fname,
      const char *proto, bool objectIsConst,
      long *poffset, EFunctionMatchMode mode /*= kConversionMatch*/,
      EInheritanceMode imode /*= kWithInheritance*/) const
{
   if (poffset) {
      *poffset = 0L;
   }
   if (!IsLoaded()) {
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }

   R__LOCKGUARD(gInterpreterMutex);

   if (fType) {
      const TypedefType *TT = llvm::dyn_cast<TypedefType>(fType);
      if (TT) {
         llvm::StringRef tname(TT->getDecl()->getName());
         if (tname.equals(fname)) {
            const NamedDecl *ndecl = llvm::dyn_cast<NamedDecl>(fDecl);
            if (ndecl && !ndecl->getName().equals(fname)) {
               // Constructor name matching the typedef type, use the decl name instead.
               return GetMethod(ndecl->getName().str().c_str(),proto,
                                objectIsConst,poffset,
                                mode,imode);
            }
         }
      }

   }
   const cling::LookupHelper& lh = fInterp->getLookupHelper();
   const FunctionDecl *fd;
   if (mode == kConversionMatch) {
      fd = lh.findFunctionProto(fDecl, fname, proto,
                                gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                : cling::LookupHelper::NoDiagnostics,
                                objectIsConst);
   } else if (mode == kExactMatch) {
      fd = lh.matchFunctionProto(fDecl, fname, proto,
                                 gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                 : cling::LookupHelper::NoDiagnostics,
                                 objectIsConst);
   } else {
      Error("TClingClassInfo::GetMethod",
            "The MatchMode %d is not supported.", mode);
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }
   if (!fd) {
      // Function not found.
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }
   if (imode == TClingClassInfo::kInThisScope) {
      // If requested, check whether fd is a member function of this class.
      // Even though this seems to be the wrong order (we should not allow the
      // lookup to even collect candidates from the base) it does the right
      // thing: if any function overload exists in the derived class, all
      // (but explicitly used) will be hidden. Thus we will only find the
      // derived class's function overloads (or used, which is fine). Only
      // if there is none will we find those from the base, in which case
      // we will reject them here:
      const clang::DeclContext* ourDC = llvm::dyn_cast<clang::DeclContext>(fDecl);
      if (!fd->getDeclContext()->Equals(ourDC)
          && !(fd->getDeclContext()->isTransparentContext()
               && fd->getDeclContext()->getParent()->Equals(ourDC)))
         return TClingMethodInfo(fInterp);

      // The offset must be 0 - the function must be ours.
      if (poffset) *poffset = 0;
   } else {
      if (poffset) {
         // We have been asked to return a this pointer adjustment.
         if (const CXXMethodDecl *md =
             llvm::dyn_cast<CXXMethodDecl>(fd)) {
            // This is a class member function.
            *poffset = GetOffset(md);
         }
      }
   }
   TClingMethodInfo tmi(fInterp);
   tmi.Init(fd);
   return tmi;
}

TClingMethodInfo TClingClassInfo::GetMethod(const char *fname,
                                            const llvm::SmallVectorImpl<clang::QualType> &proto,
                                            long *poffset, EFunctionMatchMode mode /*= kConversionMatch*/,
                                            EInheritanceMode imode /*= kWithInheritance*/) const
{
   return GetMethod(fname,proto,false,poffset,mode,imode);
}

TClingMethodInfo TClingClassInfo::GetMethod(const char *fname,
                                            const llvm::SmallVectorImpl<clang::QualType> &proto, bool objectIsConst,
                                            long *poffset, EFunctionMatchMode mode /*= kConversionMatch*/,
                                            EInheritanceMode imode /*= kWithInheritance*/) const
{
   if (poffset) {
      *poffset = 0L;
   }
   if (!IsLoaded()) {
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }

   R__LOCKGUARD(gInterpreterMutex);

   if (fType) {
      const TypedefType *TT = llvm::dyn_cast<TypedefType>(fType);
      if (TT) {
         llvm::StringRef tname(TT->getDecl()->getName());
         if (tname.equals(fname)) {
            const NamedDecl *ndecl = llvm::dyn_cast<NamedDecl>(fDecl);
            if (ndecl && !ndecl->getName().equals(fname)) {
               // Constructor name matching the typedef type, use the decl name instead.
               return GetMethod(ndecl->getName().str().c_str(),proto,objectIsConst,poffset,
                                mode,imode);
            }
         }
      }

   }
   const cling::LookupHelper& lh = fInterp->getLookupHelper();
   const FunctionDecl *fd;
   if (mode == kConversionMatch) {
      fd = lh.findFunctionProto(fDecl, fname, proto,
                                gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                : cling::LookupHelper::NoDiagnostics,
                                objectIsConst);
   } else if (mode == kExactMatch) {
      fd = lh.matchFunctionProto(fDecl, fname, proto,
                                 gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                 : cling::LookupHelper::NoDiagnostics,
                                 objectIsConst);
   } else {
      Error("TClingClassInfo::GetMethod",
            "The MatchMode %d is not supported.", mode);
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }
   if (!fd) {
      // Function not found.
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }
   if (poffset) {
      // We have been asked to return a this pointer adjustment.
      if (const CXXMethodDecl *md =
          llvm::dyn_cast<CXXMethodDecl>(fd)) {
         // This is a class member function.
         *poffset = GetOffset(md);
      }
   }
   TClingMethodInfo tmi(fInterp);
   tmi.Init(fd);
   return tmi;
}

TClingMethodInfo TClingClassInfo::GetMethodWithArgs(const char *fname,
      const char *arglist, long *poffset, EFunctionMatchMode mode /* = kConversionMatch*/,
      EInheritanceMode imode /* = kWithInheritance*/) const
{
   return GetMethodWithArgs(fname,arglist,false,poffset,mode,imode);
}

TClingMethodInfo TClingClassInfo::GetMethodWithArgs(const char *fname,
      const char *arglist, bool objectIsConst,
      long *poffset, EFunctionMatchMode /*mode = kConversionMatch*/,
      EInheritanceMode /* imode = kWithInheritance*/) const
{

   R__LOCKGUARD(gInterpreterMutex);

   if (fType) {
      const TypedefType *TT = llvm::dyn_cast<TypedefType>(fType);
      if (TT) {
         llvm::StringRef tname(TT->getDecl()->getName());
         if (tname.equals(fname)) {
            const NamedDecl *ndecl = llvm::dyn_cast<NamedDecl>(fDecl);
            if (ndecl && !ndecl->getName().equals(fname)) {
               // Constructor name matching the typedef type, use the decl name instead.
               return GetMethod(ndecl->getName().str().c_str(),arglist,
                                objectIsConst,poffset
                                /* ,mode,imode */);
            }
         }
      }

   }
   if (poffset) {
      *poffset = 0L;
   }
   if (!IsLoaded()) {
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }
   if (!strcmp(arglist, ")")) {
      // CINT accepted a single right paren as meaning no arguments.
      arglist = "";
   }
   const cling::LookupHelper &lh = fInterp->getLookupHelper();
   const FunctionDecl *fd
      = lh.findFunctionArgs(fDecl, fname, arglist,
                            gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                            : cling::LookupHelper::NoDiagnostics,
                            objectIsConst);
   if (!fd) {
      // Function not found.
      TClingMethodInfo tmi(fInterp);
      return tmi;
   }
   if (poffset) {
     // We have been asked to return a this pointer adjustment.
     if (const CXXMethodDecl *md =
           llvm::dyn_cast<CXXMethodDecl>(fd)) {
        // This is a class member function.
        *poffset = GetOffset(md);
     }
   }
   TClingMethodInfo tmi(fInterp);
   tmi.Init(fd);
   return tmi;
}

int TClingClassInfo::GetMethodNArg(const char *method, const char *proto,
                                   Bool_t objectIsConst,
                                   EFunctionMatchMode mode /*= kConversionMatch*/) const
{
   // Note: Used only by TQObject.cxx:170 and only for interpreted classes.
   if (!IsLoaded()) {
      return -1;
   }

   R__LOCKGUARD(gInterpreterMutex);

   TClingMethodInfo mi = GetMethod(method, proto, objectIsConst, 0, mode);
   int clang_val = -1;
   if (mi.IsValid()) {
      unsigned num_params = mi.GetMethodDecl()->getNumParams();
      clang_val = static_cast<int>(num_params);
   }
   return clang_val;
}

long TClingClassInfo::GetOffset(const CXXMethodDecl* md) const
{

   R__LOCKGUARD(gInterpreterMutex);

   long offset = 0L;
   const CXXRecordDecl* definer = md->getParent();
   const CXXRecordDecl* accessor =
      llvm::cast<CXXRecordDecl>(fDecl);
   if (definer != accessor) {
      // This function may not be accessible using a pointer
      // to the declaring class, get the adjustment necessary
      // to convert that to a pointer to the defining class.
      TClingBaseClassInfo bi(fInterp, const_cast<TClingClassInfo*>(this));
      while (bi.Next(0)) {
         TClingClassInfo* bci = bi.GetBase();
         if (bci->GetDecl() == definer) {
            // We have found the right base class, now get the
            // necessary adjustment.
            offset = bi.Offset();
            break;
         }
      }
   }
   return offset;
}

ptrdiff_t TClingClassInfo::GetBaseOffset(TClingClassInfo* base, void* address, bool isDerivedObject)
{

   R__LOCKGUARD(gInterpreterMutex);

   // Check for the offset in the cache.
   auto iter = fOffsetCache.find(base->GetDecl());
   if (iter != fOffsetCache.end()) {
      std::pair<ptrdiff_t, OffsetPtrFunc_t> offsetCache = (*iter).second;
      if (OffsetPtrFunc_t executableFunc = offsetCache.second) {
         if (address) {
            return (*executableFunc)(address, isDerivedObject);
         }
         else {
            Error("TClingBaseClassInfo::Offset", "The address of the object for virtual base offset calculation is not valid.");
            return -1;
         }
      }
      else {
         return offsetCache.first;
      }
   }

   // Compute the offset.
   TClingBaseClassInfo binfo(fInterp, this, base);
   return binfo.Offset(address, isDerivedObject);
}

static bool HasBody(const clang::FunctionDecl &decl, const cling::Interpreter &interp)
{
   if (decl.hasBody()) return true;

   GlobalDecl GD;
   if (const CXXConstructorDecl* Ctor = dyn_cast<CXXConstructorDecl>(&decl))
     GD = GlobalDecl(Ctor, Ctor_Complete);
   else if (const CXXDestructorDecl* Dtor = dyn_cast<CXXDestructorDecl>(&decl))
     GD = GlobalDecl(Dtor, Dtor_Deleting);
   else
     GD = GlobalDecl(&decl);
   std::string mangledName;
   cling::utils::Analyze::maybeMangleDeclName(GD, mangledName);

   void *GV = interp.getAddressOfGlobal(mangledName.c_str());
   if (GV) return true;

   return false;
}

bool TClingClassInfo::HasDefaultConstructor() const
{
   // Return true if there a constructor taking no arguments (including
   // a constructor that has defaults for all of its arguments) which
   // is callable.  Either it has a body, or it is trivial and the
   // compiler elides it.
   //
   // Note: This is could enhanced to also know about the ROOT ioctor
   // but this was not the case in CINT.
   //
   if (!IsLoaded()) {
      return false;
   }
   const CXXRecordDecl* CRD = llvm::dyn_cast<CXXRecordDecl>(fDecl);
   if (!CRD) {
      // Namespaces do not have constructors.
      return false;
   }
   if (CRD->hasTrivialDefaultConstructor()) {
      // This class has a default constructor that can be called,
      // but has no body.
      return true;
   }
   // Note: This iteration may force template instantiations!
   cling::Interpreter::PushTransactionRAII pushedT(fInterp);
   for (CXXRecordDecl::ctor_iterator I = CRD->ctor_begin(),
         E = CRD->ctor_end(); I != E; ++I) {
      if (I->getMinRequiredArguments() == 0) {
         if ((I->getAccess() == AS_public) && HasBody(**I,*fInterp)) {
            return true;
         }
         if (I->isTemplateInstantiation()) {
            const clang::FunctionDecl* FD =
               I->getInstantiatedFromMemberFunction();
            if ((FD->getAccess() == AS_public) && HasBody(*FD,*fInterp)) {
               return true;
            }
         }
      }
   }
   return false;
}

bool TClingClassInfo::HasMethod(const char *name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   if (IsLoaded() && !llvm::isa<EnumDecl>(fDecl)) {
      return fInterp->getLookupHelper()
         .hasFunction(fDecl, name,
                      gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                      : cling::LookupHelper::NoDiagnostics);
   }
   return false;
}

void TClingClassInfo::Init(const char *name)
{
   fFirstTime = true;
   fDescend = false;
   fIter = DeclContext::decl_iterator();
   fDecl = 0;
   fType = 0;
   fIterStack.clear();
   const cling::LookupHelper& lh = fInterp->getLookupHelper();
   fDecl = lh.findScope(name, gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                          : cling::LookupHelper::NoDiagnostics,
                        &fType, /* intantiateTemplate= */ true );
   if (!fDecl) {
      std::string buf = TClassEdit::InsertStd(name);
      if (buf != name) {
         fDecl = lh.findScope(buf, gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                              : cling::LookupHelper::NoDiagnostics,
                              &fType, /* intantiateTemplate= */ true );
      }
   }
   if (!fDecl && fType) {
      const TagType *tagtype =fType->getAs<TagType>();
      if (tagtype) {
         fDecl = tagtype->getDecl();
      }
   }
}

void TClingClassInfo::Init(const Decl* decl)
{
   fFirstTime = true;
   fDescend = false;
   fIter = DeclContext::decl_iterator();
   fDecl = decl;
   fType = 0;
   fIterStack.clear();
}

void TClingClassInfo::Init(int tagnum)
{
   Fatal("TClingClassInfo::Init(tagnum)", "Should no longer be called");
   return;
}

void TClingClassInfo::Init(const Type &tag)
{
   fType = &tag;

   R__LOCKGUARD(gInterpreterMutex);

   const TagType *tagtype = fType->getAs<TagType>();
   if (tagtype) {
      fDecl = tagtype->getDecl();
   }
   else {
      fDecl = 0;
   }
   if (!fDecl) {
      QualType qType(fType,0);
      static PrintingPolicy printPol(fInterp->getCI()->getLangOpts());
      printPol.SuppressScope = false;
      Error("TClingClassInfo::Init(const Type&)",
            "The given type %s does not point to a Decl",
            qType.getAsString(printPol).c_str());
   }
}

bool TClingClassInfo::IsBase(const char *name) const
{
   if (!IsLoaded()) {
      return false;
   }
   TClingClassInfo base(fInterp, name);
   if (!base.IsValid()) {
      return false;
   }

   R__LOCKGUARD(gInterpreterMutex);

   const CXXRecordDecl *CRD =
      llvm::dyn_cast<CXXRecordDecl>(fDecl);
   if (!CRD) {
      // We are an enum, namespace, or translation unit,
      // we cannot be the base of anything.
      return false;
   }
   const CXXRecordDecl *baseCRD =
      llvm::dyn_cast<CXXRecordDecl>(base.GetDecl());
   return CRD->isDerivedFrom(baseCRD);
}

bool TClingClassInfo::IsEnum(cling::Interpreter *interp, const char *name)
{
   // Note: This is a static member function.
   TClingClassInfo info(interp, name);
   if (info.IsValid() && (info.Property() & kIsEnum)) {
      return true;
   }
   return false;
}

bool TClingClassInfo::IsLoaded() const
{
   // IsLoaded in CINT was meaning is known to the interpreter
   // and has a complete definition.
   // IsValid in Cling (as in CING) means 'just' is known to the
   // interpreter.
   if (!IsValid()) {
      return false;
   }
   if (fDecl == 0) {
      return false;
   }

   R__LOCKGUARD(gInterpreterMutex);

   const CXXRecordDecl *CRD = llvm::dyn_cast<CXXRecordDecl>(fDecl);
   if ( CRD ) {
      if (!CRD->hasDefinition()) {
         return false;
      }
   } else {
      const TagDecl *TD = llvm::dyn_cast<TagDecl>(fDecl);
      if (TD && TD->getDefinition() == 0) {
         return false;
      }
   }
   // All clang classes are considered loaded.
   return true;
}

bool TClingClassInfo::IsValid() const
{
   return fDecl;
}

bool TClingClassInfo::IsValidMethod(const char *method, const char *proto,
                                    Bool_t objectIsConst,
                                    long *offset,
                                    EFunctionMatchMode mode /*= kConversionMatch*/) const
{
   // Check if the method with the given prototype exist.
   if (!IsLoaded()) {
      return false;
   }
   if (offset) {
      *offset = 0L;
   }
   TClingMethodInfo mi = GetMethod(method, proto, offset, mode);
   return mi.IsValid();
}

int TClingClassInfo::InternalNext()
{

   R__LOCKGUARD(gInterpreterMutex);

   if (!*fIter) {
      // Iterator is already invalid.
      if (fFirstTime && fDecl) {
         std::string buf;
         if (const NamedDecl* ND =
               llvm::dyn_cast<NamedDecl>(fDecl)) {
            PrintingPolicy Policy(fDecl->getASTContext().
               getPrintingPolicy());
            llvm::raw_string_ostream stream(buf);
            ND->getNameForDiagnostic(stream, Policy, /*Qualified=*/false);
         }
         Error("TClingClassInfo::InternalNext",
            "Next called but iteration not prepared for %s!", buf.c_str());
      }
      return 0;
   }
   cling::Interpreter::PushTransactionRAII pushedT(fInterp);
   while (true) {
      // Advance to next usable decl, or return if there is no next usable decl.
      if (fFirstTime) {
         // The cint semantics are strange.
         fFirstTime = false;
      }
      else {
         // Advance the iterator one decl, descending into the current decl
         // context if necessary.
         if (!fDescend) {
            // Do not need to scan the decl context of the current decl,
            // move on to the next decl.
            ++fIter;
         }
         else {
            // Descend into the decl context of the current decl.
            fDescend = false;
            //fprintf(stderr,
            //   "TClingClassInfo::InternalNext:  "
            //   "pushing ...\n");
            fIterStack.push_back(fIter);
            DeclContext *DC = llvm::cast<DeclContext>(*fIter);
            if (fIterAll)
               fIter = DC->decls_begin();
            else
               fIter = DC->noload_decls_begin();
         }
         // Fix it if we went past the end.
         while (!*fIter && fIterStack.size()) {
            //fprintf(stderr,
            //   "TClingClassInfo::InternalNext:  "
            //   "popping ...\n");
            fIter = fIterStack.back();
            fIterStack.pop_back();
            ++fIter;
         }
         // Check for final termination.
         if (!*fIter) {
            // We have reached the end of the translation unit, all done.
            fDecl = 0;
            fType = 0;
            return 0;
         }
      }
      // Return if this decl is a class, struct, union, enum, or namespace.
      Decl::Kind DK = fIter->getKind();
      if ((DK == Decl::Namespace) || (DK == Decl::Enum) ||
            (DK == Decl::CXXRecord) ||
            (DK == Decl::ClassTemplateSpecialization)) {
         const TagDecl *TD = llvm::dyn_cast<TagDecl>(*fIter);
         if (TD && !TD->isCompleteDefinition()) {
            // For classes and enums, stop only on definitions.
            continue;
         }
         if (DK == Decl::Namespace) {
            // For namespaces, stop only on the first definition.
            if (!fIter->isCanonicalDecl()) {
               // Not the first definition.
               fDescend = true;
               continue;
            }
         }
         if (DK != Decl::Enum) {
            // We do not descend into enums.
            DeclContext *DC = llvm::cast<DeclContext>(*fIter);
            if ((fIterAll && *DC->decls_begin())
                || (!fIterAll && *DC->noload_decls_begin())) {
               // Next iteration will begin scanning the decl context
               // contained by this decl.
               fDescend = true;
            }
         }
         // Iterator is now valid.
         fDecl = *fIter;
         fType = 0;
         if (fDecl) {
            if (fDecl->isInvalidDecl()) {
               Warning("TClingClassInfo::Next()","Reached an invalid decl.");
            }
            if (const RecordDecl *RD =
                  llvm::dyn_cast<RecordDecl>(fDecl)) {
               fType = RD->getASTContext().getRecordType(RD).getTypePtr();
            }
         }
         return 1;
      }
   }
}

int TClingClassInfo::Next()
{
   return InternalNext();
}

void *TClingClassInfo::New(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   // Invoke a new expression to use the class constructor
   // that takes no arguments to create an object of this class type.
   if (!IsValid()) {
      Error("TClingClassInfo::New()", "Called while invalid!");
      return 0;
   }
   if (!IsLoaded()) {
      Error("TClingClassInfo::New()", "Class is not loaded: %s",
            FullyQualifiedName(fDecl).c_str());
      return 0;
   }
   {
      R__LOCKGUARD(gInterpreterMutex);
      const CXXRecordDecl* RD = dyn_cast<CXXRecordDecl>(fDecl);
      if (!RD) {
         Error("TClingClassInfo::New()", "This is a namespace!: %s",
               FullyQualifiedName(fDecl).c_str());
         return 0;
      }
      if (!HasDefaultConstructor()) {
         // FIXME: We fail roottest root/io/newdelete if we issue this message!
         //Error("TClingClassInfo::New()", "Class has no default constructor: %s",
         //      FullyQualifiedName(fDecl).c_str());
         return 0;
      }
   } // End of Lock section.
   void* obj = 0;
   TClingCallFunc cf(fInterp,normCtxt);
   obj = cf.ExecDefaultConstructor(this, /*address=*/0, /*nary=*/0);
   if (!obj) {
      Error("TClingClassInfo::New()", "Call of default constructor "
            "failed to return an object for class: %s",
            FullyQualifiedName(fDecl).c_str());
      return 0;
   }
   return obj;
}

void *TClingClassInfo::New(int n, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   // Invoke a new expression to use the class constructor
   // that takes no arguments to create an array object
   // of this class type.
   if (!IsValid()) {
      Error("TClingClassInfo::New(n)", "Called while invalid!");
      return 0;
   }
   if (!IsLoaded()) {
      Error("TClingClassInfo::New(n)", "Class is not loaded: %s",
            FullyQualifiedName(fDecl).c_str());
      return 0;
   }

   {
      R__LOCKGUARD(gInterpreterMutex);

      const CXXRecordDecl* RD = dyn_cast<CXXRecordDecl>(fDecl);
      if (!RD) {
         Error("TClingClassInfo::New(n)", "This is a namespace!: %s",
               FullyQualifiedName(fDecl).c_str());
         return 0;
      }
      if (!HasDefaultConstructor()) {
         // FIXME: We fail roottest root/io/newdelete if we issue this message!
         //Error("TClingClassInfo::New(n)",
         //      "Class has no default constructor: %s",
         //      FullyQualifiedName(fDecl).c_str());
         return 0;
      }
   } // End of Lock section.
   void* obj = 0;
   TClingCallFunc cf(fInterp,normCtxt);
   obj = cf.ExecDefaultConstructor(this, /*address=*/0,
                                   /*nary=*/(unsigned long)n);
   if (!obj) {
      Error("TClingClassInfo::New(n)", "Call of default constructor "
            "failed to return an array of class: %s",
            FullyQualifiedName(fDecl).c_str());
      return 0;
   }
   return obj;
}

void *TClingClassInfo::New(int n, void *arena, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   // Invoke a placement new expression to use the class
   // constructor that takes no arguments to create an
   // array of objects of this class type in the given
   // memory arena.
   if (!IsValid()) {
      Error("TClingClassInfo::New(n, arena)", "Called while invalid!");
      return 0;
   }
   if (!IsLoaded()) {
      Error("TClingClassInfo::New(n, arena)", "Class is not loaded: %s",
            FullyQualifiedName(fDecl).c_str());
      return 0;
   }
   {
      R__LOCKGUARD(gInterpreterMutex);

      const CXXRecordDecl* RD = dyn_cast<CXXRecordDecl>(fDecl);
      if (!RD) {
         Error("TClingClassInfo::New(n, arena)", "This is a namespace!: %s",
               FullyQualifiedName(fDecl).c_str());
         return 0;
      }
      if (!HasDefaultConstructor()) {
         // FIXME: We fail roottest root/io/newdelete if we issue this message!
         //Error("TClingClassInfo::New(n, arena)",
         //      "Class has no default constructor: %s",
         //      FullyQualifiedName(fDecl).c_str());
         return 0;
      }
   } // End of Lock section
   void* obj = 0;
   TClingCallFunc cf(fInterp,normCtxt);
   // Note: This will always return arena.
   obj = cf.ExecDefaultConstructor(this, /*address=*/arena,
                                   /*nary=*/(unsigned long)n);
   return obj;
}

void *TClingClassInfo::New(void *arena, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   // Invoke a placement new expression to use the class
   // constructor that takes no arguments to create an
   // object of this class type in the given memory arena.
   if (!IsValid()) {
      Error("TClingClassInfo::New(arena)", "Called while invalid!");
      return 0;
   }
   if (!IsLoaded()) {
      Error("TClingClassInfo::New(arena)", "Class is not loaded: %s",
            FullyQualifiedName(fDecl).c_str());
      return 0;
   }
   {
      R__LOCKGUARD(gInterpreterMutex);

      const CXXRecordDecl* RD = dyn_cast<CXXRecordDecl>(fDecl);
      if (!RD) {
         Error("TClingClassInfo::New(arena)", "This is a namespace!: %s",
               FullyQualifiedName(fDecl).c_str());
         return 0;
      }
      if (!HasDefaultConstructor()) {
         // FIXME: We fail roottest root/io/newdelete if we issue this message!
         //Error("TClingClassInfo::New(arena)",
         //      "Class has no default constructor: %s",
         //      FullyQualifiedName(fDecl).c_str());
         return 0;
      }
   } // End of Locked section.
   void* obj = 0;
   TClingCallFunc cf(fInterp,normCtxt);
   // Note: This will always return arena.
   obj = cf.ExecDefaultConstructor(this, /*address=*/arena, /*nary=*/0);
   return obj;
}

long TClingClassInfo::Property() const
{
   if (!IsValid()) {
      return 0L;
   }

   R__LOCKGUARD(gInterpreterMutex);

   long property = 0L;
   property |= kIsCPPCompiled;
   const clang::DeclContext *ctxt = fDecl->getDeclContext();
   clang::NamespaceDecl *std_ns =fInterp->getSema().getStdNamespace();
   while (! ctxt->isTranslationUnit())  {
      if (ctxt->Equals(std_ns)) {
         property |= kIsDefinedInStd;
         break;
      }
      ctxt = ctxt->getParent();
   }
   Decl::Kind DK = fDecl->getKind();
   if ((DK == Decl::Namespace) || (DK == Decl::TranslationUnit)) {
      property |= kIsNamespace;
      return property;
   }
   // Note: Now we have class, enum, struct, union only.
   const TagDecl *TD = llvm::dyn_cast<TagDecl>(fDecl);
   if (!TD) {
      return 0L;
   }
   if (TD->isEnum()) {
      property |= kIsEnum;
      return property;
   }
   // Note: Now we have class, struct, union only.
   const CXXRecordDecl *CRD =
      llvm::dyn_cast<CXXRecordDecl>(fDecl);
   if (CRD->isClass()) {
      property |= kIsClass;
   }
   else if (CRD->isStruct()) {
      property |= kIsStruct;
   }
   else if (CRD->isUnion()) {
      property |= kIsUnion;
   }
   if (CRD->hasDefinition() && CRD->isAbstract()) {
      property |= kIsAbstract;
   }
   return property;
}

int TClingClassInfo::RootFlag() const
{
   if (!IsValid()) {
      return 0;
   }
   // FIXME: Implement this when rootcling provides the value.
   return 0;
}

int TClingClassInfo::Size() const
{
   if (!IsValid()) {
      return -1;
   }
   if (!fDecl) {
      // A forward declared class.
      return 0;
   }

   R__LOCKGUARD(gInterpreterMutex);

   Decl::Kind DK = fDecl->getKind();
   if (DK == Decl::Namespace) {
      // Namespaces are special for cint.
      return 1;
   }
   else if (DK == Decl::Enum) {
      // Enums are special for cint.
      return 0;
   }
   const RecordDecl *RD = llvm::dyn_cast<RecordDecl>(fDecl);
   if (!RD) {
      // Should not happen.
      return -1;
   }
   if (!RD->getDefinition()) {
      // Forward-declared class.
      return 0;
   }
   ASTContext &Context = fDecl->getASTContext();
   cling::Interpreter::PushTransactionRAII RAII(fInterp);
   const ASTRecordLayout &Layout = Context.getASTRecordLayout(RD);
   int64_t size = Layout.getSize().getQuantity();
   int clang_size = static_cast<int>(size);
   return clang_size;
}

long TClingClassInfo::Tagnum() const
{
   if (!IsValid()) {
      return -1L;
   }
   return reinterpret_cast<long>(fDecl);
}

const char *TClingClassInfo::FileName()
{
   if (!IsValid()) {
      return 0;
   }
   fDeclFileName = ROOT::TMetaUtils::GetFileName(*GetDecl(), *fInterp);
   return fDeclFileName.c_str();
}

void TClingClassInfo::FullName(std::string &output, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   // Return QualifiedName.
   output.clear();
   if (!IsValid()) {
      return;
   }
   if (fType) {
      QualType type(fType, 0);
      ROOT::TMetaUtils::GetNormalizedName(output, type, *fInterp, normCtxt);
   }
   else {
      if (const NamedDecl* ND =
            llvm::dyn_cast<NamedDecl>(fDecl)) {
         PrintingPolicy Policy(fDecl->getASTContext().
            getPrintingPolicy());
         llvm::raw_string_ostream stream(output);
         ND->getNameForDiagnostic(stream, Policy, /*Qualified=*/true);
      }
   }
}

const char *TClingClassInfo::Name() const
{
   // Return unqualified name.
   if (!IsValid()) {
      return 0;
   }
   // Note: This *must* be static/thread_local because we are returning a pointer inside it!
   TTHREAD_TLS_DECL( std::string, buf);

   buf.clear();
   if (const NamedDecl* ND = llvm::dyn_cast<NamedDecl>(fDecl)) {
      PrintingPolicy Policy(fDecl->getASTContext().getPrintingPolicy());
      llvm::raw_string_ostream stream(buf);
      ND->getNameForDiagnostic(stream, Policy, /*Qualified=*/false);
   }
   return buf.c_str();
}

const char *TClingClassInfo::Title()
{
   if (!IsValid()) {
      return 0;
   }
   // NOTE: We cannot cache the result, since we are really an iterator.
   // Try to get the comment either from the annotation or the header
   // file, if present.
   // Iterate over the redeclarations, we can have multiple definitions in the
   // redecl chain (came from merging of pcms).

   R__LOCKGUARD(gInterpreterMutex);

   if (const TagDecl *TD = llvm::dyn_cast<TagDecl>(GetDecl())) {
      if ( (TD = ROOT::TMetaUtils::GetAnnotatedRedeclarable(TD)) ) {
         if (AnnotateAttr *A = TD->getAttr<AnnotateAttr>()) {
            std::string attr = A->getAnnotation().str();
            if (attr.find(TMetaUtils::propNames::separator) != std::string::npos) {
               if (TMetaUtils::ExtractAttrPropertyFromName(*TD,TMetaUtils::propNames::comment,attr)) {
                  fTitle = attr;
                  return fTitle.c_str();
               }
            } else {
               fTitle = attr;
               return fTitle.c_str();
            }
         }
      }
   }
   // Try to get the comment from the header file, if present.
   // but not for decls from AST file, where rootcling would have
   // created an annotation
   const CXXRecordDecl *CRD =
      llvm::dyn_cast<CXXRecordDecl>(GetDecl());
   if (CRD && !CRD->isFromASTFile()) {
      fTitle = ROOT::TMetaUtils::GetClassComment(*CRD,0,*fInterp).str();
   }
   return fTitle.c_str();
}

const char *TClingClassInfo::TmpltName() const
{
   if (!IsValid()) {
      return 0;
   }

   R__LOCKGUARD(gInterpreterMutex);

   // Note: This *must* be static/thread_local because we are returning a pointer inside it!
   TTHREAD_TLS_DECL( std::string, buf);
   buf.clear();
   if (const NamedDecl* ND = llvm::dyn_cast<NamedDecl>(fDecl)) {
      // Note: This does *not* include the template arguments!
      buf = ND->getNameAsString();
   }
   return buf.c_str();
}

