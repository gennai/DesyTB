// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME OnlineMonDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardHit.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardCluster.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardPlane.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonConfiguration.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/CheckEOF.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *SimpleStandardHit_Dictionary();
   static void SimpleStandardHit_TClassManip(TClass*);
   static void delete_SimpleStandardHit(void *p);
   static void deleteArray_SimpleStandardHit(void *p);
   static void destruct_SimpleStandardHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimpleStandardHit*)
   {
      ::SimpleStandardHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SimpleStandardHit));
      static ::ROOT::TGenericClassInfo 
         instance("SimpleStandardHit", "SimpleStandardHit.hh", 11,
                  typeid(::SimpleStandardHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SimpleStandardHit_Dictionary, isa_proxy, 0,
                  sizeof(::SimpleStandardHit) );
      instance.SetDelete(&delete_SimpleStandardHit);
      instance.SetDeleteArray(&deleteArray_SimpleStandardHit);
      instance.SetDestructor(&destruct_SimpleStandardHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimpleStandardHit*)
   {
      return GenerateInitInstanceLocal((::SimpleStandardHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SimpleStandardHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SimpleStandardHit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SimpleStandardHit*)0x0)->GetClass();
      SimpleStandardHit_TClassManip(theClass);
   return theClass;
   }

   static void SimpleStandardHit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SimpleStandardCluster_Dictionary();
   static void SimpleStandardCluster_TClassManip(TClass*);
   static void *new_SimpleStandardCluster(void *p = 0);
   static void *newArray_SimpleStandardCluster(Long_t size, void *p);
   static void delete_SimpleStandardCluster(void *p);
   static void deleteArray_SimpleStandardCluster(void *p);
   static void destruct_SimpleStandardCluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimpleStandardCluster*)
   {
      ::SimpleStandardCluster *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SimpleStandardCluster));
      static ::ROOT::TGenericClassInfo 
         instance("SimpleStandardCluster", "SimpleStandardCluster.hh", 13,
                  typeid(::SimpleStandardCluster), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SimpleStandardCluster_Dictionary, isa_proxy, 0,
                  sizeof(::SimpleStandardCluster) );
      instance.SetNew(&new_SimpleStandardCluster);
      instance.SetNewArray(&newArray_SimpleStandardCluster);
      instance.SetDelete(&delete_SimpleStandardCluster);
      instance.SetDeleteArray(&deleteArray_SimpleStandardCluster);
      instance.SetDestructor(&destruct_SimpleStandardCluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimpleStandardCluster*)
   {
      return GenerateInitInstanceLocal((::SimpleStandardCluster*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SimpleStandardCluster*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SimpleStandardCluster_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SimpleStandardCluster*)0x0)->GetClass();
      SimpleStandardCluster_TClassManip(theClass);
   return theClass;
   }

   static void SimpleStandardCluster_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *OnlineMonConfiguration_Dictionary();
   static void OnlineMonConfiguration_TClassManip(TClass*);
   static void *new_OnlineMonConfiguration(void *p = 0);
   static void *newArray_OnlineMonConfiguration(Long_t size, void *p);
   static void delete_OnlineMonConfiguration(void *p);
   static void deleteArray_OnlineMonConfiguration(void *p);
   static void destruct_OnlineMonConfiguration(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::OnlineMonConfiguration*)
   {
      ::OnlineMonConfiguration *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::OnlineMonConfiguration));
      static ::ROOT::TGenericClassInfo 
         instance("OnlineMonConfiguration", "OnlineMonConfiguration.hh", 19,
                  typeid(::OnlineMonConfiguration), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &OnlineMonConfiguration_Dictionary, isa_proxy, 0,
                  sizeof(::OnlineMonConfiguration) );
      instance.SetNew(&new_OnlineMonConfiguration);
      instance.SetNewArray(&newArray_OnlineMonConfiguration);
      instance.SetDelete(&delete_OnlineMonConfiguration);
      instance.SetDeleteArray(&deleteArray_OnlineMonConfiguration);
      instance.SetDestructor(&destruct_OnlineMonConfiguration);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::OnlineMonConfiguration*)
   {
      return GenerateInitInstanceLocal((::OnlineMonConfiguration*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::OnlineMonConfiguration*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *OnlineMonConfiguration_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::OnlineMonConfiguration*)0x0)->GetClass();
      OnlineMonConfiguration_TClassManip(theClass);
   return theClass;
   }

   static void OnlineMonConfiguration_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SimpleStandardPlane_Dictionary();
   static void SimpleStandardPlane_TClassManip(TClass*);
   static void delete_SimpleStandardPlane(void *p);
   static void deleteArray_SimpleStandardPlane(void *p);
   static void destruct_SimpleStandardPlane(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimpleStandardPlane*)
   {
      ::SimpleStandardPlane *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SimpleStandardPlane));
      static ::ROOT::TGenericClassInfo 
         instance("SimpleStandardPlane", "SimpleStandardPlane.hh", 28,
                  typeid(::SimpleStandardPlane), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SimpleStandardPlane_Dictionary, isa_proxy, 0,
                  sizeof(::SimpleStandardPlane) );
      instance.SetDelete(&delete_SimpleStandardPlane);
      instance.SetDeleteArray(&deleteArray_SimpleStandardPlane);
      instance.SetDestructor(&destruct_SimpleStandardPlane);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimpleStandardPlane*)
   {
      return GenerateInitInstanceLocal((::SimpleStandardPlane*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SimpleStandardPlane*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SimpleStandardPlane_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SimpleStandardPlane*)0x0)->GetClass();
      SimpleStandardPlane_TClassManip(theClass);
   return theClass;
   }

   static void SimpleStandardPlane_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *BaseCollection_Dictionary();
   static void BaseCollection_TClassManip(TClass*);
   static void delete_BaseCollection(void *p);
   static void deleteArray_BaseCollection(void *p);
   static void destruct_BaseCollection(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BaseCollection*)
   {
      ::BaseCollection *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BaseCollection));
      static ::ROOT::TGenericClassInfo 
         instance("BaseCollection", "BaseCollection.hh", 44,
                  typeid(::BaseCollection), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &BaseCollection_Dictionary, isa_proxy, 1,
                  sizeof(::BaseCollection) );
      instance.SetDelete(&delete_BaseCollection);
      instance.SetDeleteArray(&deleteArray_BaseCollection);
      instance.SetDestructor(&destruct_BaseCollection);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BaseCollection*)
   {
      return GenerateInitInstanceLocal((::BaseCollection*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::BaseCollection*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *BaseCollection_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::BaseCollection*)0x0)->GetClass();
      BaseCollection_TClassManip(theClass);
   return theClass;
   }

   static void BaseCollection_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CheckEOF_Dictionary();
   static void CheckEOF_TClassManip(TClass*);
   static void *new_CheckEOF(void *p = 0);
   static void *newArray_CheckEOF(Long_t size, void *p);
   static void delete_CheckEOF(void *p);
   static void deleteArray_CheckEOF(void *p);
   static void destruct_CheckEOF(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CheckEOF*)
   {
      ::CheckEOF *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CheckEOF));
      static ::ROOT::TGenericClassInfo 
         instance("CheckEOF", "CheckEOF.hh", 14,
                  typeid(::CheckEOF), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CheckEOF_Dictionary, isa_proxy, 1,
                  sizeof(::CheckEOF) );
      instance.SetNew(&new_CheckEOF);
      instance.SetNewArray(&newArray_CheckEOF);
      instance.SetDelete(&delete_CheckEOF);
      instance.SetDeleteArray(&deleteArray_CheckEOF);
      instance.SetDestructor(&destruct_CheckEOF);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CheckEOF*)
   {
      return GenerateInitInstanceLocal((::CheckEOF*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::CheckEOF*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CheckEOF_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CheckEOF*)0x0)->GetClass();
      CheckEOF_TClassManip(theClass);
   return theClass;
   }

   static void CheckEOF_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void delete_OnlineMonWindow(void *p);
   static void deleteArray_OnlineMonWindow(void *p);
   static void destruct_OnlineMonWindow(void *p);
   static void streamer_OnlineMonWindow(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::OnlineMonWindow*)
   {
      ::OnlineMonWindow *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::OnlineMonWindow >(0);
      static ::ROOT::TGenericClassInfo 
         instance("OnlineMonWindow", ::OnlineMonWindow::Class_Version(), "OnlineMonWindow.hh", 46,
                  typeid(::OnlineMonWindow), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::OnlineMonWindow::Dictionary, isa_proxy, 16,
                  sizeof(::OnlineMonWindow) );
      instance.SetDelete(&delete_OnlineMonWindow);
      instance.SetDeleteArray(&deleteArray_OnlineMonWindow);
      instance.SetDestructor(&destruct_OnlineMonWindow);
      instance.SetStreamerFunc(&streamer_OnlineMonWindow);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::OnlineMonWindow*)
   {
      return GenerateInitInstanceLocal((::OnlineMonWindow*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::OnlineMonWindow*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr OnlineMonWindow::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *OnlineMonWindow::Class_Name()
{
   return "OnlineMonWindow";
}

//______________________________________________________________________________
const char *OnlineMonWindow::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::OnlineMonWindow*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int OnlineMonWindow::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::OnlineMonWindow*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *OnlineMonWindow::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::OnlineMonWindow*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *OnlineMonWindow::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::OnlineMonWindow*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SimpleStandardHit(void *p) {
      delete ((::SimpleStandardHit*)p);
   }
   static void deleteArray_SimpleStandardHit(void *p) {
      delete [] ((::SimpleStandardHit*)p);
   }
   static void destruct_SimpleStandardHit(void *p) {
      typedef ::SimpleStandardHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimpleStandardHit

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimpleStandardCluster(void *p) {
      return  p ? new(p) ::SimpleStandardCluster : new ::SimpleStandardCluster;
   }
   static void *newArray_SimpleStandardCluster(Long_t nElements, void *p) {
      return p ? new(p) ::SimpleStandardCluster[nElements] : new ::SimpleStandardCluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimpleStandardCluster(void *p) {
      delete ((::SimpleStandardCluster*)p);
   }
   static void deleteArray_SimpleStandardCluster(void *p) {
      delete [] ((::SimpleStandardCluster*)p);
   }
   static void destruct_SimpleStandardCluster(void *p) {
      typedef ::SimpleStandardCluster current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimpleStandardCluster

namespace ROOT {
   // Wrappers around operator new
   static void *new_OnlineMonConfiguration(void *p) {
      return  p ? new(p) ::OnlineMonConfiguration : new ::OnlineMonConfiguration;
   }
   static void *newArray_OnlineMonConfiguration(Long_t nElements, void *p) {
      return p ? new(p) ::OnlineMonConfiguration[nElements] : new ::OnlineMonConfiguration[nElements];
   }
   // Wrapper around operator delete
   static void delete_OnlineMonConfiguration(void *p) {
      delete ((::OnlineMonConfiguration*)p);
   }
   static void deleteArray_OnlineMonConfiguration(void *p) {
      delete [] ((::OnlineMonConfiguration*)p);
   }
   static void destruct_OnlineMonConfiguration(void *p) {
      typedef ::OnlineMonConfiguration current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::OnlineMonConfiguration

namespace ROOT {
   // Wrapper around operator delete
   static void delete_SimpleStandardPlane(void *p) {
      delete ((::SimpleStandardPlane*)p);
   }
   static void deleteArray_SimpleStandardPlane(void *p) {
      delete [] ((::SimpleStandardPlane*)p);
   }
   static void destruct_SimpleStandardPlane(void *p) {
      typedef ::SimpleStandardPlane current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimpleStandardPlane

namespace ROOT {
   // Wrapper around operator delete
   static void delete_BaseCollection(void *p) {
      delete ((::BaseCollection*)p);
   }
   static void deleteArray_BaseCollection(void *p) {
      delete [] ((::BaseCollection*)p);
   }
   static void destruct_BaseCollection(void *p) {
      typedef ::BaseCollection current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BaseCollection

namespace ROOT {
   // Wrappers around operator new
   static void *new_CheckEOF(void *p) {
      return  p ? new(p) ::CheckEOF : new ::CheckEOF;
   }
   static void *newArray_CheckEOF(Long_t nElements, void *p) {
      return p ? new(p) ::CheckEOF[nElements] : new ::CheckEOF[nElements];
   }
   // Wrapper around operator delete
   static void delete_CheckEOF(void *p) {
      delete ((::CheckEOF*)p);
   }
   static void deleteArray_CheckEOF(void *p) {
      delete [] ((::CheckEOF*)p);
   }
   static void destruct_CheckEOF(void *p) {
      typedef ::CheckEOF current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CheckEOF

//______________________________________________________________________________
void OnlineMonWindow::Streamer(TBuffer &R__b)
{
   // Stream an object of class OnlineMonWindow.

   TGMainFrame::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_OnlineMonWindow(void *p) {
      delete ((::OnlineMonWindow*)p);
   }
   static void deleteArray_OnlineMonWindow(void *p) {
      delete [] ((::OnlineMonWindow*)p);
   }
   static void destruct_OnlineMonWindow(void *p) {
      typedef ::OnlineMonWindow current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_OnlineMonWindow(TBuffer &buf, void *obj) {
      ((::OnlineMonWindow*)obj)->::OnlineMonWindow::Streamer(buf);
   }
} // end of namespace ROOT for class ::OnlineMonWindow

namespace {
  void TriggerDictionaryInitialization_OnlineMonDict_Impl() {
    static const char* headers[] = {
"/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh",
"/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardHit.hh",
"/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardCluster.hh",
"/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardPlane.hh",
"/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonConfiguration.hh",
"/Users/gennai/DESY/euda53/monitors/onlinemon/include/CheckEOF.hh",
0
    };
    static const char* includePaths[] = {
"/Users/gennai/DESY/euda53/bd",
"/Users/gennai/DESY/euda53/./main/include",
"/Users/gennai/DESY/euda53/./extern/jsoncons-0.93/src",
"/Users/gennai/DESY/euda53/monitors/onlinemon/.",
"/Users/gennai/DESY/euda53/monitors/onlinemon/include",
"/Users/gennai/root/include",
"/Users/gennai/root/include",
"/Users/gennai/DESY/euda53/bd/monitors/onlinemon/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "OnlineMonDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$include/SimpleStandardHit.hh")))  __attribute__((annotate("$clingAutoload$/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh")))  SimpleStandardHit;
class __attribute__((annotate("$clingAutoload$include/SimpleStandardCluster.hh")))  __attribute__((annotate("$clingAutoload$/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh")))  SimpleStandardCluster;
class __attribute__((annotate("$clingAutoload$include/OnlineMonConfiguration.hh")))  __attribute__((annotate("$clingAutoload$/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh")))  OnlineMonConfiguration;
class __attribute__((annotate("$clingAutoload$include/SimpleStandardPlane.hh")))  __attribute__((annotate("$clingAutoload$/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh")))  SimpleStandardPlane;
class __attribute__((annotate("$clingAutoload$BaseCollection.hh")))  __attribute__((annotate("$clingAutoload$/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh")))  BaseCollection;
class __attribute__((annotate("$clingAutoload$CheckEOF.hh")))  __attribute__((annotate("$clingAutoload$/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh")))  CheckEOF;
class __attribute__((annotate("$clingAutoload$/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh")))  OnlineMonWindow;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "OnlineMonDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonWindow.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardHit.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardCluster.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/SimpleStandardPlane.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/OnlineMonConfiguration.hh"
#include "/Users/gennai/DESY/euda53/monitors/onlinemon/include/CheckEOF.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"BaseCollection", payloadCode, "@",
"CheckEOF", payloadCode, "@",
"OnlineMonConfiguration", payloadCode, "@",
"OnlineMonWindow", payloadCode, "@",
"SimpleStandardCluster", payloadCode, "@",
"SimpleStandardHit", payloadCode, "@",
"SimpleStandardPlane", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("OnlineMonDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_OnlineMonDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_OnlineMonDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_OnlineMonDict() {
  TriggerDictionaryInitialization_OnlineMonDict_Impl();
}
