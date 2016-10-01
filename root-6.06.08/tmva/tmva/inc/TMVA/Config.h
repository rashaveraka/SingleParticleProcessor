// @(#)root/tmva $Id$   
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Config                                                                *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      GLobal configuration settings (singleton)                                 *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker    <Andreas.Hocker@cern.ch>     - CERN, Switzerland       *
 *      Joerg Stelzer      <Joerg.Stelzer@cern.ch>      - CERN, Switzerland       *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch> - Iowa State U., USA      *
 *      Helge Voss         <Helge.Voss@cern.ch>         - MPI-K Heidelberg, GER   *
 *                                                                                *
 * Copyright (c) 2006:                                                            *
 *      CERN, Switzerland                                                         *
 *      Iowa State U., USA                                                        *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://mva.sourceforge.net/license.txt)                                       *
 **********************************************************************************/

#ifndef ROOT_TMVA_Config
#define ROOT_TMVA_Config

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Config                                                               //
//                                                                      //
// Singleton class for global configuration settings used by TMVA       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#if __cplusplus > 199711L
#include <atomic>
#endif
#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif

namespace TMVA {

   class MsgLogger;

   class Config {
               
   public:

      static Config& Instance();
      static void    DestroyInstance();

      Bool_t UseColor() const { return fUseColoredConsole; }
      void   SetUseColor( Bool_t uc ) { fUseColoredConsole = uc; }

      Bool_t IsSilent() const { return fSilent; }
      void   SetSilent( Bool_t s ) { fSilent = s; }

      Bool_t WriteOptionsReference() const { return fWriteOptionsReference; }
      void   SetWriteOptionsReference( Bool_t w ) { fWriteOptionsReference = w; }

      Bool_t DrawProgressBar() const { return fDrawProgressBar; }
      void   SetDrawProgressBar( Bool_t d ) { fDrawProgressBar = d; }

   public:

      class VariablePlotting;
      class IONames;

      VariablePlotting& GetVariablePlotting() { return fVariablePlotting; }
      IONames&          GetIONames()          { return fIONames; }

      // publicly accessible global settings
      class VariablePlotting {
         // data collection class to configure plotting of variables
      public:

         Float_t fTimesRMS;
         Int_t   fNbins1D;
         Int_t   fNbins2D;
         Int_t   fMaxNumOfAllowedVariablesForScatterPlots;
         Int_t   fNbinsMVAoutput;
         Int_t   fNbinsXOfROCCurve;
         Bool_t  fUsePaperStyle;

      } fVariablePlotting; // Customisable plotting properties

      // for file names and similar
      class IONames {

      public:

         TString fWeightFileDir;
         TString fWeightFileExtension;
         TString fOptionsReferenceFileDir;
      } fIONames; // Customisable weight file properties
         
      
   private:

      // private constructor
      Config();
      Config( const Config& );
      Config& operator=( const Config&);
      virtual ~Config();
#if __cplusplus > 199711L
      static std::atomic<Config*> fgConfigPtr;
#else
      static Config* fgConfigPtr;
#endif                  
   private:

#if __cplusplus > 199711L
      std::atomic<Bool_t> fUseColoredConsole;     // coloured standard output
      std::atomic<Bool_t> fSilent;                // no output at all
      std::atomic<Bool_t> fWriteOptionsReference; // if set true: Configurable objects write file with option reference
      std::atomic<Bool_t> fDrawProgressBar;       // draw progress bar to indicate training evolution
#else
      Bool_t fUseColoredConsole;     // coloured standard output
      Bool_t fSilent;                // no output at all
      Bool_t fWriteOptionsReference; // if set true: Configurable objects write file with option reference
      Bool_t fDrawProgressBar;       // draw progress bar to indicate training evolution
#endif
      mutable MsgLogger* fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }
         
      ClassDef(Config,0) // Singleton class for global configuration settings
   };

   // global accessor
   Config& gConfig();
}

#endif
