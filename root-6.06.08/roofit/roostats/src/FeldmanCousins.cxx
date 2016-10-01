// @(#)root/roostats:$Id$
// Author: Kyle Cranmer   January 2009

/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef RooStats_FeldmanCousins
#include "RooStats/FeldmanCousins.h"
#endif

#ifndef RooStats_RooStatsUtils
#include "RooStats/RooStatsUtils.h"
#endif

#ifndef RooStats_PointSetInterval
#include "RooStats/PointSetInterval.h"
#endif

#include "RooStats/ModelConfig.h"

#include "RooStats/SamplingDistribution.h"

#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/NeymanConstruction.h"
#include "RooStats/RooStatsUtils.h"

#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGlobalFunc.h"
#include "RooMsgService.h"
#include "TFile.h"
#include "TTree.h"

ClassImp(RooStats::FeldmanCousins) ;

using namespace RooFit;
using namespace RooStats;
using namespace std;


/*
////////////////////////////////////////////////////////////////////////////////

FeldmanCousins::FeldmanCousins() : 
  //  fModel(NULL),
   fData(0),
   fTestStatSampler(0),
   fPointsToTest(0),
   fAdaptiveSampling(false), 
   fNbins(10), 
   fFluctuateData(true),
   fDoProfileConstruction(true),
   fSaveBeltToFile(false),
   fCreateBelt(false)
{
   // default constructor
//   fWS = new RooWorkspace("FeldmanCousinsWS");
//   fOwnsWorkspace = true;
//   fDataName = "";
//   fPdfName = "";
}
*/

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

FeldmanCousins::FeldmanCousins(RooAbsData& data, ModelConfig& model) : 
  fSize(0.05), 
  fModel(model),
  fData(data),
  fTestStatSampler(0),
  fPointsToTest(0),
  fPOIToTest(0),
  fConfBelt(0),
  fAdaptiveSampling(false), 
  fAdditionalNToysFactor(1.),
  fNbins(10), 
  fFluctuateData(true),
  fDoProfileConstruction(true),
  fSaveBeltToFile(false),
  fCreateBelt(false)
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor
///if(fOwnsWorkspace && fWS) delete fWS;

FeldmanCousins::~FeldmanCousins() {
  if(fPointsToTest) delete fPointsToTest;
  if(fPOIToTest) delete fPOIToTest;
  if(fTestStatSampler) delete fTestStatSampler;
}


////////////////////////////////////////////////////////////////////////////////
/// set the model

void FeldmanCousins::SetModel(const ModelConfig & model) { 
  fModel = model;
}

////////////////////////////////////////////////////////////////////////////////

TestStatSampler*  FeldmanCousins::GetTestStatSampler() const{
  if(!fTestStatSampler)
    this->CreateTestStatSampler();
  return fTestStatSampler; 
}

////////////////////////////////////////////////////////////////////////////////
/// specify the Test Statistic and create a ToyMC test statistic sampler

void FeldmanCousins::CreateTestStatSampler() const{
  // use the profile likelihood ratio as the test statistic
  ProfileLikelihoodTestStat* testStatistic = new ProfileLikelihoodTestStat(*fModel.GetPdf());
  
  // create the ToyMC test statistic sampler
  fTestStatSampler = new ToyMCSampler(*testStatistic,int(fAdditionalNToysFactor*50./fSize)) ;
  fTestStatSampler->SetParametersForTestStat(*fModel.GetParametersOfInterest() );
  if(fModel.GetObservables())
    fTestStatSampler->SetObservables(*fModel.GetObservables());
  fTestStatSampler->SetPdf(*fModel.GetPdf());
  
  if(!fAdaptiveSampling){
    ooccoutP(&fModel,Generation) << "FeldmanCousins: ntoys per point = " << (int) (fAdditionalNToysFactor*50./fSize) << endl;
  } else{
    ooccoutP(&fModel,Generation) << "FeldmanCousins: ntoys per point: adaptive" << endl;
  }
  if(fFluctuateData){
    ooccoutP(&fModel,Generation) << "FeldmanCousins: nEvents per toy will fluctuate about  expectation" << endl;
  } else{
    ooccoutP(&fModel,Generation) << "FeldmanCousins: nEvents per toy will not fluctuate, will always be " << fData.numEntries() << endl;
    fTestStatSampler->SetNEventsPerToy(fData.numEntries());
  }
}


////////////////////////////////////////////////////////////////////////////////
/// specify the parameter points to perform the construction.
/// allow ability to profile on some nuisance paramters

void FeldmanCousins::CreateParameterPoints() const{
  // get ingredients
  RooAbsPdf* pdf   = fModel.GetPdf(); 
  if (!pdf ){
    ooccoutE(&fModel,Generation) << "FeldmanCousins: ModelConfig has no PDF" << endl;
    return;
  }

  // get list of all paramters
  RooArgSet* parameters = new RooArgSet(*fModel.GetParametersOfInterest());
  if(fModel.GetNuisanceParameters())
    parameters->add(*fModel.GetNuisanceParameters());
  
  
  if( fModel.GetNuisanceParameters() && ! fModel.GetParametersOfInterest()->equals(*parameters) && fDoProfileConstruction) {
    // if parameters include nuisance parameters, do profile construction
    ooccoutP(&fModel,Generation) << "FeldmanCousins: Model has nuisance parameters, will do profile construction" << endl;
    
    // set nbins for the POI
    TIter it2 = fModel.GetParametersOfInterest()->createIterator();
    RooRealVar *myarg2; 
    while ((myarg2 = dynamic_cast<RooRealVar*>(it2.Next()))) { 
      myarg2->setBins(fNbins);
    }
    
    // get dataset for POI scan
    //     RooDataHist* parameterScan = NULL;
    RooAbsData* parameterScan = NULL;
    if(fPOIToTest)
      parameterScan = fPOIToTest;
    else
      parameterScan = new RooDataHist("parameterScan", "", *fModel.GetParametersOfInterest());


    ooccoutP(&fModel,Generation) << "FeldmanCousins: # points to test = " << parameterScan->numEntries() << endl;
    // make profile construction
    RooFit::MsgLevel previous  = RooMsgService::instance().globalKillBelow();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL) ;
    RooAbsReal* nll = pdf->createNLL(fData,RooFit::CloneData(false));
    RooAbsReal* profile = nll->createProfile(*fModel.GetParametersOfInterest());
    
    RooDataSet* profileConstructionPoints = new RooDataSet("profileConstruction",
							   "profileConstruction",
							   *parameters);
    
    
    for(Int_t i=0; i<parameterScan->numEntries(); ++i){
      // here's where we figure out the profiled value of nuisance parameters
      *parameters = *parameterScan->get(i);
      profile->getVal();
      profileConstructionPoints->add(*parameters);
    }   
    RooMsgService::instance().setGlobalKillBelow(previous) ;
    delete profile; 
    delete nll;
    if(!fPOIToTest) delete parameterScan;

    // done
    fPointsToTest = profileConstructionPoints;

  } else{
    // Do full construction
    ooccoutP(&fModel,Generation) << "FeldmanCousins: Model has no nuisance parameters" << endl;

    TIter it = parameters->createIterator();
    RooRealVar *myarg; 
    while ((myarg = dynamic_cast<RooRealVar*>(it.Next()))) { 
      myarg->setBins(fNbins);
    }

    RooDataHist* parameterScan = new RooDataHist("parameterScan", "", *parameters);
    ooccoutP(&fModel,Generation) << "FeldmanCousins: # points to test = " << parameterScan->numEntries() << endl;
    
    fPointsToTest = parameterScan;
  }
  
  delete parameters;
  
}


////////////////////////////////////////////////////////////////////////////////
/// Main interface to get a RooStats::ConfInterval.  
/// It constructs a RooStats::PointSetInterval.

PointSetInterval* FeldmanCousins::GetInterval() const {
  // local variables
  //  RooAbsData* data = fData; //fWS->data(fDataName);

  // fill in implied variables given data
  fModel.GuessObsAndNuisance(fData);
  
  // create the test statistic sampler (private data member fTestStatSampler)
  if(!fTestStatSampler)
    this->CreateTestStatSampler();

  fTestStatSampler->SetObservables(*fModel.GetObservables());

  if(!fFluctuateData)
    fTestStatSampler->SetNEventsPerToy(fData.numEntries());

  // create paramter points to perform construction (private data member fPointsToTest)
  this->CreateParameterPoints();

  // Create a Neyman Construction
  NeymanConstruction nc(fData,fModel);
  // configure it
  nc.SetTestStatSampler(*fTestStatSampler);
  nc.SetTestSize(fSize); // set size of test
  //  nc.SetParameters( fModel.GetParametersOfInterest);
  nc.SetParameterPointsToTest( *fPointsToTest );
  nc.SetLeftSideTailFraction(0.); // part of definition of Feldman-Cousins
  nc.SetData(fData);
  nc.UseAdaptiveSampling(fAdaptiveSampling);
  nc.AdditionalNToysFactor(fAdditionalNToysFactor);
  nc.SaveBeltToFile(fSaveBeltToFile);
  nc.CreateConfBelt(fCreateBelt);
  fConfBelt = nc.GetConfidenceBelt();
  // use it
  return nc.GetInterval();
}
