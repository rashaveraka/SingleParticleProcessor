
#include "HitResiduals.h"

#include <iostream>
#include <algorithm>    // std::sort

#include <EVENT/LCCollection.h>
#include <EVENT/LCObject.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/IMarlinTrack.h"

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include <marlin/AIDAProcessor.h>

#include "fpcompare.h"

#include "LinkDef.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <cmath>

//


using namespace std;
using namespace lcio ;
using namespace marlin ;

TH1F* SPAnalysisX ;
TCanvas* totalenergy = new TCanvas("total energy","total energy",800,800);
 
 
SPAnalysis aSPAnalysis ;


SPAnalysis::SPAnalysi() : Processor("HitResiduals") {

    // modify processor description
    _description = "HitResiduals plots te residual between the track fit and the hit in the local coordinate system u,v,w." ;


    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::TRACK,
			      "TrackCollectionName" , 
			      "Name of the input track collection"  ,
			      _inputTrkColName ,
			      std::string("HCalBarrelHits") ) ;
	
    registerProcessorParameter( "outFileName",
				"Name of the output root file",
				_outFileName,
				std::string("Singleparticle.root")
				);

    registerProcessorParameter( "treeName",
				"Name of the tree",
				_treeName,
				std::string("tree")
				);
/* Commenting out for development test - Ross
    registerProcessorParameter("MultipleScatteringOn",
			       "Use MultipleScattering in Fit",
			       _MSOn,
			       bool(true));
  
    registerProcessorParameter("EnergyLossOn",
			       "Use Energy Loss in Fit",
			       _ElossOn,
			       bool(true));
  
    registerProcessorParameter("SmoothOn",
			       "Smooth All Mesurement Sites in Fit",
			       _SmoothOn,
			       bool(false));
  
    registerProcessorParameter("MaxChi2Increment",
			       "Maximum increment allowed for the chi2",
			       _Max_Chi2_Incr,
			       double(1000.));
*/
}



void HitResiduals::init() { 

  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  _out  = new TFile(_outFileName.c_str(), "RECREATE");
  _tree = new TTree(_treeName.c_str(), _treeName.c_str());

  int bufsize = 32000; //default buffer size 32KB

  _tree->Branch("nRun", &_nRun,"nRun/I");
  _tree->Branch("nEvt", &_nEvt,"nEvt/I");
  _tree->Branch("resX", "std::vector<double >",&_resX,bufsize,0);
  _tree->Branch("resY", "std::vector<double >",&_resY,bufsize,0);
  _tree->Branch("subdet", "std::vector<int >",&_subdet,bufsize,0);
  _tree->Branch("layer", "std::vector<int >",&_layer,bufsize,0);

  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();

  const double pos[3]={0,0,0}; 
  double bFieldVec[3]={0,0,0}; 
  lcdd.field().magneticField(pos,bFieldVec); // get the magnetic field vector from DD4hep
  _bField = bFieldVec[2]/dd4hep::tesla; // z component at (0,0,0)

  // trksystem for marlin track
  _trksystem =  MarlinTrk::Factory::createMarlinTrkSystem( "DDKalTest" , marlin::Global::GEAR , "" ) ;
  
  if( _trksystem == 0 ) 
    throw EVENT::Exception( std::string("  Cannot initialize MarlinTrkSystem of Type: ") + std::string("DDKalTest") );
  
  _trksystem->setOption( IMarlinTrkSystem::CFG::useQMS,        _MSOn ) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::usedEdx,       _ElossOn) ;
  _trksystem->setOption( IMarlinTrkSystem::CFG::useSmoothing,  _SmoothOn) ;
  _trksystem->init() ;  
  
 } 


void HitResiduals::processRunHeader( LCRunHeader* run) { 

    _nRun++ ;
} 



void HitResiduals::processEvent( LCEvent * evt ) { 

  // this gets called for every event 
  std:cout << "Processing event " << _nEvt << std::endl;
  
  //streamlog_out(DEBUG2) << "----- _nEvt = " << _nEvt << std::endl;

  //clear vectors
  _resX.clear();
  _resY.clear();
  _subdet.clear();
  _layer.clear();

  UTIL::BitField64 cellid_decoder( SiDecoderString ) ; 	    
  UTIL::BitField64 encoder(SiDecoderString) ; 	    
  encoder.reset() ;  // reset to 0
  int layerID = encoder.lowWord() ;  
  int elementID = 0 ;    

  LCCollection* inputTrkCol = this->GetCollection( evt, _inputTrkColName );

  if (inputTrkCol!=0) {

    for( int i= 0; i < inputTrkCol->getNumberOfElements(); i++ ) {
      
      // std::cout << "Processing track " << i << std::endl;

      Track* track = dynamic_cast<Track*>( inputTrkCol->getElementAt(i) );
      MarlinTrk::IMarlinTrack* marlin_trk = _trksystem->createTrack();
      
      EVENT::TrackerHitVec trkHits = track->getTrackerHits() ;
      // sort(trkHits.begin(), trkHits.end(), RSort() );
      // reverse(trkHits.begin(), trkHits.end());
	
      for( EVENT::TrackerHitVec::iterator it = trkHits.begin() ; it != trkHits.end() ; ++it ){
      	marlin_trk->addHit(*it);	
      }
    
      int init_status = FitInitFromLCIOTrackState(track, marlin_trk);          
      
      double resX = 0.;
      double resY = 0.;
      
      for( EVENT::TrackerHitVec::iterator it = trkHits.begin(); it != trkHits.end(); ++it ){
	
	encoder.reset() ;  // reset to 0
	int layerID = encoder.lowWord() ;  	   
	int elementID = 0;
	DD4hep::long64 id = (*it)->getCellID0() ;
	cellid_decoder.setValue( id ) ;
	streamlog_out(DEBUG1) << "id = " << id << std::endl;

	int layer = cellid_decoder["layer"].value();
	int subdet = cellid_decoder["subdet"].value();
	streamlog_out(DEBUG1) << "layer = " << layer << std::endl;
	streamlog_out(DEBUG1) << "subdet = " << subdet << std::endl;

	//note: SiDecoderString="subdet:6,side:3,layer:4,module:12,sensor:1";
	encoder[lcio::ILDCellID0::subdet] = subdet;  
	
	encoder[lcio::ILDCellID0::layer]  = layer;   
	layerID = encoder.lowWord();  

	streamlog_out(DEBUG1) << "layerID = " << layerID << std::endl;

	TrackStateImpl trkState;
	double chi2 = 0 ;
	int ndf = 0 ;
	
	if ( marlin_trk->propagateToLayer( layerID, trkState, chi2, ndf, elementID, IMarlinTrack::modeClosest) == MarlinTrk::IMarlinTrack::success) {

	  // cout << "This is subdet number " << subdet << endl;
	  const float* pivot = trkState.getReferencePoint();
	  double fitX  = pivot[0];
	  double fitY  = pivot[1];
	  streamlog_out(DEBUG2) << " ----- fit x[mm], y[mm], r[mm] = " << fitX <<"   "<< fitY <<"   "<< sqrt(pow(fitX,2)+pow(fitY,2)) <<std::endl;
	  const double* hit_pos = (*it)->getPosition();
	  double hitX = hit_pos[0];
	  double hitY = hit_pos[1];
	  streamlog_out(DEBUG2) << " ----- hit x[mm], y[mm], r[mm] and layerID  =  " << hitX <<"   "<< hitY <<"   "<< sqrt(pow(hitX,2)+pow(hitY,2))<<"  "<<layerID<<std::endl;
	  
	  double resX = hitX - fitX ;
	  double resY = hitY -fitY ;

	  // cout<<"layerID  "<<layerID<<endl;
	  
    	  _resX.push_back(resX);
	  _resY.push_back(resY);
	  _subdet.push_back(subdet);
	  _layer.push_back(layer);

 	} else streamlog_out(DEBUG4) << "FAIL" << std::endl;

      }//end loop on hits

    }//end loop on tracks

  }

  _tree->Fill();
 
  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		       << "   in run:  " << evt->getRunNumber() << std::endl;

  _nEvt ++ ;

}



void HitResiduals::check( LCEvent * evt ) { 
    // nothing to check here 
}


void HitResiduals::end(){ 

  _tree->Write();
  _out->Write();
  _out->Close();
  
  std::cout << "HitResiduals::end()  " << name() 
    	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    	    << std::endl ;
}


LCCollection* HitResiduals::GetCollection( LCEvent*& evt, std::string colName ){

  LCCollection* col = NULL;  
  try{
    col = evt->getCollection( colName.c_str() ) ;
    streamlog_out( DEBUG3 ) << " --> " << colName.c_str() << " track collection found in event = " << col << " number of elements " << col->getNumberOfElements() << std::endl;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG3 ) << " --> " << colName.c_str() <<  " collection absent in event" << std::endl;     
  }
  return col; 
}

/*int HitResiduals::FitInit2( Track*& track, MarlinTrk::IMarlinTrack*& _marlinTrk ){
  
  TrackStateImpl trackState( TrackState::AtOther, 
			     track->getD0(), 
			     track->getPhi(), 
			     track->getOmega(), 
			     track->getZ0(), 
			     track->getTanLambda(), 
			     track->getCovMatrix(), 
			     track->getReferencePoint()
			     );
  
  _marlinTrk->initialise( trackState, _bField, IMarlinTrack::forward ) ;
  
  return IMarlinTrack::success ;   
  }*/


int HitResiduals::FitInitFromLCIOTrackState( Track*& track, MarlinTrk::IMarlinTrack*& _marlinTrk ){
  
  TrackStateImpl trackState( *(track->getTrackState(TrackState::AtFirstHit)) );  
  
  _marlinTrk->initialise( trackState, _bField, IMarlinTrack::forward ) ;
    
  return IMarlinTrack::success ;   
}
