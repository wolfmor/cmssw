// -*- C++ -*-
//
// Package:    V0Producer
// Class:      V0Fitter
//
/**\class V0Fitter V0Fitter.cc RecoVertex/V0Producer/src/V0Fitter.cc

   Description: <one line class summary>

   Implementation:
   <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
//
//

#include "V0Fitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <typeinfo>
#include <memory>
#include <cmath>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include <TLorentzVector.h>

// pdg mass constants
namespace {
  const double piMass = 0.13957018;
  const double piMassSquared = piMass*piMass;
  const double protonMass = 0.938272046;
  const double protonMassSquared = protonMass*protonMass;
  const double kShortMass = 0.497614;
  const double lambdaMass = 1.115683;
}

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

V0Fitter::V0Fitter(const edm::ParameterSet& theParameters, edm::ConsumesCollector && iC)
{
  token_beamSpot = iC.consumes<reco::BeamSpot>(theParameters.getParameter<edm::InputTag>("beamSpot"));
  useVertex_ = theParameters.getParameter<bool>("useVertex");
  sameSign = theParameters.getParameter<bool>("sameSign");
  token_vertices = iC.consumes<std::vector<reco::Vertex>>(theParameters.getParameter<edm::InputTag>("vertices"));

  token_muons = iC.consumes<std::vector<reco::Muon>>(theParameters.getParameter<edm::InputTag>("muons"));
  token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  //label_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  vertexFitter_ = theParameters.getParameter<bool>("vertexFitter");
  useRefTracks_ = theParameters.getParameter<bool>("useRefTracks");
   
  // whether to reconstruct KShorts
  doKShorts_ = theParameters.getParameter<bool>("doKShorts");
  // whether to reconstruct Lambdas
  doLambdas_ = theParameters.getParameter<bool>("doLambdas");

  // cuts on initial track selection
  trackQualities_ = theParameters.getParameter<std::string>("trackQualities");
  tkChi2Cut_ = theParameters.getParameter<double>("tkChi2Cut");
  tkNHitsCut_ = theParameters.getParameter<int>("tkNHitsCut");
  tkPtCut_ = theParameters.getParameter<double>("tkPtCut");
  tkPtMax_ = theParameters.getParameter<double>("tkPtMax");
  tkIPSigXYCut_ = theParameters.getParameter<double>("tkIPSigXYCut");
  tkIPSigZCut_ = theParameters.getParameter<double>("tkIPSigZCut");
   
  // cuts on vertex
  vtxChi2Cut_ = theParameters.getParameter<double>("vtxChi2Cut");
  vtxDecaySigXYZCut_ = theParameters.getParameter<double>("vtxDecaySigXYZCut");
  vtxDecaySigXYCut_ = theParameters.getParameter<double>("vtxDecaySigXYCut");
  // miscellaneous cuts
  tkDCACut_ = theParameters.getParameter<double>("tkDCACut");
  mPiPiCut_ = theParameters.getParameter<double>("mPiPiCut");
  innerHitPosCut_ = theParameters.getParameter<double>("innerHitPosCut");
  cosThetaXYCut_ = theParameters.getParameter<double>("cosThetaXYCut");
  cosThetaXYZCut_ = theParameters.getParameter<double>("cosThetaXYZCut");
  // cuts on the V0 candidate mass
  kShortMassCut_ = theParameters.getParameter<double>("kShortMassCut");
  lambdaMassCut_ = theParameters.getParameter<double>("lambdaMassCut");
  
  //double_t nEvents = 0;
}

// method containing the algorithm for vertex reconstruction
void V0Fitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup,
		      reco::VertexCompositeCandidateCollection & theKshorts, std::vector<float> & dca_theKshorts , reco::VertexCompositeCandidateCollection & theLambdas, std::vector<float> & dca_theLambdas, std::vector<int> & iprotons )
{
  using std::vector;
  bool debug = false; 
  //nEvents ++;

  edm::Handle<reco::TrackCollection> theTrackHandle;
  iEvent.getByToken(token_tracks, theTrackHandle);
  //if (iEvent.eventAuxiliary().event() != 43) return;
  if (debug) std::cout << "------------New event ------------------------------------" <<  std::endl;
  if (debug) std::cout << "Entering the Fitter track muon...for a event with tracks: " << theTrackHandle->size() << std::endl;
  if (debug) std::cout << "--------------------- ------------------------------------" <<  std::endl;

  std::cout << "------------New event ------------------------------------" <<  std::endl;
  std::cout << "Entering the Fitter track muon...for a event with tracks: " << theTrackHandle->size() << std::endl;
  std::cout << "Entering the Fitter track muon...for event : " << iEvent.eventAuxiliary().event() << std::endl;
  std::cout << "--------------------- ------------------------------------" <<  std::endl;
    
  edm::Handle<std::vector<reco::Muon>> theMuonHandle;
  iEvent.getByToken(token_muons, theMuonHandle);
  
  if (theTrackHandle->empty()) return;
  const reco::TrackCollection* theTrackCollection = theTrackHandle.product();   
  //const reco::MuonCollection* theMuonCollection = theMuonHandle.product();   
  if (theTrackCollection->size()<2) return; // Quickfix should go into inputs at some time
  edm::Handle<reco::BeamSpot> theBeamSpotHandle;
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);
  const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
  math::XYZPoint referencePos(theBeamSpot->position());
  
  reco::Vertex referenceVtx;
  if (useVertex_) {
    edm::Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(token_vertices, vertices);
    referenceVtx = vertices->at(0);
    referencePos = referenceVtx.position();
  }

  edm::ESHandle<MagneticField> theMagneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagneticFieldHandle);
  const MagneticField* theMagneticField = theMagneticFieldHandle.product();

  std::vector<reco::TrackRef> theTrackRefs;
  std::vector<reco::TransientTrack> theTransTracks;

  std::vector<unsigned int> goodtrackindices;
  int itrack = -1;
  // fill vectors of TransientTracks and TrackRefs after applying preselection cuts
  for (reco::TrackCollection::const_iterator iTk = theTrackCollection->begin(); iTk != theTrackCollection->end(); ++iTk) {
    //
    itrack+=1;
    //if (itrack != 103 and itrack != 210) continue;
    //if (itrack != 94 and itrack != 95 and itrack != 51 and itrack != 1 and itrack != 0 and itrack != 2 and itrack != 125) continue;
    const reco::Track* tmpTrack = &(*iTk);
    if (false) std::cout << "Entering the track loop... for track " << itrack << " with pt " << tmpTrack->pt() << std::endl;
    double ipsigXY = std::abs(tmpTrack->dxy(*theBeamSpot)/tmpTrack->dxyError());
    if (false) std::cout << " ipsigXY is good..." << itrack << std::endl;
    if (useVertex_) ipsigXY = std::abs(tmpTrack->dxy(referencePos)/tmpTrack->dxyError()); // currently useVertex = false
    double ipsigZ = std::abs(tmpTrack->dz(referencePos)/tmpTrack->dzError());
    if (false) std::cout << " ipsigZ is good..." << itrack << std::endl;
    if (false) std::cout << " Checking: track Qual " << (tmpTrack->quality(tmpTrack->qualityByName(trackQualities_))) << 
    " trk Chi2 " << (tmpTrack->normalizedChi2() < tkChi2Cut_) << " trk Hits " << (tmpTrack->numberOfValidHits() >= tkNHitsCut_)<<
     " trk pt " << (tmpTrack->pt() > tkPtCut_) << " ipsigXY " << (ipsigXY > tkIPSigXYCut_) << " ipsigZ " << (ipsigZ > tkIPSigZCut_) << std::endl;
  
    if (tmpTrack->quality(tmpTrack->qualityByName(trackQualities_)) &&
	tmpTrack->normalizedChi2() < tkChi2Cut_ &&
	tmpTrack->numberOfValidHits() >= tkNHitsCut_ &&
	tmpTrack->pt() > tkPtCut_ && tmpTrack->pt() < tkPtMax_ && ipsigXY > tkIPSigXYCut_ && ipsigZ > tkIPSigZCut_) {
	  if (false) std::cout << "pass preselection.." << itrack << std::endl;
      reco::TrackRef tmpRef(theTrackHandle, std::distance(theTrackCollection->begin(), iTk));
      theTrackRefs.push_back(std::move(tmpRef));
      reco::TransientTrack tmpTransient(*tmpRef, theMagneticField);
      theTransTracks.push_back(std::move(tmpTransient));
      goodtrackindices.push_back(itrack);
    }
  }
  //std::cout << theTransTracks.size() << " good tracks have now been selected for vertexing..." <<  std::endl;
  // good tracks have now been selected for vertexing

  // loop over tracks and vertex good charged track pairs
  protonidx = -1;
  int positiveIdx = -1;
  int negativeIdx = -1;
  int posCharge = 999;
  int negCharge = 999;
  if (debug) std::cout << "----------------Pairing of preselected tracks ----------------"<< std::endl;
  for (unsigned int trdx1 = 0; trdx1 < theTrackRefs.size(); ++trdx1) {
    //if ((trdx1 != 197)) continue;
    if (false) std::cout << "Take track " << trdx1 << std::endl;
    for (unsigned int trdx2 = trdx1 + 1; trdx2 < theTrackRefs.size(); ++trdx2) {
      //if ((trdx2 != 56)) continue;
      reco::TrackRef positiveTrackRef;
      reco::TrackRef negativeTrackRef;
      reco::TransientTrack* posTransTkPtr = nullptr;
      reco::TransientTrack* negTransTkPtr = nullptr;
      if (false) std::cout << "Take track " << trdx1 << std::endl;
      if (false) std::cout << "Pair it with track " << trdx2 <<  "has pt " << theTrackRefs[trdx2]->pt() << " and charge " << theTrackRefs[trdx2]->charge() << std::endl;
      if (false) std::cout << "have opposite charge " << ((theTrackRefs[trdx1]->charge() < 0. && theTrackRefs[trdx2]->charge() > 0.) or (theTrackRefs[trdx1]->charge() > 0. && theTrackRefs[trdx2]->charge() < 0.))<<  std::endl;
	  if (sameSign){
		if (theTrackRefs[trdx1]->charge() < 0. && theTrackRefs[trdx2]->charge() > 0.)  {
			continue;
		  } else if (theTrackRefs[trdx1]->charge() > 0. && theTrackRefs[trdx2]->charge() < 0.)  {
			continue;
		  } else {
			negativeTrackRef = theTrackRefs[trdx2];
			positiveTrackRef = theTrackRefs[trdx1];
			negTransTkPtr = &theTransTracks[trdx2];
			posTransTkPtr = &theTransTracks[trdx1];
			positiveIdx = goodtrackindices.at(trdx1);
			negativeIdx = goodtrackindices.at(trdx2);
		  }

		if (false) std::cout << "Setting the pos. and neg. TrackRef " << positiveTrackRef->charge() << " " << negativeTrackRef->charge()  << std::endl;
	  }
	  else{
		  if (theTrackRefs[trdx1]->charge() < 0. && theTrackRefs[trdx2]->charge() > 0.) {
		negativeTrackRef = theTrackRefs[trdx1];
		positiveTrackRef = theTrackRefs[trdx2];
		negTransTkPtr = &theTransTracks[trdx1];
		posTransTkPtr = &theTransTracks[trdx2];	
		positiveIdx = goodtrackindices.at(trdx2);
		negativeIdx = goodtrackindices.at(trdx1);
		  } else if (theTrackRefs[trdx1]->charge() > 0. && theTrackRefs[trdx2]->charge() < 0.) {
		negativeTrackRef = theTrackRefs[trdx2];
		positiveTrackRef = theTrackRefs[trdx1];
		negTransTkPtr = &theTransTracks[trdx2];
		posTransTkPtr = &theTransTracks[trdx1];
		positiveIdx = goodtrackindices.at(trdx1);
		negativeIdx = goodtrackindices.at(trdx2);
		  } else {
		continue;
		  }
	}
	posCharge = positiveTrackRef->charge();
	negCharge = negativeTrackRef->charge();
      
	// put muon loop here
	TLorentzVector tlv1;
	TLorentzVector tlv2;
	tlv1.SetPxPyPzE(theTrackRefs[trdx1]->px(), theTrackRefs[trdx1]->py(),theTrackRefs[trdx1]->pz(),theTrackRefs[trdx1]->pt()*cosh(theTrackRefs[trdx1]->eta()));
	tlv2.SetPxPyPzE(theTrackRefs[trdx2]->px(), theTrackRefs[trdx2]->py(),theTrackRefs[trdx2]->pz(),theTrackRefs[trdx2]->pt()*cosh(theTrackRefs[trdx2]->eta()));
	//int imuon = -1;
	
	double_t dRmin1 = 9999;	
	double_t dRmin2 = 9999;	
	for( const auto& muon : *theMuonHandle){
		TLorentzVector tlvmu;
		if (false) std::cout << " in muon matching, muon has charge" << muon.charge() << std::endl;
		if (muon.charge() == 0) continue;
		//if (muon.charge()*theTrackRefs[trdx1]->charge()>0) continue;
		if (muon.charge()*theTrackRefs[trdx1]->charge()>0){
			tlvmu.SetPxPyPzE(muon.px(), muon.py(),muon.pz(),muon.pt()*cosh(muon.eta()));
			if (tlvmu.DeltaR(tlv1)<dRmin1){
				dRmin1 = tlvmu.DeltaR(tlv1);
			}
		}
		if (muon.charge()*theTrackRefs[trdx2]->charge()>0){
			tlvmu.SetPxPyPzE(muon.px(), muon.py(),muon.pz(),muon.pt()*cosh(muon.eta()));
			if (tlvmu.DeltaR(tlv2)<dRmin2){
				dRmin2 = tlvmu.DeltaR(tlv2);
			}
		}
	}

	if  (((theTrackRefs[trdx1]->charge() < 0. && theTrackRefs[trdx2]->charge() < 0.)) or ((theTrackRefs[trdx1]->charge() > 0. && theTrackRefs[trdx2]->charge() > 0.))){
		if (false) std::cout << "hAVE SAME chARGE " << std::endl;
		//std::cout << "hAVE SAME chARGE " << trdx1 << std::endl;
		//std::cout << "Pair it with track " << trdx2 <<  "has pt " << theTrackRefs[trdx2]->pt() << " and charge " << theTrackRefs[trdx2]->charge() << std::endl;
		//std::cout << "muon matching drmin1 " << dRmin1 << " dRmin2 " << dRmin2 << std::endl;
	}
    if (false) std::cout << "muon matching drmin1 " << dRmin1 << " dRmin2 " << dRmin2 << std::endl;
    if ((dRmin1>0.01) && (dRmin2>0.01)) continue;
  	if (debug) std::cout << "Take track 1" << trdx1 <<  "has pt " << theTrackRefs[trdx1]->pt() << " and charge " << theTrackRefs[trdx1]->charge() << std::endl;
	if (debug) std::cout << "Pair it with track 2 " << trdx2 <<  "has pt " << theTrackRefs[trdx2]->pt() << " and charge " << theTrackRefs[trdx2]->charge() << std::endl;
	if (debug) std::cout << "muon matching drmin1 " << dRmin1 << " dRmin2 " << dRmin2 << std::endl;  
    
    
    
      // measure distance between tracks at their closest approach


      //these two variables are needed to 'pin' the temporary value returned to the stack
      // in order to keep posState and negState from pointing to destructed objects
      auto const& posImpact = posTransTkPtr->impactPointTSCP();

      

      //std::cout << "found this nifty lil guy " << posImpact.pt() << std::endl;//ok this lit up 
      auto const& negImpact = negTransTkPtr->impactPointTSCP();
      if (!posImpact.isValid() || !negImpact.isValid()) continue;
      if (debug) std::cout << "1. good tracks have valid impactPointTSCP ..." <<  std::endl;
      FreeTrajectoryState const & posState = posImpact.theState(); // 6-dimensional state vector of a helix
      FreeTrajectoryState const & negState = negImpact.theState();
      

      //std::cout << " 6-dimensional state vector of a helix " << posState.momentum() << negState.momentum() <<  std::endl;
      ClosestApproachInRPhi cApp;
      cApp.calculate(posState, negState);
      //std::cout << "ClosestApproachInRPhi:: status"<< cApp.status() <<  std::endl; // cApp.calculate(posState, negState) should be 1 (if == 0  calculation fails (e.g. concentric circles), return value = 0;)
      //std::cout << "Given two trajectory states, compute the two points of closest approach. Calculation yields..."<< cApp.calculate(posState, negState) << " Dpcas" << cApp.distance()<< std::endl; // cApp.calculate(posState, negState) should be 1 (if == 0  calculation fails (e.g. concentric circles), return value = 0;)
      
      //std::cout << "status of this: " << cApp.status() << std::endl;

      if (!cApp.status()) continue;
      
      

      float dca = std::abs(cApp.distance());
	  if (debug) std::cout << "3. dca <  tkDCACut_...; dca = " << dca <<  std::endl;
      if (dca > tkDCACut_) continue;
      if (debug) std::cout << "3. dca <  tkDCACut_...; dca = " << dca <<  std::endl;
      //std::cout << "dca " << dca <<  std::endl;
      // the POCA should at least be in the sensitive volume
      GlobalPoint cxPt = cApp.crossingPoint();
      if (sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) continue;
      if (debug) std::cout << "4. POCA is at least in the sensitive volume ..." <<  std::endl;
	  
      // the tracks should at least point in the same quadrant
      TrajectoryStateClosestToPoint const & posTSCP = posTransTkPtr->trajectoryStateClosestToPoint(cxPt);
      TrajectoryStateClosestToPoint const & negTSCP = negTransTkPtr->trajectoryStateClosestToPoint(cxPt);
      if (!posTSCP.isValid() || !negTSCP.isValid()) continue; 
      if (debug) std::cout << "5. the tracks trajectories at poca are valid..." <<  std::endl;
      //if (posTSCP.momentum().dot(negTSCP.momentum())  < 0) continue;
      if (debug) std::cout << "6. the tracks at least point in the same quadrant..." <<  std::endl;
      
      // calculate mPiPi
      double totalE = sqrt(posTSCP.momentum().mag2() + piMassSquared) + sqrt(negTSCP.momentum().mag2() + piMassSquared);
      double totalESq = totalE*totalE;
      double totalPSq = (posTSCP.momentum() + negTSCP.momentum()).mag2();
      double mass = sqrt(totalESq - totalPSq);
      if (mass > mPiPiCut_) continue;
      if (debug) std::cout << "7. passed  mPiPi cut..." <<  std::endl;
      // Fill the vector of TransientTracks to send to KVF
      std::vector<reco::TransientTrack> transTracks;
      transTracks.reserve(2);
      transTracks.push_back(*posTransTkPtr);
      transTracks.push_back(*negTransTkPtr);

      // create the vertex fitter object and vertex the tracks
      TransientVertex theRecoVertex;
      if (vertexFitter_) {
	KalmanVertexFitter theKalmanFitter(useRefTracks_ == 0 ? false : true);
	theRecoVertex = theKalmanFitter.vertex(transTracks);
      } else if (!vertexFitter_) {
	useRefTracks_ = false;
	AdaptiveVertexFitter theAdaptiveFitter;
	theRecoVertex = theAdaptiveFitter.vertex(transTracks);
      }
      if (!theRecoVertex.isValid()) continue;
      if (debug) std::cout << "8. created valid vertex..." <<  std::endl;
      reco::Vertex theVtx = theRecoVertex;
      if (debug) std::cout << " theVtx.normalizedChi2():" << theVtx.normalizedChi2() <<  std::endl;
      if (theVtx.normalizedChi2() > vtxChi2Cut_) continue;
      if (debug) std::cout << "9. passed  chi2 cut..." <<  std::endl;
      GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

      // 2D decay significance
      SMatrixSym3D totalCov = theBeamSpot->rotatedCovariance3D() + theVtx.covariance();
      if (useVertex_) totalCov = referenceVtx.covariance() + theVtx.covariance();
      SVector3 distVecXY(vtxPos.x()-referencePos.x(), vtxPos.y()-referencePos.y(), 0.);
      double distMagXY = ROOT::Math::Mag(distVecXY);
      double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
      if (distMagXY/sigmaDistMagXY < vtxDecaySigXYCut_) continue;
      if (debug) std::cout << "10. passed  2dsig cut..." <<  std::endl;
	  
      // 3D decay significance
      if (vtxDecaySigXYZCut_ > 0.) {
	SVector3 distVecXYZ(vtxPos.x()-referencePos.x(), vtxPos.y()-referencePos.y(), vtxPos.z()-referencePos.z());
	double distMagXYZ = ROOT::Math::Mag(distVecXYZ);
	double sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;
	if (distMagXYZ/sigmaDistMagXYZ < vtxDecaySigXYZCut_) continue;
	if (debug) std::cout << "11. passed  23dsig cut..." <<  std::endl;
      }

      // make sure the vertex radius is within the inner track hit radius
      if (innerHitPosCut_ > 0. && positiveTrackRef->innerOk()) {
	reco::Vertex::Point posTkHitPos = positiveTrackRef->innerPosition();
	double posTkHitPosD2 =  (posTkHitPos.x()-referencePos.x())*(posTkHitPos.x()-referencePos.x()) +
	  (posTkHitPos.y()-referencePos.y())*(posTkHitPos.y()-referencePos.y());
	if (sqrt(posTkHitPosD2) < (distMagXY - sigmaDistMagXY*innerHitPosCut_)) continue;
	if (debug) std::cout << "12. make sure the vertex radius is within the inner track hit radius..." <<  std::endl;
      }
      if (innerHitPosCut_ > 0. && negativeTrackRef->innerOk()) {
	reco::Vertex::Point negTkHitPos = negativeTrackRef->innerPosition();
	double negTkHitPosD2 = (negTkHitPos.x()-referencePos.x())*(negTkHitPos.x()-referencePos.x()) +
	  (negTkHitPos.y()-referencePos.y())*(negTkHitPos.y()-referencePos.y());
	if (sqrt(negTkHitPosD2) < (distMagXY - sigmaDistMagXY*innerHitPosCut_)) continue;
	if (debug) std::cout << "13. make sure the vertex radius is within the inner track hit radius..." <<  std::endl;
      }
      
      std::unique_ptr<TrajectoryStateClosestToPoint> trajPlus;
      std::unique_ptr<TrajectoryStateClosestToPoint> trajMins;
      std::vector<reco::TransientTrack> theRefTracks;
      if (theRecoVertex.hasRefittedTracks()) {
	theRefTracks = theRecoVertex.refittedTracks();
      }

      if (useRefTracks_ && theRefTracks.size() > 1) {
	reco::TransientTrack* thePositiveRefTrack = nullptr;
	reco::TransientTrack* theNegativeRefTrack = nullptr;
	for (std::vector<reco::TransientTrack>::iterator iTrack = theRefTracks.begin(); iTrack != theRefTracks.end(); ++iTrack) {
	  if (iTrack->track().charge() > 0.) {
	    thePositiveRefTrack = &*iTrack;
	  } else if (iTrack->track().charge() < 0.) {
	    theNegativeRefTrack = &*iTrack;
	  }
	}
	if(sameSign){
		//if (theRefTracks.size() > 2) continue;
		for (std::vector<reco::TransientTrack>::iterator iTrack = theRefTracks.begin(); iTrack != theRefTracks.end(); ++iTrack) {
			if (iTrack== theRefTracks.begin()){thePositiveRefTrack = &*iTrack;}
			else {theNegativeRefTrack = &*iTrack;}
		}
	}
	if (thePositiveRefTrack == nullptr || theNegativeRefTrack == nullptr) continue;
	trajPlus.reset(new TrajectoryStateClosestToPoint(thePositiveRefTrack->trajectoryStateClosestToPoint(vtxPos)));
	trajMins.reset(new TrajectoryStateClosestToPoint(theNegativeRefTrack->trajectoryStateClosestToPoint(vtxPos)));
      } else {
	trajPlus.reset(new TrajectoryStateClosestToPoint(posTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
	trajMins.reset(new TrajectoryStateClosestToPoint(negTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
      }

      if (trajPlus.get() == nullptr || trajMins.get() == nullptr || !trajPlus->isValid() || !trajMins->isValid()) continue;
      if (debug) std::cout << "14. make sure the vertex radius is within the inner track hit radius..." <<  std::endl;
      GlobalVector positiveP(trajPlus->momentum());
      GlobalVector negativeP(trajMins->momentum());
      GlobalVector totalP(positiveP + negativeP);

      // 2D pointing angle
      double dx = theVtx.x()-referencePos.x();
      double dy = theVtx.y()-referencePos.y();
      double px = totalP.x();
      double py = totalP.y();
      double angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
      if (angleXY < cosThetaXYCut_) continue;
      if (debug) std::cout << "15. pass 2D pointing angle cut..." <<  std::endl;
      // 3D pointing angle
      if (cosThetaXYZCut_ > -1.) {
	double dz = theVtx.z()-referencePos.z();
	double pz = totalP.z();
	double angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
	if (angleXYZ < cosThetaXYZCut_) continue;
	if (debug) std::cout << "16. cosThetaXYZCut_ > -1: pass 3D pointing angle cut..." <<  std::endl;
      }

      // calculate total energy of V0 3 ways: assume it's a kShort, a Lambda, or a LambdaBar.
      double piPlusE = sqrt(positiveP.mag2() + piMassSquared);
      double piMinusE = sqrt(negativeP.mag2() + piMassSquared);
      double protonE = sqrt(positiveP.mag2() + protonMassSquared);
      double antiProtonE = sqrt(negativeP.mag2() + protonMassSquared);
      double kShortETot = piPlusE + piMinusE;
      double lambdaEtot = protonE + piMinusE;
      double lambdaBarEtot = antiProtonE + piPlusE;

      // Create momentum 4-vectors for the 3 candidate types
      const reco::Particle::LorentzVector kShortP4(totalP.x(), totalP.y(), totalP.z(), kShortETot);
      const reco::Particle::LorentzVector lambdaP4(totalP.x(), totalP.y(), totalP.z(), lambdaEtot);
      const reco::Particle::LorentzVector lambdaBarP4(totalP.x(), totalP.y(), totalP.z(), lambdaBarEtot);

      reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
      const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
      double vtxChi2(theVtx.chi2());
      double vtxNdof(theVtx.ndof());

      // Create the VertexCompositeCandidate object that will be stored in the Event
      reco::VertexCompositeCandidate* theKshort = nullptr;
      reco::VertexCompositeCandidate* theLambda = nullptr;
      reco::VertexCompositeCandidate* theLambdaBar = nullptr;

      if (doKShorts_) {
	theKshort = new reco::VertexCompositeCandidate(0, kShortP4, vtx, vtxCov, vtxChi2, vtxNdof);
      }
      if (doLambdas_) {
	if (positiveP.mag2() > negativeP.mag2()) {
	  theLambda = new reco::VertexCompositeCandidate(0, lambdaP4, vtx, vtxCov, vtxChi2, vtxNdof);
	} else {
	  theLambdaBar = new reco::VertexCompositeCandidate(0, lambdaBarP4, vtx, vtxCov, vtxChi2, vtxNdof);
	}
      }

      // Create daughter candidates for the VertexCompositeCandidates
      
	  reco::RecoChargedCandidate thePiPlusCand(
					       posCharge, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), piPlusE), vtx);
      thePiPlusCand.setTrack(positiveTrackRef);
      if (debug) std::cout << "                    Built thePiPlusCand PiMinusCand... with charge" << posCharge <<  " " << negCharge << std::endl;
      
      reco::RecoChargedCandidate thePiMinusCand(
						negCharge, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), piMinusE), vtx);
      thePiMinusCand.setTrack(negativeTrackRef);
      
      reco::RecoChargedCandidate theProtonCand(
					       1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), protonE), vtx);
      theProtonCand.setTrack(positiveTrackRef);

      reco::RecoChargedCandidate theAntiProtonCand(
						   -1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), antiProtonE), vtx);
      theAntiProtonCand.setTrack(negativeTrackRef);
      
      //if(sameSign){
		  
		//std::cout << "                    Build thePiPlusCand PiMinusCand... with charge" << charge <<  std::endl;
	  //reco::RecoChargedCandidate thePiPlusCand(
						   //1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), piPlusE), vtx);
	  //thePiPlusCand.setTrack(positiveTrackRef);
	  
	  //reco::RecoChargedCandidate thePiMinusCand(
						//1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), piMinusE), vtx);
	  //thePiMinusCand.setTrack(negativeTrackRef);
	  
      //reco::RecoChargedCandidate theProtonCand(
					       //1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), protonE), vtx);
      //theProtonCand.setTrack(positiveTrackRef);

      //reco::RecoChargedCandidate theAntiProtonCand(
						   //-1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), antiProtonE), vtx);
      //theAntiProtonCand.setTrack(negativeTrackRef);

		//}
	//else{
		//reco::RecoChargedCandidate thePiPlusCand(
					       //1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), piPlusE), vtx);
      //thePiPlusCand.setTrack(positiveTrackRef);
      
      //reco::RecoChargedCandidate thePiMinusCand(
						//-1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), piMinusE), vtx);
      //thePiMinusCand.setTrack(negativeTrackRef);
      
      //reco::RecoChargedCandidate theProtonCand(
					       //1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), protonE), vtx);
      //theProtonCand.setTrack(positiveTrackRef);

      //reco::RecoChargedCandidate theAntiProtonCand(
						   //-1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), antiProtonE), vtx);
      //theAntiProtonCand.setTrack(negativeTrackRef);
		
		//}


      AddFourMomenta addp4;
      // Store the daughter Candidates in the VertexCompositeCandidates if they pass mass cuts
      if (doKShorts_) {
	theKshort->addDaughter(thePiPlusCand);
	if (debug) std::cout << "                    Adding thePiPlusCand daughter... with charge" << (theKshort->daughter(0))->charge() <<  std::endl;
	//theKshort->addDaughter(thePiMinusCand);
	theKshort->addDaughter(thePiMinusCand);
	if (debug) std::cout << "                    Adding thePiMinusCand daughter... with charge" << (theKshort->daughter(1))->charge() <<  std::endl;
	theKshort->setPdgId(310);
	addp4.set(*theKshort);
	if (debug) std::cout << "17.doing kshorts..." <<  std::endl;
	if (theKshort->mass() < kShortMass + kShortMassCut_ && theKshort->mass() > kShortMass - kShortMassCut_) {
	  if (debug) std::cout << "18. passed kShorMassCut..." <<  std::endl;
	  if (debug) std::cout << "18. passed kShorMassCut..." <<  std::endl;
	  theKshorts.push_back(std::move(*theKshort));
	  dca_theKshorts.push_back(std::move(dca));
	}
      }
      if (doLambdas_ && theLambda) {
	theLambda->addDaughter(theProtonCand);
	theLambda->addDaughter(thePiMinusCand);
	theLambda->setPdgId(3122);
	addp4.set( *theLambda );
	if (theLambda->mass() < lambdaMass + lambdaMassCut_ && theLambda->mass() > lambdaMass - lambdaMassCut_) {
	  theLambdas.push_back(std::move(*theLambda));
	  iprotons.push_back(std::move(positiveIdx));
	  dca_theLambdas.push_back(std::move(dca));
	}
      } else if (doLambdas_ && theLambdaBar) {
	theLambdaBar->addDaughter(theAntiProtonCand);
	theLambdaBar->addDaughter(thePiPlusCand);
	theLambdaBar->setPdgId(-3122);
	addp4.set(*theLambdaBar);
	if (theLambdaBar->mass() < lambdaMass + lambdaMassCut_ && theLambdaBar->mass() > lambdaMass - lambdaMassCut_) {
	  theLambdas.push_back(std::move(*theLambdaBar));
	  iprotons.push_back(std::move(negativeIdx));
	  dca_theLambdas.push_back(std::move(dca));
	}
      }

      delete theKshort;
      delete theLambda;
      delete theLambdaBar;
      theKshort = theLambda = theLambdaBar = nullptr;

    }
  }
}

