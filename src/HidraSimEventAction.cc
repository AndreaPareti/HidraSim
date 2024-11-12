//**************************************************
// \file HidraSimEventAction.cc
// \brief: Implementation of HidraSimEventAction 
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "HidraSimEventAction.hh"
#include "HidraSimRunAction.hh"
#include "HidraSimDetectorConstruction.hh"
#include "G4HCofThisEvent.hh"
//#include "HidraSimGeoPar.hh"

//includers for calorimeter hits
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "HidraSimCalorimeterHit.hh"
#include "HidraSimCalorimeterSD.hh"
#include "G4THitsMap.hh"

//Includers from Geant4
//
#include "g4root.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

//Includers from C++
//
#include <iomanip>
#include <vector>
#include<numeric>

// Include SiPM EndOfEvent action
//#include "sipm/SiPMProperties.h"
//#include "sipm/SiPMAnalogSignal.h"




//Define constructor
//
HidraSimEventAction::HidraSimEventAction()
    : G4UserEventAction(),
    EnergyScin(0.),
    EnergyCher(0.),
    NofCherDet(0),
    NofScinDet(0),
    EnergyTot(0.),
    PrimaryPDGID(0),
    PrimaryParticleEnergy(0.),
    PrimaryX(0),
    PrimaryY(0),
    EscapedEnergy(0.),
    EscapedEnergyl(0.),
    EscapedEnergyd(0.),
    PSEnergy(0.),
    VectorSignals(0.),
    VectorSignalsCher(0.),
    VecSPMT(0.),
    VecCPMT(0.),
    VecTowerE(0.),
    VecLeakCounter(0.),
    fSfiberHCID(-1),
    fCfiberHCID(-1)
    {
}

//Define de-constructor
//
HidraSimEventAction::~HidraSimEventAction() {}



HidraSimCalorimeterHitsCollection* 
HidraSimEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<HidraSimCalorimeterHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("HidraSimEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    


//Define BeginOfEventAction() and EndOfEventAction() methods
//
void HidraSimEventAction::BeginOfEventAction(const G4Event*) {  
    
    //Initialize data memebers at begin of each event
    //
    EnergyScin = 0.;
    EnergyCher = 0.;
    NofCherDet = 0;
    NofScinDet = 0;
    EnergyTot = 0.;
    PrimaryPDGID = 0;
    PrimaryX = 0;
    PrimaryY = 0;
    PrimaryParticleEnergy = 0.;
    EscapedEnergy = 0.;
    EscapedEnergyl = 0.;
    EscapedEnergyd = 0.;
    PSEnergy = 0.;

    VectorSignals.clear();
    VectorSignalsCher.clear();
    VecSPMT.clear();
    VecCPMT.clear();
    VecTowerE.clear();
    VecLeakCounter.clear();
    fHitPheSvector.clear();
    fHitZcoordSvector.clear();
    fHitSiPMIDSvector.clear();
    fHitPheCvector.clear();
    fHitZcoordCvector.clear();
    fHitSiPMIDCvector.clear();


    VectorSignals.assign(NoFibersTower*NofModulesSiPM, 0.);
    VectorSignalsCher.assign(NoFibersTower*NofModulesSiPM, 0.);
    VecSPMT.assign(NoModulesActive, 0.);
    VecCPMT.assign(NoModulesActive, 0.);
    VecTowerE.assign(NoModulesActive, 0.);
    VecLeakCounter.assign(4*NofLeakCounterLayers+1, 0.);
    fHitPheSvector.assign(fHitPheSvector.size(), 0.);
    fHitZcoordSvector.assign(fHitZcoordSvector.size(), 0.);
    fHitSiPMIDSvector.assign(fHitSiPMIDSvector.size(), 0.);
    fHitPheSvector.assign(fHitPheCvector.size(), 0.);
    fHitZcoordSvector.assign(fHitZcoordCvector.size(), 0.);
    fHitSiPMIDSvector.assign(fHitSiPMIDCvector.size(), 0.);
}









void HidraSimEventAction::EndOfEventAction(const G4Event* event) {
 
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Get hits collections IDs (only once)
  if ( fSfiberHCID == -1 ) {
    fSfiberHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("SfiberHitsCollection");
    fCfiberHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("CfiberHitsCollection");
  }

  //get hits collections IDs
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  //if (!hce) return;

  if (hce) {
        // Get the number of hit collections
        G4int numberOfHitCollections = hce->GetNumberOfCollections();
        G4cout << "Number of hit collections: " << numberOfHitCollections << G4endl;

        // Optionally, loop over hit collections and print their names
        //for (G4int i = 0; i < numberOfHitCollections; i++) {
        //    G4VHitsCollection* hitsCollection = hce->GetHC(i);
        //    if (hitsCollection) {
        //        G4cout << "Hit collection " << i << " name: " << hitsCollection->GetName() << G4endl;
        //    }
        //}
    } else {
       G4cout << "No hit collections in this event." << G4endl;
    }



  G4double rindexScin = 1.59;
  G4double rindexCher = 1.49;
  G4double c = 0.2998; // speed of light in mm / ps
  G4double vS = c/rindexScin;   // photon velocity in S fibers
  G4double vC = c/rindexCher;

  // Get hits collections
  auto SfiberHC = GetHitsCollection(fSfiberHCID, event);
  auto CfiberHC = GetHitsCollection(fCfiberHCID, event);
  int SnOfHitsCollections = SfiberHC->entries();
  int CnOfHitsCollections = CfiberHC->entries();

  std::vector<double> HitPheSvector, HitZcoordSvector, HitSiPMIDSvector;
  std::vector<double> HitPheCvector, HitZcoordCvector, HitSiPMIDCvector;

  // Get individual hits in S fibers
  for(int current_Shit=0; current_Shit<SnOfHitsCollections; current_Shit++)
  {
    auto ShitCollection = (*SfiberHC)[current_Shit];
    fHitPheSvector.push_back(ShitCollection->GetPhe());
    fHitZcoordSvector.push_back(ShitCollection->GetZcoord());
    fHitSiPMIDSvector.push_back(ShitCollection->GetSiPMID());
    
    std::vector<G4double> fiberPheVec = ShitCollection->GetPheVec();
    std::vector<G4double> fiberZVec = ShitCollection->GetZVec();
    std::vector<G4double> SfiberTimes;

    if(fiberPheVec.size()>0){
      G4cout << "S Fiber n " << current_Shit << "\tID: " << ShitCollection->GetSiPMID() << "\t number of packets: " << fiberPheVec.size() << "\tTotal number of phe: " << std::accumulate(fiberPheVec.begin(), fiberPheVec.end(), 0) << G4endl;
      for(int i=0; i<fiberPheVec.size(); i++){
        G4int NofPhe = fiberPheVec.at(i);
        for(int j=0; j<NofPhe; j++){
        G4double distance_to_sipm = moduleZ/2 - fiberZVec.at(i);
        G4double time = distance_to_sipm/vS;
        G4cout << "Position: " << fiberZVec.at(i) << "\tDistance to SiPM: " << distance_to_sipm << "\tTime: " << time <<  "\n";
        SfiberTimes.push_back(time);  // Array of photon timings to input to SimSiPM
        }
      }
      G4cout << G4endl;
      
    }

  }
  assert(fHitPheSvector.size() == fHitZcoordSvector.size());

  // Get individual hits in C fibers
  for(int current_Chit=0; current_Chit<CnOfHitsCollections; current_Chit++)
  {
    auto ChitCollection = (*CfiberHC)[current_Chit];
    fHitPheCvector.push_back(ChitCollection->GetPhe());
    fHitZcoordCvector.push_back(ChitCollection->GetZcoord());
    fHitSiPMIDCvector.push_back(ChitCollection->GetSiPMID());

    std::vector<G4double> fiberPheVec = ChitCollection->GetPheVec();
    std::vector<G4double> fiberZVec = ChitCollection->GetZVec();
    std::vector<G4double> CfiberTimes;

    if(fiberPheVec.size()>0){
      G4cout << "C Fiber n " << current_Chit << "\tID: " << ChitCollection->GetSiPMID() << "\t number of packets: " << fiberPheVec.size() << "\tTotal number of phe: " << std::accumulate(fiberPheVec.begin(), fiberPheVec.end(), 0) << G4endl;
        for(int i=0; i<fiberPheVec.size(); i++){
          G4int NofPhe = fiberPheVec.at(i);
          for(int j=0; j<NofPhe; j++){
          G4double distance_to_sipm = moduleZ/2 - fiberZVec.at(i);
          G4double time = distance_to_sipm/vC;
          G4cout << "Position: " << fiberZVec.at(i) << "\tDistance to SiPM: " << distance_to_sipm << "\tTime: " << time <<  "\n";
          CfiberTimes.push_back(time);  // Array of photon timings to input to SimSiPM
          }
        }
        G4cout << G4endl;
    
    }
  }
  assert(fHitPheCvector.size() == fHitZcoordCvector.size());
  



  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  //auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
  G4cout << "---> End of event: " << eventID << G4endl;     



  //Add all p.e. in Scin and Cher fibers before calibration
  //
  for (auto& n : VectorSignals) NofScinDet += n;
  for (auto& n : VecSPMT) NofScinDet += n;
  for (auto& n : VectorSignalsCher) NofCherDet += n;
  for (auto& n : VecCPMT) NofCherDet += n;

  //Fill ntuple event by event
  //entries with vectors are automatically filled
  //
  G4cout << "Filling histos" << G4endl;
  analysisManager->FillNtupleDColumn(1, 0, EnergyScin);
  analysisManager->FillNtupleDColumn(1, 1, EnergyCher);
  analysisManager->FillNtupleDColumn(1, 2, NofCherDet);
  analysisManager->FillNtupleDColumn(1, 3, NofScinDet);
  analysisManager->FillNtupleDColumn(1, 4, EnergyTot);
  analysisManager->FillNtupleDColumn(1, 5, PrimaryParticleEnergy);
  analysisManager->FillNtupleIColumn(1, 6, PrimaryPDGID);
  analysisManager->FillNtupleDColumn(1, 7, EscapedEnergyl);
  analysisManager->FillNtupleDColumn(1, 8, EscapedEnergyd);
  analysisManager->FillNtupleDColumn(1, 9, PSEnergy);
  analysisManager->FillNtupleDColumn(1, 10, PrimaryX);
  analysisManager->FillNtupleDColumn(1, 11, PrimaryY);


  analysisManager->AddNtupleRow(1);    // Remember this otherwise data is not printed on file

}

//**************************************************
