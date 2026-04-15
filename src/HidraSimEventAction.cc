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
//#include "g4root.hh"
#include "G4AnalysisManager.hh"
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
//#include <SiPMProperties.h>
//#include <sipm/SiPMSensor.h>
//#include <sipm/SiPMAnalogSignal.h>
#include "SiPM.h"
#include "SiPMProperties.h"
using namespace sipm;

// Create a SiPMSensor object
//SiPMSensor mySensor(myProperties);



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
    EmEnergy(0.),
    NeutronEkin(0.),
    PionCount(0),
    NeutronCount(0),
    EventID(0),
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
void HidraSimEventAction::BeginOfEventAction(const G4Event* event) {  
    
    //Initialize data memebers at begin of each event
    //
    EventID = event->GetEventID();
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
    EmEnergy = 0.;
    NeutronEkin = 0.;
    PionCount = 0;
    NeutronCount = 0;
    NeutronEnergies.clear(); // reset per-event list

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
        //G4cout << "Number of hit collections: " << numberOfHitCollections << G4endl;

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

    
    //double intStart = 8.33;
    double intStart = 0.1;
    double intGate = 100.;
    //G4double integrationTime = 30;

  // Scintillating SiPMs
  // Create a SiPMProperties object
  SiPMProperties mySciProperties;
  mySciProperties.setDcr(0);
  mySciProperties.setProperty("Pitch", 10);
  mySciProperties.setProperty("Size", 1.00);
  mySciProperties.setSampling(0.10); // 100 ps
  mySciProperties.setRecoveryTime(20.);
  mySciProperties.setRiseTime(1.);
  mySciProperties.setFallTimeFast(30.);

  // Create a SiPMSensor object
  SiPMSensor mySciSensor(mySciProperties);
  // Change parameters
  //mySensor.properties().setAp(0.01);    // Using proper getter/setter
  //mySensor.setProperty("Pitch", 25);    // Using parameter name

  // Cerenkov SiPMs
  // Create a SiPMProperties object
  SiPMProperties myCerProperties;
  myCerProperties.setDcr(0);
  myCerProperties.setProperty("Pitch", 15);
  myCerProperties.setProperty("Size", 1.00);
  myCerProperties.setSampling(0.10);  // 100 ps
  myCerProperties.setRecoveryTime(20.);
  myCerProperties.setRiseTime(1.);
  myCerProperties.setFallTimeFast(30.);
  //myCerProperties.setRecoveryTime(10);
  // Create a SiPMSensor object
  SiPMSensor myCerSensor(myCerProperties);
  // Change parameters
  //mySensor.properties().setAp(0.01);    // Using proper getter/setter
  //mySensor.setProperty("Pitch", 25);    // Using parameter name





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
    
    std::vector<G4double> fiberPheVec = ShitCollection->GetPheVec();
    std::vector<G4double> fiberZVec = ShitCollection->GetZVec();
    std::vector<G4double> SfiberTimes;

    if(fiberPheVec.size()>0)
    {
        // Each packet shares the same Z position along fibre ~ same time
        //G4cout << "S Fiber n " << current_Shit << "\tID: " << ShitCollection->GetSiPMID() << "\t number of packets: " << fiberPheVec.size() << "\tTotal number of phe: " << std::accumulate(fiberPheVec.begin(), fiberPheVec.end(), 0) << G4endl;
        double nPheS = std::accumulate(fiberPheVec.begin(), fiberPheVec.end(), 0);
        for(int i=0; i<fiberPheVec.size(); i++){
            G4int NofPhe = fiberPheVec.at(i);
            for(int j=0; j<NofPhe; j++){
                G4double phTime = fiberZVec.at(i)/vS;
                SfiberTimes.push_back(phTime/1000.);  // Array of photon timings to input to SimSiPM
            }
        }

        mySciSensor.resetState();
        mySciSensor.addPhotons(SfiberTimes);    // Sets photon times (times are in ns) (not appending)
        mySciSensor.runEvent();           // Runs the simulation
        SiPMAnalogSignal mySciSignal = mySciSensor.signal();
        double integral = mySciSignal.integral(intStart,intGate,0.5);   // (intStart, intGate, threshold)
        double peak = mySciSignal.peak(intStart,intGate,0.5);   // (intStart, intGate, threshold)
        double toa = mySciSignal.toa(intStart,intGate,0.5);   // (intStart, intGate, threshold)
        double tot = mySciSignal.tot(intStart,intGate,0.5);   // (intStart, intGate, threshold)
        //fHitPheSvector.push_back(integral);     // fill with SiPM signal intergral
        fHitPheSvector.push_back(nPheS);     // fill with total number of phe
        fHitZcoordSvector.push_back(toa+intStart);   // fill with arrival time of first phe
        //fHitZcoordSvector.push_back(toa);   // fill with arrival time of first phe
        fHitSiPMIDSvector.push_back(current_Shit);        // fiber ID
    }  


  }

  // Get individual hits in C fibers
  for(int current_Chit=0; current_Chit<CnOfHitsCollections; current_Chit++)
  {
    auto ChitCollection = (*CfiberHC)[current_Chit];
    std::vector<G4double> fiberPheVec = ChitCollection->GetPheVec();
    std::vector<G4double> fiberZVec = ChitCollection->GetZVec();
    std::vector<G4double> CfiberTimes;

    if(fiberPheVec.size()>0)
    {
      double nPheC = std::accumulate(fiberPheVec.begin(), fiberPheVec.end(), 0);
        for(int i=0; i<fiberPheVec.size(); i++)
        {
            G4int NofPhe = fiberPheVec.at(i);
            for(int j=0; j<NofPhe; j++){
              G4double phTime = fiberZVec.at(i)/vC;
              CfiberTimes.push_back(phTime/1000.);  // Array of photon timings to input to SimSiPM
            }
        }

      myCerSensor.resetState();
      myCerSensor.addPhotons(CfiberTimes);    // Sets photon times (times are in ns) (not appending)
      myCerSensor.runEvent();           // Runs the simulation
      SiPMAnalogSignal myCerSignal = myCerSensor.signal();

      double integral = myCerSignal.integral(intStart,intGate,0.5);   // (intStart, intGate, threshold)
      double peak = myCerSignal.peak(intStart,intGate,0.5);   // (intStart, intGate, threshold)
      double toa = myCerSignal.toa(intStart,intGate,0.5);   // (intStart, intGate, threshold)
      double tot = myCerSignal.tot(intStart,intGate,0.5);   // (intStart, intGate, threshold)
      //fHitPheCvector.push_back(integral);     // fill with SiPM signal integral
      fHitPheCvector.push_back(nPheC);     // fill with total number of phe
      fHitZcoordCvector.push_back(toa+intStart);   // fill with arrival time of first phe
      fHitSiPMIDCvector.push_back(current_Chit);        // fiber ID
    }  


  }
  



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
  analysisManager->FillNtupleIColumn(1, 0, EventID);
  analysisManager->FillNtupleDColumn(1, 1, EnergyScin);
  analysisManager->FillNtupleDColumn(1, 2, EnergyCher);
  analysisManager->FillNtupleDColumn(1, 3, NofCherDet);
  analysisManager->FillNtupleDColumn(1, 4, NofScinDet);
  analysisManager->FillNtupleDColumn(1, 5, EnergyTot);
  analysisManager->FillNtupleDColumn(1, 6, PrimaryParticleEnergy);
  analysisManager->FillNtupleIColumn(1, 7, PrimaryPDGID);
  analysisManager->FillNtupleDColumn(1, 8, EscapedEnergyl);
  analysisManager->FillNtupleDColumn(1, 9, EscapedEnergyd);
  analysisManager->FillNtupleDColumn(1, 10, PSEnergy);
  analysisManager->FillNtupleDColumn(1, 11, PrimaryX);
  analysisManager->FillNtupleDColumn(1, 12, PrimaryY);
  analysisManager->FillNtupleDColumn(1, 13, EmEnergy);
  analysisManager->FillNtupleDColumn(1, 14, NeutronEkin);
  analysisManager->FillNtupleIColumn(1, 15, PionCount);
  analysisManager->FillNtupleIColumn(1, 16, NeutronCount);


  analysisManager->AddNtupleRow(1);    // Remember this otherwise data is not printed on file

}

void HidraSimEventAction::AddPionCount() {
    PionCount++;
}

void HidraSimEventAction::AddNeutronCount() {
    NeutronCount++;
}

//**************************************************
