//**************************************************
// \file DREMTubesEventAction.cc
// \brief: Implementation of DREMTubesEventAction 
//         class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "DREMTubesEventAction.hh"
#include "DREMTubesRunAction.hh"
#include "DREMTubesDetectorConstruction.hh"
#include "G4HCofThisEvent.hh"

//includers for calorimeter hits
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "DREMTubesCalorimeterHit.hh"
#include "DREMTubesCalorimeterSD.hh"
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




//Define constructor
//
DREMTubesEventAction::DREMTubesEventAction()
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
    fSfiberHCID(-1),
    fCfiberHCID(-1)
    {
}

//Define de-constructor
//
DREMTubesEventAction::~DREMTubesEventAction() {}



DREMTubesCalorimeterHitsCollection* 
DREMTubesEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<DREMTubesCalorimeterHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("DREMTubesEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    



//void DREMTubesEventAction::PrintEventStatistics(
//                              G4double SPhe, G4double SZcoord, G4int SSiPMID,
//                              G4double CPhe, G4double CZcoord, G4int CSiPMID) const
//{
  // print event statistics
  /*G4cout
     << "   Sfiber Total Photoelectrons: " 
     << std::setw(7) << G4BestUnit(SPhe, "eplus")
     << "       Sfiber Z coord: " 
     << std::setw(7) << G4BestUnit(SZcoord, "Length")
     << G4endl
     << "        S fiber ID: " 
     << std::setw(7) << G4BestUnit(SSiPMID, "eplus")
     << "       Cfiber total eergty: " 
     << std::setw(7) << G4BestUnit(CPhe, "eplus")
     << G4endl;*/
//     G4cout << std::setw(7) << "Total photoelectrons in S fibers: " << SPhe << G4endl;
//     G4cout << std::setw(7) << "Depth: " << SZcoord << G4endl;
//}






//Define BeginOfEventAction() and EndOfEventAction() methods
//
void DREMTubesEventAction::BeginOfEventAction(const G4Event*) {  
    
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
    fHitPheSvector.assign(fHitPheSvector.size(), 0.);
    fHitZcoordSvector.assign(fHitZcoordSvector.size(), 0.);
    fHitSiPMIDSvector.assign(fHitSiPMIDSvector.size(), 0.);
    fHitPheSvector.assign(fHitPheCvector.size(), 0.);
    fHitZcoordSvector.assign(fHitZcoordCvector.size(), 0.);
    fHitSiPMIDSvector.assign(fHitSiPMIDCvector.size(), 0.);

    //NofHitSInEvt = 0;
    //HitEventID = 0;
}









void DREMTubesEventAction::EndOfEventAction(const G4Event* event) {
 
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();



  
  // Get hits collections IDs (only once)
  if ( fSfiberHCID == -1 ) {
    fSfiberHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("SfiberHitsCollection");
    fCfiberHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("CfiberHitsCollection");
  }



  // Get hits collections
  auto SfiberHC = GetHitsCollection(fSfiberHCID, event);
  auto CfiberHC = GetHitsCollection(fCfiberHCID, event);
  int SnOfHitsCollections = SfiberHC->entries();
  int CnOfHitsCollections = CfiberHC->entries();


  //DREMTubesCalorimeterHitsCollection* myShitCollection = 0;
  //if(SfiberHC){
  //  myShitCollection = (DREMTubesCalorimeterhitsCollection*)(SfiberHC->GetHC(fSfiberHCID));
  //}

  //G4cout << "S fibers: number of hit collections: " << SnOfHitsCollections << G4endl;
  //G4cout << "C fibers: number of hit collections: " << CnOfHitsCollections << G4endl;

 // Get hit with total values
  //auto SfiberHit = (*SfiberHC)[SfiberHC->entries()-1];
  //auto CfiberHit = (*CfiberHC)[CfiberHC->entries()-1];
 
  //if ( fSfiberHCID >= 0)
  //{
  //  G4THitsMap<G4double>* hitMapS = dynamic_cast<G4THitsMap<G4double>*>(SfiberHC);
  //}

  std::vector<double> HitPheSvector, HitZcoordSvector, HitSiPMIDSvector;
  std::vector<double> HitPheCvector, HitZcoordCvector, HitSiPMIDCvector;

  //G4cout << "Number of entries  in SfiberHitCollection: " << SfiberHC->size() << G4endl;
  //G4double TotPheS = 0;
  // Get individual hits in S fibers
  for(int current_Shit=0; current_Shit<SnOfHitsCollections; current_Shit++)
  {
    auto ShitCollection = (*SfiberHC)[current_Shit];
    //G4double HitTower = ShitCollection->GetTowerID();
    G4double HitPhe = ShitCollection->GetPhe();
    G4int HitFiber = ShitCollection->GetSiPMID();
    G4double HitZ = ShitCollection->GetZcoord();
    fHitPheSvector.push_back(HitPhe);
    fHitZcoordSvector.push_back(HitZ);
    fHitSiPMIDSvector.push_back(HitFiber);
    //G4cout << "Hit number: " << current_Shit << "\tphe" << HitPhe << "\tZcoord: " << HitZ << "\tSiPM: " << HitFiber << G4endl;
  }
  assert(fHitPheSvector.size() == fHitZcoordSvector.size());

  // Get individual hits in C fibers
  for(int current_Chit=0; current_Chit<CnOfHitsCollections; current_Chit++)
  {
    auto ChitCollection = (*CfiberHC)[current_Chit];
    //G4double HitTower = ChitCollection->GetTowerID();
    G4double HitPhe = ChitCollection->GetPhe();
    G4int HitFiber = ChitCollection->GetSiPMID();
    G4double HitZ = ChitCollection->GetZcoord();
    fHitPheCvector.push_back(HitPhe);
    fHitZcoordCvector.push_back(HitZ);
    fHitSiPMIDCvector.push_back(HitFiber);
    //G4cout << "Hit number: " << current_Chit << "\tphe" << HitPhe << "\tZcoord: " << HitZ << "\tCiPM: " << HitFiber << G4endl;
  }
  assert(fHitPheCvector.size() == fHitZcoordCvector.size());

  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  //auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  //if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;     

    //PrintEventStatistics(
    //  SfiberHit->GetPhe(), SfiberHit->GetZcoord(), SfiberHit->GetSiPMID(),
    //  CfiberHit->GetPhe(), CfiberHit->GetZcoord(), CfiberHit->GetSiPMID());
  //}  
  

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
    //analysisManager->FillNtupleIColumn(1, 12, eventID);
    //analysisManager->FillNtupleDColumn(1, 13, TotPheS);
    //analysisManager->FillNtupleIColumn(1, 14, NofHitSInEvt);


    analysisManager->AddNtupleRow(1);    // Remember this otherwise data is not printed on file
 
    //analysisManager->FillNtupleIColumn(0, HitEventID);
    //analysisManager->AddNtupleRow(1);
 
    //analysisManager->FillNtupleDColumn(12, hitMapS);


}

//**************************************************
