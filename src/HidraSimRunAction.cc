//**************************************************
// \file HidraSimRunAction.cc 
// \brief: Implementation of HidraSimRunAction class 
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "HidraSimRunAction.hh"
#include "HidraSimEventAction.hh"

//Includers from Geant4
//
#include "g4root.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//Includers from C++
//
#include <string>

//Define constructor
//
HidraSimRunAction::HidraSimRunAction( HidraSimEventAction* eventAction )
    : G4UserRunAction(),
      fEventAction( eventAction ){ 
  
    //print event number per each event (default, can be overwritten with macro)
    //
    G4RunManager::GetRunManager()->SetPrintProgress(1);     

    //Instantiate analysis manager
    //
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel( 1 );
    analysisManager->SetNtupleMerging( 1 );

    // Set ntuple ID to one (if more than one ntuple is saved)
    analysisManager->SetFirstNtupleId(1);

    //Using ROOT as analysisManager type, print it
    //
    G4cout << "HidraSim-> Using " << analysisManager->GetType() << G4endl;

    //Define ntuple structure
    //
    analysisManager->CreateNtuple("HidraSimout", "HidraSimoutput");
    analysisManager->CreateNtupleDColumn("EnergyScin");                     //0
    analysisManager->CreateNtupleDColumn("EnergyCher");                     //1
    analysisManager->CreateNtupleDColumn("NofCherDet");                     //2
    analysisManager->CreateNtupleDColumn("NofScinDet");                     //3
    analysisManager->CreateNtupleDColumn("EnergyTot");                      //4
    analysisManager->CreateNtupleDColumn("PrimaryParticleEnergy");          //5
    analysisManager->CreateNtupleIColumn("PrimaryPDGID");                   //6
    analysisManager->CreateNtupleDColumn("EscapedEnergyl");                  //7
    analysisManager->CreateNtupleDColumn("EscapedEnergylup");                  //8
    analysisManager->CreateNtupleDColumn("EscapedEnergyldown");                  //9
    analysisManager->CreateNtupleDColumn("EscapedEnergylright");                  //10
    analysisManager->CreateNtupleDColumn("EscapedEnergylleft");                  //11    
    analysisManager->CreateNtupleDColumn("EscapedEnergyd");                  //12
    analysisManager->CreateNtupleDColumn("PSEnergy");                       //13
    analysisManager->CreateNtupleDColumn("PrimaryX");                       //14
    analysisManager->CreateNtupleDColumn("PrimaryY");                       //15

    analysisManager->CreateNtupleDColumn("VectorSignals", fEventAction->GetVectorSignals());
    analysisManager->CreateNtupleDColumn("VectorSignalsCher", fEventAction->GetVectorSignalsCher());
    analysisManager->CreateNtupleDColumn("VecTowerE", fEventAction->GetVecTowerE());
    analysisManager->CreateNtupleDColumn("VecSPMT", fEventAction->GetVecSPMT());
    analysisManager->CreateNtupleDColumn ("VecCPMT", fEventAction->GetVecCPMT());
    analysisManager->CreateNtupleDColumn("HitPheSvector", fEventAction->GetHitPheSvector());
    analysisManager->CreateNtupleDColumn("HitZcoordSvector", fEventAction->GetHitZcoordSvector());
    analysisManager->CreateNtupleIColumn("HitSiPMIDSvector", fEventAction->GetHitSiPMIDSvector());
    analysisManager->CreateNtupleDColumn("HitPheCvector", fEventAction->GetHitPheCvector());
    analysisManager->CreateNtupleDColumn("HitZcoordCvector", fEventAction->GetHitZcoordCvector());
    analysisManager->CreateNtupleIColumn("HitSiPMIDCvector", fEventAction->GetHitSiPMIDCvector());
    
    analysisManager->FinishNtuple(1);

    //analysisManager->CreateNtuple("SfiberHits", "SfiberHits");
    //analysisManager->CreateNtupleIColumn("HitEventID");                    //0
    //analysisManager->CreateNtupleDColumn("HitZcoord");                     //1
    //analysisManager->CreateNtupleIColumn("SiPMID");                        //2
    //analysisManager->CreateNtupleIColumn("TowerID");                       //3
    //analysisManager->CreateNtupleDColumn("pheS");                          //4
    //analysisManager->CreateNtupleDColumn("TrueX");                         //5
    //analysisManager->CreateNtupleDColumn("TrueY");                         //6
    //analysisManager->FinishNtuple(2);

    //analysisManager->CreateNtuple("CfiberHits", "CfiberHits");
    //analysisManager->CreateNtupleIColumn("HitEventID");                    //0
    //analysisManager->CreateNtupleDColumn("HitZcoord");                     //1
    //analysisManager->CreateNtupleIColumn("SiPMID");                        //2
    //analysisManager->CreateNtupleIColumn("TowerID");                       //3
    //analysisManager->CreateNtupleDColumn("pheC");                          //4
    //analysisManager->CreateNtupleDColumn("TrueX");                         //5
    //analysisManager->CreateNtupleDColumn("TrueY");                         //6
    //analysisManager->FinishNtuple(3);



}

//Define de-constructor
//
HidraSimRunAction::~HidraSimRunAction(){
   
    //Delete only instance of G4AnalysisManager
    //
    delete G4AnalysisManager::Instance();  

}

//Define BeginOfRunAction() and EndOfRunAction() methods
//
void HidraSimRunAction::BeginOfRunAction( const G4Run* Run )  { 
    
    //Save random seeds (optional)
    //
    //G4RunManager::GetRunManager()->SetRandomNumberStore( true );
    
    //Open output file, one per Run
    //
    auto analysisManager = G4AnalysisManager::Instance();
    std::string runnumber = std::to_string( Run->GetRunID() );
    G4String outputfile = "HidraSimout_Run"+runnumber;
    analysisManager->OpenFile( outputfile );
}

void HidraSimRunAction::EndOfRunAction( const G4Run* ) {
  
    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->Write();
    analysisManager->CloseFile();

}

//**************************************************
