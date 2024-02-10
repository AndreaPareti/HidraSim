//
/// \file HidraSimCalorimeterSD.cc
/// \brief Implementation of the HidraSimCalorimeterSD class

#include "HidraSimCalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "HidraSimDetectorConstruction.hh"
#include "HidraSimSignalHelper.hh"
#include <set>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HidraSimCalorimeterSD::HidraSimCalorimeterSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofActiveFibers)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fNofActiveFibers(nofActiveFibers)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HidraSimCalorimeterSD::~HidraSimCalorimeterSD() 
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HidraSimCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  //G4cout <<  "Initializing Hits collection " << collectionName[0] << " with Sensitive detector: " << SensitiveDetectorName << G4endl;
  // Create hits collection
  fHitsCollection = new HidraSimCalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  auto hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
  // Create hits
  // fNofActiveFibers for ActiveFibers + one more for total sums 
  //for (G4int i=0; i<fNofActiveFibers+1; i++ ) {
  //  fHitsCollection->insert(new HidraSimCalorimeterHit());
  //}
  //for (G4int i=0; i<fNofActiveFibers; i++ ) {
  //  fHitsCollection->insert(new HidraSimCalorimeterHit());
  //}
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HidraSimCalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{
	//G4cout << "ProcessHits process has been called" << G4endl;
  // Get step info
  G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4double edep = step->GetTotalEnergyDeposit();
  G4double steplength = step->GetStepLength();
  std::string Fiber;
  std::string S_fiber = "S_fiber";
  std::string C_fiber = "C_fiber";
  Fiber = volume->GetName();     
  G4int TowerID;
  G4int SiPMTower;
  G4double signalhit;


  // Scintillating fiber
  if ( strstr( Fiber.c_str(), S_fiber.c_str() ) )
  { 
  	if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return false;  } //not ionizing particle
    TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));
    SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
  	signalhit = fSignalHelper->SmearSSignal( fSignalHelper->ApplyBirks( edep, steplength ) );
    if (signalhit == 0.) {return false;}                       // if no photoelectron is produced, do not save hit
    if (SiPMTower >-1)                                        // save hits only in SiPM-mounted modules
    {
      G4StepPoint *preStepPoint = step->GetPreStepPoint();
      G4ThreeVector position = preStepPoint->GetPosition();

      G4int SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));

      // Get calorimeter cell id 
      //auto hit = (*fHitsCollection)[SiPMID+NofFibersrow*NofFiberscolumn*SiPMTower/2];
      auto hit  = new HidraSimCalorimeterHit();

      if ( ! hit ) {
        G4ExceptionDescription msg;
        msg << "Cannot access hit " ; 
        G4Exception("HidraSimCalorimeterSD::ProcessHits()",
          "MyCode0004", FatalException, msg);
      }         


      hit->SetZcoord(position.z());
      hit->SetTowerID(TowerID);
      hit->SetSiPMID(SiPMID+NofFibersrow*NofFiberscolumn*SiPMTower/2);      // Count fibers only in SiPM-mounted towers
      //hit->SetSiPMID(SiPMID+NofFibersrow*NofFiberscolumn*TowerID/2);          // Count fibers in all towers
      hit->SetPhe(signalhit);
    

      // Add values
      //hit->Add(signalhit);

      // Get hit for total accounting
      //auto hitTotal 
      //  = (*fHitsCollection)[fHitsCollection->entries()-1];
      //hitTotal->Add(edep); 
      //hitTotal->Add(signalhit); 
      //fEventAction->AddNewHit();
      //fEventAction->InsertHitPheS(signalhit);
      //G4cout << "New Hit added" << G4endl;

      //G4int HitEventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
      //auto man = G4AnalysisManager::Instance();
      //man->FillNtupleIColumn(2, 0, HitEventID);
      //man->FillNtupleDColumn(2, 1, position.z());
      //man->FillNtupleIColumn(2, 2, SiPMID);
      //man->FillNtupleIColumn(2, 3, TowerID);
      //man->FillNtupleDColumn(2, 4, signalhit);
      //man->FillNtupleDColumn(2, 5, position.x());
      //man->FillNtupleDColumn(2, 6, position.y());
      //man->AddNtupleRow(2);   // add row to the hit ntuple
      fHitsCollection->insert(hit);

      return true;
    }
  } //end Scintillating fiber case


  //Cherenkov fiber/tube
  else if ( strstr( Fiber.c_str(), C_fiber.c_str() ) )     
  { 
    // Select photons that have total internal reflection in C fibers
    if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() )
    {
      G4OpBoundaryProcessStatus theStatus = Undefined;
      G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
      if (OpManager) 
      {
        G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
        G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

        for ( G4int i=0; i<MAXofPostStepLoops; i++) 
        {
            G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
            fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
            if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break; }
        }
      }
      // Total Internal Reflection Requirement case
      switch ( theStatus )
      {   
          case TotalInternalReflection:
          {       
            G4double c_signal = fSignalHelper->SmearCSignal( ); // return random variable with poissonian distribution around 0.153
            //G4int c_signal = 1;                                // No Smearing 
            if(c_signal==0){return false;}                        // if no photoelectron is produced, do not save Hit
            TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));		
            SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
            if(SiPMTower > -1)
            { 
                G4StepPoint *preStepPoint = step->GetPreStepPoint();
                G4ThreeVector position = preStepPoint->GetPosition();                
                //G4int HitEventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
                //G4double z = fDetConstruction->GetSiPMID(step->GetTrack()->GetPosition().z() );
                G4int SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
                
                //auto man = G4AnalysisManager::Instance();
                //man->FillNtupleIColumn(3, 0, HitEventID);
                //man->FillNtupleDColumn(3, 1, position.z());
                //man->FillNtupleIColumn(3, 2, SiPMID);
                //man->FillNtupleIColumn(3, 3, TowerID);
                //man->FillNtupleDColumn(3, 4, c_signal);
                //man->AddNtupleRow(3);   // add row to the hit ntuple

                auto hit  = new HidraSimCalorimeterHit();

                if ( ! hit ) {
                  G4ExceptionDescription msg;
                  msg << "Cannot access hit " ; 
                  G4Exception("HidraSimCalorimeterSD::ProcessHits()",
                    "MyCode0004", FatalException, msg);
                }         

                hit->SetZcoord(position.z());
                hit->SetTowerID(TowerID);
                hit->SetSiPMID(SiPMID+NofFibersrow*NofFiberscolumn*SiPMTower/2);      // Count fibers only in SiPM-mounted towers
                //hit->SetSiPMID(SiPMID+NofFibersrow*NofFiberscolumn*TowerID/2);          // Count fibers in all towers
                hit->SetPhe(c_signal);
                // Add values
                fHitsCollection->insert(hit);
                return true;
            }
          }
          default:
          {
            return false;
          }
      
      } //end of swich cases
    } //end of optical photon
  } //end of Cherenkov fiber

  else
    return false;         // if volume is not S or C fiber, do not save anything
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HidraSimCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{



  if ( verboseLevel>1 ) { 
     G4int nofHitsCollections = fHitsCollection->entries();
     G4cout
       << G4endl 
       << "-------->Hits Collection: in this event there are " << nofHitsCollections 
       << " hits collections: " << G4endl;
     for ( G4int i=0; i<nofHitsCollections; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
