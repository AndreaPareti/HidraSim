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
  // Create hits collection
  fHitsCollection = new HidraSimCalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce
  auto hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

  // Create hits
  // fNofActiveFibers for ActiveFibers + one more for total sums 
  //for (G4int i=0; i<fNofActiveFibers+1; i++ ) {
  //fHitsCollection->insert(new HidraSimCalorimeterHit());
  //}
  //for (G4int i=0; i<fNofActiveFibers; i++ ) {
  //fHitsCollection->insert(new HidraSimCalorimeterHit());
  //}
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HidraSimCalorimeterSD::ProcessHits(G4Step* step, 
                                     G4TouchableHistory*)
{
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


  /*****************************************************************/
  /************************* Cerenkov ******************************/
  /*****************************************************************/

  if( (volume->GetName() == "Core_C_fiber") || (volume->GetName() == "Clad_C_fiber") )   // C fibers
  {
    // Select photons that have total internal reflection in C fibers
    if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) // check if particle is optical photon
    {
          //int a = 0; 
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
          } // return "theStatus" variable of optical photon

 
          // Total Internal Reflection Requirement case
          switch ( theStatus )
          {   
              case TotalInternalReflection:
              {          
                G4int c_signal = fSignalHelper->SmearCSignal( ); // return random variable with poissonian distribution around 0.153
                G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);
                // Attenuate Signal
                c_signal = fSignalHelper->AttenuateCSignal(c_signal, distance_to_sipm);            
                if(c_signal==0){return false;}                        // if no photoelectron is produced, do not save Hit
                TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));		
                SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
                
                if(SiPMTower > -1)
                { 
                    G4StepPoint *preStepPoint = step->GetPreStepPoint();
                    G4ThreeVector position = preStepPoint->GetPosition();                
                    G4int SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
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
              default:    // If not total internal reflection, exit 
                return false;

          } //end of switch cases
        

    } //end of check optical photon
    return false;   // if not optical photon, exit
  } // end of C fibers
  



  /*****************************************************************/
  /*********************** Scintillation ***************************/
  /*****************************************************************/
  //if ( strstr( Fiber.c_str(), S_fiber.c_str() ) )
  else if ( (volume->GetName() == "Core_S_fiber") || (volume->GetName() == "Clad_S_fiber") )
  { 
  	//if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return false;  } //not ionizing particle
    if( (step->GetTrack()->GetDefinition()->GetPDGCharge() > 0) && (step->GetStepLength() > 0.) )
    {      
          TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));
          SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
          G4int signalhit = fSignalHelper->SmearSSignal( fSignalHelper->ApplyBirks( edep, steplength ) );
          // Attenuate Signal
          G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);
          signalhit = fSignalHelper->AttenuateSSignal(signalhit, distance_to_sipm);

          if (signalhit == 0.) {return false;}                       // if no photoelectron is produced, do not save hit
          

          
          if (SiPMTower > -1)                                        // save hits only in SiPM-mounted modules
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
          

            fHitsCollection->insert(hit);
            return true;
          }
          else{return false;}
          
    }
    else
    {
      return false;
    }

    



   // if scintillating fiber, but none of previous cases has been run, do not fill anything and return
    return false; 
  } //end Scintillating fiber case

  

  // if volume is not S or C fiber, do not save anything
  else
  {
    return false;
  } 
  
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
