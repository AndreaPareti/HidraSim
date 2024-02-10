//**************************************************
// \file HidraSimSteppingAction.cc
// \brief: Implementation of 
//         HidraSimSteppingAction.cc
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Includers from project files
//
#include "HidraSimSteppingAction.hh"
#include "HidraSimEventAction.hh"
#include "HidraSimDetectorConstruction.hh"

//Includers from Geant4
//
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"

//Define constructor
//
HidraSimSteppingAction::HidraSimSteppingAction( HidraSimEventAction* eventAction,
						  const HidraSimDetectorConstruction* detConstruction)
    : G4UserSteppingAction(),
    fEventAction(eventAction),
    fDetConstruction(detConstruction){
		
        fSignalHelper = HidraSimSignalHelper::Instance(); 
		
}

//Define de-constructor
//
HidraSimSteppingAction::~HidraSimSteppingAction() {}

//Define UserSteppingAction() method
//
void HidraSimSteppingAction::UserSteppingAction( const G4Step* step ) {
    
    //Save auxiliary information
    //
    AuxSteppingAction( step );

    //Save fast signal information
    //
    FastSteppingAction( step );
    //}

}

//Define AuxSteppingAction() method
//
void HidraSimSteppingAction::AuxSteppingAction( const G4Step* step ) {

    // Get step info
    //
    G4VPhysicalVolume* volume 
        = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4double edep = step->GetTotalEnergyDeposit();

    //--------------------------------------------------
    //Store auxiliary information from event steps
    //--------------------------------------------------

    //    if ( volume == fDetConstruction->GetLeakCntPV() ){
        //Take care operator== works with pointers only
	//if there is a single placement of the volume
	//use names or cpNo if not the case
	//
    if ( volume->GetName() == "leakageabsorberl"){
        fEventAction->AddEscapedEnergyl(step->GetTrack()->GetKineticEnergy());
        step->GetTrack()->SetTrackStatus(fStopAndKill);
    } 
    if ( volume->GetName() == "leakageabsorberd" ){
        fEventAction->AddEscapedEnergyd(step->GetTrack()->GetKineticEnergy());
        step->GetTrack()->SetTrackStatus(fStopAndKill);
    } 

    if ( volume->GetName() == "Clad_S_fiber" ||
         volume->GetName() == "Core_S_fiber" ||
	 volume->GetName() == "Abs_Scin_fiber"  ||
	 volume->GetName() == "Clad_C_fiber" ||
	 volume->GetName() == "Core_C_fiber" ||
         volume->GetName() == "Abs_Cher_fiber"  ) {
        fEventAction->AddVecTowerE(fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3)),
				  edep );
    }
    	
    if ( volume->GetName() == "Preshower_scin" || volume->GetName() == "Preshower_pb" ){
        fEventAction->AddPSEnergy( edep );
    }
    
    if ( volume != fDetConstruction->GetWorldPV() &&
         volume != fDetConstruction->GetLeakCntlPV() &&
         volume != fDetConstruction->GetLeakCntdPV() &&
         volume->GetName() != "Preshower_scin" &&
         volume->GetName() != "Preshower_pb" ) { fEventAction->Addenergy(edep); }
   
    if ( step->GetTrack()->GetTrackID() == 1 &&
        step->GetTrack()->GetCurrentStepNumber() == 1){
        //Save primary particle energy and name
        //
        fEventAction->SavePrimaryPDGID(step->GetTrack()->GetDefinition()->GetPDGEncoding());
        fEventAction->SavePrimaryEnergy(step->GetTrack()->GetVertexKineticEnergy());
        fEventAction->SavePrimaryXY(step->GetTrack()->GetPosition().x(),
                                    step->GetTrack()->GetPosition().y());

    }
}

//Define FastSteppingAction() method
//
void HidraSimSteppingAction::FastSteppingAction( const G4Step* step ) { 
		
    // Get step info
    //
    G4VPhysicalVolume* volume 
        = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4double edep = step->GetTotalEnergyDeposit();
    G4double steplength = step->GetStepLength();
    
    //--------------------------------------------------
    //Store information from Scintillation and Cherenkov
    //signals
    //--------------------------------------------------
   
    std::string Fiber;
    std::string S_fiber = "S_fiber";
    std::string C_fiber = "C_fiber";
    Fiber = volume->GetName(); 
    G4int TowerID;
    G4int SiPMID = 900;
    G4int SiPMTower;
    G4int signalhit = 0;
    //G4double zdep = 0.;

    
    if ( strstr( Fiber.c_str(), S_fiber.c_str() ) )         //scintillating fiber/tube
    { 

        if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ) {
            step->GetTrack()->SetTrackStatus( fStopAndKill ); 
        }

        if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return; } //not ionizing particle
        //    G4VPhysicalVolume* modvolume 
        //        = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume(3);
        //	 std::cout << " grandmother name " << modvolume->GetName() << " number " << modvolume->GetCopyNo() << std::endl;
        //        std::cout << " grandmother nunber " << step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3) << std::endl;			 
            
        G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);
        
        TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));
        SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
        fEventAction->AddScin(edep);
        signalhit = fSignalHelper->SmearSSignal( fSignalHelper->ApplyBirks( edep, steplength ) );
        // Attenuate Signal
        signalhit = fSignalHelper->AttenuateSSignal(signalhit, distance_to_sipm);
        fEventAction->AddVecSPMT( TowerID, signalhit ); 


        if(SiPMTower > -1){ 
                SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
                //zdep = fDetConstruction->GetSiPMID(step->GetTrack()->GetPosition().z() );
                fEventAction->AddVectorScin( signalhit, SiPMTower*NoFibersTower+SiPMID ); 
                //fEventAction->AddVecSciZdep( SiPMTower*NoFibersTower+SiPMID, signalhit*zdep);

        }
    }   // end of scintillating fiber




    else if ( strstr( Fiber.c_str(), C_fiber.c_str() ) )     //Cherenkov fiber/tube
    { 
        fEventAction->AddCher(edep);

        if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() ){
                        
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
            switch ( theStatus ){
                                    
            case TotalInternalReflection:
            {
                G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);

                G4int c_signal = fSignalHelper->SmearCSignal( ); // return random variable with poissonian distribution around 0.153
                //G4int c_signal = 1;                            // in case of no smearing study 
                // Attenuate Signal
                c_signal = fSignalHelper->AttenuateCSignal(c_signal, distance_to_sipm);

                TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));		
                SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
                fEventAction->AddVecCPMT( TowerID, c_signal );
                //		    if ( TowerID != 0 ) { fEventAction->AddVecCPMT( TowerID, c_signal ); }

                if(SiPMTower > -1)
                { 
                    SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
                    //zdep = fDetConstruction->GetSiPMID(step->GetTrack()->GetPosition().z() );
                    fEventAction->AddVectorCher(SiPMTower*NoFibersTower+SiPMID, c_signal);
                    //fEventAction->AddVecCerZdep( SiPMTower*NoFibersTower+SiPMID, c_signal*zdep);
                }
                step->GetTrack()->SetTrackStatus( fStopAndKill );
            }
            default:
                ;
                //step->GetTrack()->SetTrackStatus( fStopAndKill );
            } //end of swich cases
            


        } //end of optical photon
        else return;

    } //end of Cherenkov fiber
    else return;
   
}

//**************************************************
