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
    //G4cout << "Fast Stepping Action called, inside volume " << volume->GetName() << G4endl;

    //--------------------------------------------------
    //Store information from Scintillation and Cherenkov
    //signals
    //--------------------------------------------------
   
    std::string Fiber;
    std::string S_fiber = "S_fiber";
    std::string C_fiber = "C_fiber";
    Fiber = volume->GetName(); 
    G4int TowerID;
    G4int SiPMID = 999999;
    G4int SiPMTower;
    G4int signalhit = 0;
    //G4double zdep = 0.;

    /**************************/
    /** SCINTILLATING FIBRES **/
    /**************************/

    if ( strstr( Fiber.c_str(), S_fiber.c_str() ) )         
    { 

        if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() )
        {
            // Check if optical photon has total internal reflection            
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
                    double phEne = step->GetTrack()->GetKineticEnergy();
                    double lambda = 1.24/(phEne*10e6)*10e3;
                    //G4cout << lambda << G4endl;
                    if(lambda > 470.)
                    {  // KODAK Wratten 3 Optical Filter
                        G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);
                        
                        TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));
                        SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
                        fEventAction->AddScin(edep);
                        //signalhit = fSignalHelper->SmearCSignal( fSignalHelper->ApplyBirks( edep, steplength ) );
                        //signalhit = 1;  // single optical photon
                        signalhit = fSignalHelper->SmearSSignal_new( );  // single optical photon

                        // Attenuate Signal (independent of PMT/SiPM module)
                        double post_signalhit = fSignalHelper->AttenuateSSignal_WL(signalhit, distance_to_sipm, lambda);
                        //double post_signalhit = fSignalHelper->AttenuateSSignal(signalhit, distance_to_sipm);


                        // For PMTs
                        double PMTpde = fSignalHelper->GetSPMTpde(lambda);
                        fEventAction->AddVecSPMT( TowerID, post_signalhit*PMTpde); 
                        //G4cout << lambda << "\t" << distance_to_sipm << "\t" << signalhit << "\t" << post_signalhit << "\t" << pde << G4endl;


                        if(SiPMTower > -1){ 
                                SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
                                double SiPMpde = fSignalHelper->GetSsipmpde(lambda);
                                //G4cout << "S fibre: " << signalhit << "\tafter attenuation: " << post_signalhit << "\tpde: " << SiPMpde << "\tResulting signal: " << SiPMpde*post_signalhit << G4endl;

                                fEventAction->AddVectorScin(SiPMTower*NoFibersTower+SiPMID, SiPMpde*post_signalhit); 
                        }
                    }

                    step->GetTrack()->SetTrackStatus( fStopAndKill ); // (do not propagate optical photons)

                }   // end of total internal reflection
                default:
                    ;
            }   // end of optical photon


            //step->GetTrack()->SetTrackStatus( fStopAndKill ); 
        }

        
        if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return; } //not ionizing particle
        fEventAction->AddScin(edep);

        /*
        G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);
        
        TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));
        SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
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

        */
    }   // ### END OF SCINTILLATING FIBRES ###



    /*********************/
    /** CERENKOV FIBRES **/
    /*********************/
    else if ( strstr( Fiber.c_str(), C_fiber.c_str() ) )     //Cherenkov fiber/tube
    { 
        fEventAction->AddCher(edep);

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
            switch ( theStatus ){
                                    
            case TotalInternalReflection:
            {
                double phEne = step->GetTrack()->GetKineticEnergy();
                double lambda = 1.24/(phEne*10e6)*10e3;
                if(lambda > 300.)
                {
                    //G4cout << lambda << G4endl;
                    G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);

                    //G4int c_signal = fSignalHelper->SmearCSignal( ); // return random variable with poissonian distribution around 0.153
                    //G4int c_signal = 1;                            // in case of no smearing study 
                    // Attenuate Signal
                    //c_signal = fSignalHelper->AttenuateCSignal(c_signal, distance_to_sipm);

                    TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));		
                    SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
                    //fEventAction->AddVecCPMT( TowerID, c_signal );
                    //		    if ( TowerID != 0 ) { fEventAction->AddVecCPMT( TowerID, c_signal ); }
                    G4int c_signal = fSignalHelper->SmearCSignal_new( );  // single optical photon

                    // Attenuate Signal (independent of PMT/SiPM module)
                    G4int post_signalhit = fSignalHelper->AttenuateCSignal_WL(c_signal, distance_to_sipm, lambda);
                    //G4cout << "Signal before att: " << c_signal << "\t after: " << post_signalhit << G4endl;

                    // For PMTs
                    G4double PMTpde = fSignalHelper->GetCPMTpde(lambda);
                    fEventAction->AddVecCPMT( TowerID, post_signalhit*PMTpde); 
                    //G4cout << lambda << "\t" << distance_to_sipm << "\t" << signalhit << "\t" << post_signalhit << "\t" << pde << G4endl;


                    if(SiPMTower > -1){ 
                            SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
                            G4double SiPMpde = fSignalHelper->GetSsipmpde(lambda);
                            //G4cout << "C fibre: " << c_signal << "\tafter attenuation: " << post_signalhit << "\tpde: " << SiPMpde << "\tResulting signal: " << SiPMpde*post_signalhit << G4endl;
                            //G4cout << "Wavelength: " << lambda << "\tPDE: " << SiPMpde << G4endl;
                            //fEventAction->AddVectorCher( SiPMTower*NoFibersTower+SiPMID, SiPMpde*post_signalhit); 
                            fEventAction->AddVectorCher(SiPMTower*NoFibersTower+SiPMID, SiPMpde*post_signalhit); 

                    }
                }
                step->GetTrack()->SetTrackStatus( fStopAndKill ); // (do not propagate optical photons)


                //if(SiPMTower > -1)
                //{ 
                //    SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
                //    //zdep = fDetConstruction->GetSiPMID(step->GetTrack()->GetPosition().z() );
                //    fEventAction->AddVectorCher(SiPMTower*NoFibersTower+SiPMID, c_signal);
                //    //fEventAction->AddVecCerZdep( SiPMTower*NoFibersTower+SiPMID, c_signal*zdep);
                //}
                //step->GetTrack()->SetTrackStatus( fStopAndKill );   // Kill photon to not propagate it inside fibres
            }
            default:
                ;
                //step->GetTrack()->SetTrackStatus( fStopAndKill );
            } //end of swich cases
            


        } //end of optical photon
        else return;

    } //end of Cherenkov fiber
    else return;
    //G4cout << "Fast Stepping Action ended well! " << G4endl;


}


//**************************************************
