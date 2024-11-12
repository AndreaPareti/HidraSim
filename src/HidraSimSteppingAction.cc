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

    if (    volume->GetName() == "Clad_S_fiber" ||
            volume->GetName() == "Core_S_fiber" ||
            volume->GetName() == "Abs_Scin_fiber"  ||
            volume->GetName() == "Clad_C_fiber" ||
            volume->GetName() == "Core_C_fiber" ||
            volume->GetName() == "Abs_Cher_fiber"  ) {
            fEventAction->AddVecTowerE(fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3)),
				  edep );
    }
    	
    if ( volume->GetName() == "Preshower_scin" || volume->GetName() == "Preshower_pb"){
        fEventAction->AddPSEnergy( edep );
    }

    if(volume->GetName() == "leakbox"){
        G4int LCID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
        fEventAction->AddVecLeakCounter(LCID, edep); 

    }
    
    if ( volume != fDetConstruction->GetWorldPV() &&
         volume != fDetConstruction->GetLeakCntlPV() &&
         volume != fDetConstruction->GetLeakCntdPV() &&
         volume->GetName() != "Preshower_scin" &&
         volume->GetName() != "Preshower_pb" &&
         volume->GetName() != "leakbox") { fEventAction->Addenergy(edep);}
   
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
    G4double signalhit = 0;
    //G4double zdep = 0.;

    /**************************/
    /** SCINTILLATING FIBRES **/
    /**************************/

    if ( strstr( Fiber.c_str(), S_fiber.c_str() ) )
    {
        if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() )
        {
            step->GetTrack()->SetTrackStatus( fStopAndKill ); 
	    }

        if ( step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0. ) { return; } //not ionizing particle		 
        TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));
        SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
        fEventAction->AddScin(edep);
        signalhit = fSignalHelper->SmearSSignal( fSignalHelper->ApplyBirks( edep, steplength ) );

        G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);

        signalhit = fSignalHelper->AttenuateSSignal(signalhit, distance_to_sipm);
        fEventAction->AddVecSPMT( TowerID, signalhit ); 

        if(SiPMTower > -1)
        { 
            SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
            fEventAction->AddVectorScin( SiPMTower*NoFibersTower + SiPMID , signalhit ); 
            //fEventAction->AddVectorScin( SiPMID+NofFibersrow*NofFiberscolumn*SiPMTower/2, signalhit ); 
            //fEventAction->AddVectorScin( SiPMID , signalhit ); 



        }
    }
    // End Scintillating Fibers case


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
            switch ( theStatus )
            {
                            
                case TotalInternalReflection:
                {
                    G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);
                    G4double c_signal = fSignalHelper->SmearCSignal( );
                    // Attenuate Signal
                    c_signal = fSignalHelper->AttenuateCSignal(c_signal, distance_to_sipm);

                    G4int TowerID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);		
                    SiPMTower=fDetConstruction->GetSiPMTower(TowerID);

                    fEventAction->AddVecCPMT( TowerID, c_signal );

                    if(SiPMTower > -1)
                    { // in sipm-readout tower
                        G4int SiPMID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
                        fEventAction->AddVectorCher(SiPMTower*NoFibersTower+SiPMID, c_signal);
                        //fEventAction->AddVectorCher(SiPMID+NofFibersrow*NofFiberscolumn*SiPMTower/2, c_signal);
                        //fEventAction->AddVectorCher(SiPMID , c_signal);
                        
                    }
                    step->GetTrack()->SetTrackStatus( fStopAndKill );
                }
                default: 
                    step->GetTrack()->SetTrackStatus( fStopAndKill );
	        } //end of swich cases
        } //end of optical photon

        else return;

    } //end of Cherenkov fiber
    else return;


}


//******************************************/
/*       SteppingAction.cc ends here       */
//******************************************/










// Store here temporarily stepping action including photon wavelenght behaviour 
/*
    if ( strstr( Fiber.c_str(), S_fiber.c_str() ) )         
    { 

        if ( step->GetTrack()->GetParticleDefinition() == G4OpticalPhoton::Definition() )
        {
            // Maybe set optical filter before asking for total internal reflection?
            
            // Check if optical photon has total internal reflection            
            G4OpBoundaryProcessStatus theStatus = Undefined;

            G4ProcessManager* OpManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

            if (OpManager) 
            {
                G4int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
                G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

                // maybe can ask for ~few internal reflections instead of whole fiber
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
                    if(  (lambda > 300.) || (lambda < 600.) )
                    {  // KODAK Wratten 3 Optical Filter
                        G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);
                        
                        TowerID = fDetConstruction->GetTowerID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3));
                        SiPMTower=fDetConstruction->GetSiPMTower(TowerID);
                        fEventAction->AddScin(edep);
                        //signalhit = fSignalHelper->SmearCSignal( fSignalHelper->ApplyBirks( edep, steplength ) );

                        // Generate signal with smearing
                        signalhit = fSignalHelper->SmearSSignalOpticalPhoton( );  // single optical photon

                        // Attenuate signal depending on current longitudinal position  (independent of PMT/SiPM module)
                        double attenuated_signalhit = fSignalHelper->AttenuateSSignalOverWL(signalhit, distance_to_sipm, lambda);
                        //double post_signalhit = fSignalHelper->AttenuateSSignal(signalhit, distance_to_sipm);

                        // If no signal is produced, skip optical filter & optical det efficiencies and exit
                        if(attenuated_signalhit == 0){step->GetTrack()->SetTrackStatus( fStopAndKill ); return;}


                        // For PMTs
                        //double PMTpde = fSignalHelper->GetSPMTpde(lambda);
                        // Correction factor to correct for optical filter and PMT efficiency
                        double PMTcorrection = fSignalHelper->GetSpmtCorrection(lambda);
                        fEventAction->AddVecSPMT( TowerID, attenuated_signalhit*PMTcorrection); 
                        //G4cout << lambda << "\t" << distance_to_sipm << "\t" << signalhit << "\t" << attenuated_signalhit << "\t" << attenuated_signalhit*PMTcorrection << G4endl;


                        if(SiPMTower > -1){ 
                                SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
                                double SiPMcorrection = fSignalHelper->GetSsipmCorrection(lambda);
                                //G4cout << "S fibre: " << signalhit << "\tafter attenuation: " << post_signalhit << "\tpde: " << SiPMpde << "\tResulting signal: " << SiPMpde*post_signalhit << G4endl;
                                fEventAction->AddVectorScin(SiPMTower*NoFibersTower+SiPMID, attenuated_signalhit*SiPMcorrection); 
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



    }   // ### END OF SCINTILLATING FIBRES ###

*/




// Store here temporarily C photon stepping action including photon wavelength behaviour
/*

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
                    // Generate signal with smearing 
                    G4int c_signal = fSignalHelper->SmearCSignalOpticalPhoton( );  // single optical photon

                    // Attenuate signal depending on current longitudinal position (independent of PMT/SiPM module)
                    G4int attenuated_signalhit = fSignalHelper->AttenuateCSignalOverWL(c_signal, distance_to_sipm, lambda);
                    //G4cout << "Signal before att: " << c_signal << "\t after: " << post_signalhit << G4endl;

                    // If no signal is produced, skip optical filter & optical det efficiencies and exit
                    if(attenuated_signalhit == 0){step->GetTrack()->SetTrackStatus( fStopAndKill ); return;}

                    // For PMTs
                    //double PMTpde = fSignalHelper->GetSPMTpde(lambda);
                    // Correction factor to correct for optical filter and PMT efficiency
                    double PMTcorrection = fSignalHelper->GetCpmtCorrection(lambda);
                    fEventAction->AddVecCPMT( TowerID, attenuated_signalhit*PMTcorrection); 
                    //G4cout << lambda << "\t" << distance_to_sipm << "\t" << signalhit << "\t" << post_signalhit << "\t" << pde << G4endl;


                    if(SiPMTower > -1){ 
                            SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
                            G4double SiPMcorrection = fSignalHelper->GetSsipmCorrection(lambda);
                            //G4cout << "C fibre: " << c_signal << "\tafter attenuation: " << post_signalhit << "\tpde: " << SiPMpde << "\tResulting signal: " << SiPMpde*post_signalhit << G4endl;
                            //G4cout << "Wavelength: " << lambda << "\tPDE: " << SiPMpde << G4endl;
                            //fEventAction->AddVectorCher( SiPMTower*NoFibersTower+SiPMID, SiPMpde*post_signalhit); 
                            fEventAction->AddVectorCher(SiPMTower*NoFibersTower+SiPMID, SiPMcorrection*attenuated_signalhit); 

                    }
                }
                step->GetTrack()->SetTrackStatus( fStopAndKill ); // (do not propagate optical photons)

            }
            default:
                ;
                //step->GetTrack()->SetTrackStatus( fStopAndKill );
            } //end of swich cases
            


        } //end of optical photon
        */