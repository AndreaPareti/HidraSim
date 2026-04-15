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

#include "G4PhysicsModelCatalog.hh"   // optional if using GetCreatorModelName()
#include "G4VProcess.hh"

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
            G4Track* track = step->GetTrack();
            G4String particleName = track->GetDefinition()->GetParticleName();
            if (particleName == "e-" || particleName == "e+" || particleName == "gamma" ) {
                G4double EmEdep = step->GetTotalEnergyDeposit();
                fEventAction->AddEmEnergy(EmEdep);
            }
            // previously we added kinetic energy of neutrons entering the
            // fibre volumes.  this is commented out to avoid double counting
            // now that the birth-energy code (below) handles all neutron
            // contributions.
            /*
            if (particleName == "neutron") {
                if (track->GetCurrentStepNumber() == 1) {
                    G4double neutronEkin = track->GetKineticEnergy();
                    fEventAction->AddNeutronEkin(neutronEkin);
                }
            }
            */
                  
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

    // Count produced pions (secondaries only)
    if (step->GetTrack()->GetCurrentStepNumber() == 1 && step->GetTrack()->GetParentID() != 0) {
        G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();
        //if (particleName == "pi+" || particleName == "pi-" || particleName == "pi0") {
        if (particleName == "pi+" || particleName == "pi-" ) {
            fEventAction->AddPionCount();
        }
    }

    // Count produced neutrons (secondaries only)
    if (step->GetTrack()->GetCurrentStepNumber() == 1 && step->GetTrack()->GetParentID() != 0) {
        G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();
        if (particleName == "neutron") {
            fEventAction->AddNeutronCount();
        }
    }


    /*
    // Score kinetic energy of secondary neutrons created in this step.
    // This measures the energy "put into neutrons" by hadronic interactions
    // and uses the existing AddNeutronEkin accumulator for output.
    // Optionally one could filter by creator process or creation volume here.
    const auto* secondaries = step->GetSecondaryInCurrentStep();
    if (secondaries) {
        for (auto secIter = secondaries->begin(); secIter != secondaries->end(); ++secIter) {
            const G4Track* secTrack = *secIter;
            if (secTrack->GetDefinition()->GetPDGEncoding() == 2112) { // neutron PDG
                G4double ekin = secTrack->GetKineticEnergy();
                fEventAction->AddNeutronEkin(ekin);
                fEventAction->AddNeutronEnergy(ekin); // record individual neutron energy
                // print creator process:
                // const G4VProcess* proc = secTrack->GetCreatorProcess();
                // G4String procName = proc ? proc->GetProcessName() : "Unknown";
                // G4cout << "Secondary neutron created by " << procName
                //        << " Ekin=" << ekin << G4endl;
            }
        }
    }*/



    const auto* secondaries = step->GetSecondaryInCurrentStep();
    if (secondaries) {
        for (auto secIter = secondaries->begin(); secIter != secondaries->end(); ++secIter) {
            const G4Track* secTrack = *secIter;
            if (secTrack->GetDefinition()->GetPDGEncoding() == 2112) { // neutron
                G4double ekin = secTrack->GetKineticEnergy();
                //fEventAction->AddNeutronEkin(ekin);
                //fEventAction->AddNeutronEnergy(ekin);

                G4String modelName = secTrack->GetCreatorModelName();
                //G4double ekin = secTrack->GetKineticEnergy();
                //G4double time = secTrack->GetGlobalTime();

                G4String stage = "other";
                G4String tag   = "other";

                if (modelName == "model_G4EvaporationChannel" ||
                    modelName == "model_PRECO" ||
                    modelName == "model_G4FermiBreakUpVI" ||
                    modelName == "model_GammaNPreco") {
                    stage = "deexcitation";
                }
                else if (modelName == "model_BertiniCascade" ||
                        modelName == "model_FTFP") {
                    stage = "cascade";
                }
                else if (modelName == "model_hBertiniCaptureAtRest_NuclearCapture") {
                    stage = "capture";
                }

                if (stage == "deexcitation") {
                    tag = "evap_or_preco";
                }
                else if (stage == "cascade") {
                    if (ekin > 50*MeV) tag = "leading_like";
                    else               tag = "cascade_soft";
                }


                const G4VProcess* proc = secTrack->GetCreatorProcess();
                
                //G4String procName = proc ? proc->GetProcessName() : "primary/unknown";
                //G4int modelID = secTrack->GetCreatorModelID();
                //G4String modelName = secTrack->GetCreatorModelName();

                if (modelName == "model_FTFP") {

                G4cout << "Event " << G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()
                    //<< ": Neutron produced by process = " << procName
                    << ", model = " << modelName
                    << ", Ekin = " << ekin/MeV << " MeV"
                    << G4endl;
                }


                // add neutron energy if production model is not "model_FTFP" -> Esclude hard neutrons
                if (modelName != "model_FTFP") {
                    fEventAction->AddNeutronEkin(ekin);
                    fEventAction->AddNeutronEnergy(ekin);
                }







            }
        }
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
    G4int SiPMID = -1;
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
        
        G4int sipmID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
        //signalhit = fSignalHelper->ApplyPMTdishomogeneity(signalhit, sipmID);
        G4double distance_to_sipm = fSignalHelper->GetDistanceToSiPM(step);

        signalhit = fSignalHelper->AttenuateSSignal(signalhit, distance_to_sipm);
        fEventAction->AddVecSPMT( TowerID, signalhit ); 

        if(SiPMTower > -1)
        { 
            SiPMID = fDetConstruction->GetSiPMID(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1));
            fEventAction->AddVectorScin( SiPMTower*NoFibersTower + SiPMID , signalhit ); 
            
            
            //fEventAction->AddVectorScin( TowerID*NoFibersTower + SiPMID , signalhit ); \\ causing seg fault

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
                    G4int sipmID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);

                    // Apply tower dishomogeneity observed during TB24
                    //c_signal = fSignalHelper->ApplyPMTdishomogeneity(c_signal, sipmID);

                    // Attenuate Signal
                    c_signal = fSignalHelper->AttenuateCSignal(c_signal, distance_to_sipm);

                    G4int TowerID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(3);		
                    SiPMTower=fDetConstruction->GetSiPMTower(TowerID);

                    fEventAction->AddVecCPMT( TowerID, c_signal );
                    if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "leakbox"){G4cout << "Be careful: Signal is produced in leakbox!" << G4endl;}

                    if(SiPMTower > -1)
                    { // in sipm-readout tower
                        G4int SiPMID = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
                        //G4cout << step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << "\tSiPMID: " << SiPMID << "\tPhe: " << c_signal << "Tower: " << TowerID << G4endl;
                        //G4cout << "Hit fibre " << SiPMID << " in tower " << SiPMTower << G4endl;
                        fEventAction->AddVectorCher(SiPMTower*NoFibersTower+SiPMID, c_signal);
                        
                        
                        //fEventAction->AddVectorCher(TowerID*NoFibersTower+SiPMID, c_signal); // causing seg fault

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