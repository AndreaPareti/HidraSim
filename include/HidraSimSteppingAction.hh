//**************************************************
// \file HidraSimSteppingAction.hh
// \brief: Definition of HidraSimSteppingAction.hh
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) 
//          @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef HidraSimSteppingAction_h
#define HidraSimSteppingAction_h 1

//Includers from Geant4
//
#include "G4UserSteppingAction.hh"
#include "G4Types.hh"

//Forward declarations from Geant4
//
class G4OpBoundaryProcess;

//Forward declarations from project files
//
class HidraSimDetectorConstruction;
class HidraSimEventAction;

//Includers from project files
//
#include "HidraSimSignalHelper.hh"

class HidraSimSteppingAction : public G4UserSteppingAction {
    
    public:
        //Constructor
        //
        HidraSimSteppingAction(HidraSimEventAction* eventAction,
				const HidraSimDetectorConstruction* detConstruction,
                                const G4bool FullOptic );
        //De-constructor
        //
        virtual ~HidraSimSteppingAction();
        
        //User impementation of SteppingAction
        //
        virtual void UserSteppingAction( const G4Step* step );

        //Retrieve auxialiry info from Step
        //
        void AuxSteppingAction( const G4Step* step );

        //Fast signal simulation (no optical photon propagation)
        //fFullOptic == false
        //
        void FastSteppingAction( const G4Step* step ); 

        //Slow signal simulation (optical photon propagation)
        //fFullOptic == true
        //
        void SlowSteppingAction( const G4Step* step );
    
    private:

        HidraSimEventAction*  fEventAction;  

        G4OpBoundaryProcess* fOpProcess;
                
	//Pointer to HidraSimDetectorConstruction
	//
        const HidraSimDetectorConstruction* fDetConstruction;
				
	G4bool fFullOptic;

        //Pointer to only existing implementation (singleton)
    	//of HidraSimTowerHelper
    	//
        HidraSimSignalHelper* fSignalHelper;

};

#endif

//**************************************************
