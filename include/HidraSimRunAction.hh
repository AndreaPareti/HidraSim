//**************************************************
// \file HidraSimRunAction.hh 
// \brief: Definition of HidraSimBRunAction class 
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef HidraSimRunAction_h
#define HidraSimRunAction_h 1

//Includers from Geant4
//
#include "G4UserRunAction.hh"
#include "globals.hh"

class HidraSimEventAction;
class G4Run;

class HidraSimRunAction : public G4UserRunAction {
    
    public:
        //Constructor
        //
        HidraSimRunAction( HidraSimEventAction* eventAction );
        //De-constructor
        //
        virtual ~HidraSimRunAction();

        //Methods
        //
        virtual void BeginOfRunAction(const G4Run*);
        virtual void EndOfRunAction(const G4Run*);

    private:
        HidraSimEventAction* fEventAction;

};

#endif

//**************************************************
