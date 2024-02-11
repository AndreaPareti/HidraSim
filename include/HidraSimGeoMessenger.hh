//**************************************************
// \file HidraSimGeoMessenger.hh
// \brief: Definition of HidraSimGeomMessenger class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 10 August 2023
//**************************************************

#ifndef HidraSimGeoMessenger_h
#define HidraSimGeoMessenger_h 1

//Includers from Geant4
//
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//Includers from project files
//
class HidraSimDetectorConstruction;

class HidraSimGeoMessenger final : public G4UImessenger {

    public:
        //Constructor and destructor
        HidraSimGeoMessenger(HidraSimDetectorConstruction *DetConstruction);
        ~HidraSimGeoMessenger();
        
        //Virtual methods from base class
        void SetNewValue(G4UIcommand *command, G4String newValue) override;

    private:
        //Members
        HidraSimDetectorConstruction *fDetConstruction;
        G4UIdirectory *fMsgrDirectory;
        G4UIcmdWithADoubleAndUnit *fXshiftcmd;
        G4UIcmdWithADoubleAndUnit *fYshiftcmd;
        G4UIcmdWithADoubleAndUnit *fOrzrotcmd;
        G4UIcmdWithADoubleAndUnit *fVerrotcmd;

};

#endif //HidraSimGeoMessenger_h

//**************************************************
