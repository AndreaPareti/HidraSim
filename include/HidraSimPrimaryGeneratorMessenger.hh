#ifndef HidraSimPrimaryGeneratorMessenger_h
#define HidraSimPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class HidraSimPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcommand;

class HidraSimPrimaryGeneratorMessenger : public G4UImessenger {
public:
  explicit HidraSimPrimaryGeneratorMessenger(HidraSimPrimaryGeneratorAction* gen);
  ~HidraSimPrimaryGeneratorMessenger() override;

  void SetNewValue(G4UIcommand* cmd, G4String newValue) override;

private:
  HidraSimPrimaryGeneratorAction* fGen = nullptr;
  G4UIdirectory* fDir = nullptr;
  G4UIcmdWithAString* fModeCmd = nullptr;
};

#endif