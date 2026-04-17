//**************************************************
// \file HidraSimActionInitialization.hh
// \brief: Definition of HidraSimActionInitialization class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef HidraSimActionInitialization_h
#define HidraSimActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "G4Types.hh"

class HidraSimDetectorConstruction;
class G4GenericMessenger;

class HidraSimActionInitialization : public G4VUserActionInitialization {
public:
  explicit HidraSimActionInitialization(HidraSimDetectorConstruction*);
  ~HidraSimActionInitialization() override;

  void BuildForMaster() const override;
  void Build() const override;

private:
  G4bool fFullOptic = false;
  HidraSimDetectorConstruction* fDetConstruction = nullptr;
  G4GenericMessenger* fMessenger = nullptr;
};

#endif

//**************************************************
