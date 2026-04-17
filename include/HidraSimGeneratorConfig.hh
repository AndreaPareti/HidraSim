// ============================================================================
// FILE: include/HidraSimGeneratorConfig.hh
// ============================================================================
#ifndef HidraSimGeneratorConfig_h
#define HidraSimGeneratorConfig_h 1

#include "globals.hh"
#include <vector>

class HidraSimGeneratorConfig {
public:
  static G4String& Mode() {
    static G4String mode = "gps";
    return mode;
  }

  // Store final state particle information from Pythia events
  static G4int& NumFinalStateParticles() {
    static G4int numParticles = 0;
    return numParticles;
  }

  static std::vector<G4double>& FinalStateEnergy() {
    static std::vector<G4double> energies;
    return energies;
  }

  static std::vector<G4double>& FinalStatePx() {
    static std::vector<G4double> px;
    return px;
  }

  static std::vector<G4double>& FinalStatePy() {
    static std::vector<G4double> py;
    return py;
  }

  static std::vector<G4double>& FinalStatePz() {
    static std::vector<G4double> pz;
    return pz;
  }

  static std::vector<G4int>& FinalStatePDGID() {
    static std::vector<G4int> pdgids;
    return pdgids;
  }

  // Clear all final state particle data
  static void ClearFinalStateData() {
    NumFinalStateParticles() = 0;
    FinalStateEnergy().clear();
    FinalStatePx().clear();
    FinalStatePy().clear();
    FinalStatePz().clear();
    FinalStatePDGID().clear();
  }
};

#endif