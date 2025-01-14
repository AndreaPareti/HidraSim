//
/// \file HidraSimCalorimeterSD.hh
/// \brief Definition of the HidraSimCalorimeterSD class

#ifndef HidraSimCalorimeterSD_h
#define HidraSimCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "HidraSimCalorimeterHit.hh"
//Includers from project files
#include "HidraSimSignalHelper.hh"
#include <vector>
//#include "g4root.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Material.hh"
class G4Step;
class G4HCofThisEvent;
class HidraSimDetectorConstruction;

/// Calorimeter sensitive detector class
///
/// In Initialize(), it creates one hit for each calorimeter layer and one more
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class HidraSimCalorimeterSD : public G4VSensitiveDetector
{
  public:
    HidraSimCalorimeterSD(const G4String& name, 
                     const G4String& hitsCollectionName, 
                     G4int NofActiveFibers);
    virtual ~HidraSimCalorimeterSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);  


  private:
    HidraSimCalorimeterHitsCollection* fHitsCollection;

    G4int collectionID;
    G4int fNofActiveFibers;
    HidraSimSignalHelper* fSignalHelper;
    //HidraSimEventAction*  fEventAction;  
    const HidraSimDetectorConstruction* fDetConstruction;
    G4OpBoundaryProcess* fOpProcess;
  

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

