// NOT IN ORIGINAL SIMULATION
// TRY TO SAVE HIT INFORMATION FOR TIMING STUDIES

#ifndef DREMTubesCalorimeterHit_h
#define DREMTubesCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"


/// Calorimeter hit class
/// It defines data members to store the photoelectrons associated to the hit and Z coordinate of the hit

class DREMTubesCalorimeterHit : public G4VHit
{
  public:
    DREMTubesCalorimeterHit();
    DREMTubesCalorimeterHit(const DREMTubesCalorimeterHit&);
    virtual ~DREMTubesCalorimeterHit();

    // operators
    const DREMTubesCalorimeterHit& operator=(const DREMTubesCalorimeterHit&);
    G4bool operator==(const DREMTubesCalorimeterHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void Add(G4double dphe);
    void SetSiPMID(G4int dSiPMID);
    void SetTowerID(G4int dTowerID);
    void SetZcoord(G4double dz);
    void SetPhe(G4double dPhe);
    void SetEdep(G4double DEdep);


    // get methods
    G4int GetSiPMID();
    G4int GetTowerID();
    G4double GetZcoord();
    G4double GetPhe();
    G4double GetEdep();
      
  private:
    G4int fSiPMID;        ///<ID of the hit fiber
    G4int fTowerID;        ///<ID of the tower whom hit fiber belong to
    G4double fZcoord;        ///< Z coordinate of the Hit
    G4double fPhe;           // Photoelectrons associated to the hit
    G4double fEdep;          // Truth energy deposited in the hit
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using DREMTubesCalorimeterHitsCollection = G4THitsCollection<DREMTubesCalorimeterHit>;
extern G4ThreadLocal G4Allocator<DREMTubesCalorimeterHit>* DREMTubesCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* DREMTubesCalorimeterHit::operator new(size_t)
{
  if (!DREMTubesCalorimeterHitAllocator) {
    DREMTubesCalorimeterHitAllocator = new G4Allocator<DREMTubesCalorimeterHit>;
  }
  void *hit;
  hit = (void *) DREMTubesCalorimeterHitAllocator->MallocSingle();
  return hit;
}

inline void DREMTubesCalorimeterHit::operator delete(void *hit)
{
  if (!DREMTubesCalorimeterHitAllocator) {
    DREMTubesCalorimeterHitAllocator = new G4Allocator<DREMTubesCalorimeterHit>;
  }
  DREMTubesCalorimeterHitAllocator->FreeSingle((DREMTubesCalorimeterHit*) hit);
}

inline void DREMTubesCalorimeterHit::Add(G4double dphe) {
  fPhe += dphe; 
}

// Define Getters

inline G4double DREMTubesCalorimeterHit::GetZcoord() { 
  return fZcoord; 
}

inline G4double DREMTubesCalorimeterHit::GetPhe() { 
  return fPhe; 
}

inline G4int DREMTubesCalorimeterHit::GetSiPMID() { 
  return fSiPMID; 
}

inline G4int DREMTubesCalorimeterHit::GetTowerID() {
  return fTowerID;
}

inline G4double DREMTubesCalorimeterHit::GetEdep() { 
  return fEdep; 
}

// Define Setters
inline void DREMTubesCalorimeterHit::SetSiPMID(G4int dSiPMID) {
  fSiPMID = dSiPMID;
}

inline void DREMTubesCalorimeterHit::SetTowerID(G4int dTowerID) {
  fTowerID = dTowerID;
}

inline void DREMTubesCalorimeterHit::SetZcoord(G4double dz) {
  fZcoord = dz;
}

inline void DREMTubesCalorimeterHit::SetPhe(G4double dPhe) {
  fPhe = dPhe;
}

inline void DREMTubesCalorimeterHit::SetEdep(G4double dEdep) {
  fEdep = dEdep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
