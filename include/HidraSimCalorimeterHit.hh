// NOT IN ORIGINAL SIMULATION
// TRY TO SAVE HIT INFORMATION FOR TIMING STUDIES

#ifndef HidraSimCalorimeterHit_h
#define HidraSimCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"


/// Calorimeter hit class
/// It defines data members to store the photoelectrons associated to the hit and Z coordinate of the hit

class HidraSimCalorimeterHit : public G4VHit
{
  public:
    HidraSimCalorimeterHit();
    HidraSimCalorimeterHit(const HidraSimCalorimeterHit&);
    virtual ~HidraSimCalorimeterHit();

    // operators
    const HidraSimCalorimeterHit& operator=(const HidraSimCalorimeterHit&);
    G4bool operator==(const HidraSimCalorimeterHit&) const;

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
    void AppendHit(G4double dPhe, G4double dZ);


    // get methods
    G4int GetSiPMID();
    G4int GetTowerID();
    G4double GetZcoord();
    G4double GetPhe();
    G4double GetEdep();
    std::vector<G4double> GetPheVec();
    std::vector<G4double> GetZVec();
      
  private:
    G4int fSiPMID;        ///<ID of the hit fiber
    G4int fTowerID;        ///<ID of the tower whom hit fiber belong to
    G4double fZcoord;        ///< Z coordinate of the Hit
    G4double fPhe;           // Photoelectrons associated to the hit
    G4double fEdep;          // Truth energy deposited in the hit
    std::vector<G4double> fPheVec; // vector of photoelectrons produced in fibers
    std::vector<G4double> fZVec; // vector of z positions at which photoelectrons are produced in the fiber
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using HidraSimCalorimeterHitsCollection = G4THitsCollection<HidraSimCalorimeterHit>;
extern G4ThreadLocal G4Allocator<HidraSimCalorimeterHit>* HidraSimCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* HidraSimCalorimeterHit::operator new(size_t)
{
  if (!HidraSimCalorimeterHitAllocator) {
    HidraSimCalorimeterHitAllocator = new G4Allocator<HidraSimCalorimeterHit>;
  }
  void *hit;
  hit = (void *) HidraSimCalorimeterHitAllocator->MallocSingle();
  return hit;
}

inline void HidraSimCalorimeterHit::operator delete(void *hit)
{
  if (!HidraSimCalorimeterHitAllocator) {
    HidraSimCalorimeterHitAllocator = new G4Allocator<HidraSimCalorimeterHit>;
  }
  HidraSimCalorimeterHitAllocator->FreeSingle((HidraSimCalorimeterHit*) hit);
}

inline void HidraSimCalorimeterHit::Add(G4double dphe) {
  fPhe += dphe; 
}

// Define Getters

inline G4double HidraSimCalorimeterHit::GetZcoord() { 
  return fZcoord; 
}

inline G4double HidraSimCalorimeterHit::GetPhe() { 
  return fPhe; 
}

inline G4int HidraSimCalorimeterHit::GetSiPMID() { 
  return fSiPMID; 
}

inline G4int HidraSimCalorimeterHit::GetTowerID() {
  return fTowerID;
}

inline G4double HidraSimCalorimeterHit::GetEdep() { 
  return fEdep; 
}

inline std::vector<G4double> HidraSimCalorimeterHit::GetPheVec(){
  return fPheVec;
}

inline std::vector<G4double> HidraSimCalorimeterHit::GetZVec(){
  return fZVec;
}

// Define Setters
inline void HidraSimCalorimeterHit::SetSiPMID(G4int dSiPMID) {
  fSiPMID = dSiPMID;
}

inline void HidraSimCalorimeterHit::SetTowerID(G4int dTowerID) {
  fTowerID = dTowerID;
}

inline void HidraSimCalorimeterHit::SetZcoord(G4double dz) {
  fZcoord = dz;
}

inline void HidraSimCalorimeterHit::SetPhe(G4double dPhe) {
  fPhe = dPhe;
}

inline void HidraSimCalorimeterHit::SetEdep(G4double dEdep) {
  fEdep = dEdep;
}

inline void HidraSimCalorimeterHit::AppendHit(G4double dPhe, G4double dZ){
  fPheVec.push_back(dPhe);
  fZVec.push_back(dZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
