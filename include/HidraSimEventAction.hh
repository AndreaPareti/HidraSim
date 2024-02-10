//**************************************************
// \file HidraSimEventAction.hh
// \brief: Definition of HidraSimEventAction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multiple times
//
#ifndef HidraSimEventAction_h
#define HidraSimEventAction_h 1

//Includers from Geant4
//
#include "G4UserEventAction.hh"
#include "globals.hh"

#include "HidraSimCalorimeterHit.hh"
#include "G4THitsMap.hh"
#include "HidraSimEventAction.hh"

//Includers from C++
//
#include <vector>
#include <set>
#include <tuple>


class HidraSimEventAction : public G4UserEventAction {
    
    public:
        //Constructor
        //
        HidraSimEventAction();
        //De-constructor
        //
        virtual ~HidraSimEventAction();

        virtual void  BeginOfEventAction(const G4Event* event);
        virtual void    EndOfEventAction(const G4Event* event);
    
        void AddScin(G4double de);//Add energy in scintillating fibers
        void AddCher(G4double de);//Add energy in Cherenkov fibers
        void Addenergy(G4double de);//Add energy depositedin calo
        void SavePrimaryPDGID(G4int pdgid);
        void SavePrimaryXY(G4double x, G4double y);
        void SaveAbsorberMaterial(G4String AbsorberMaterialName);
        void SavePrimaryEnergy(G4double primaryparticleenergy);
        void AddEscapedEnergy(G4double escapedenergy);
        void AddEscapedEnergyl(G4double escapedenergy);
        void AddEscapedEnergyd(G4double escapedenergy);
        void AddPSEnergy(G4double de);


        //Save vectors in ntuple
	//
        std::vector<G4double>& GetVectorSignals() {return VectorSignals;} 
        std::vector<G4double>& GetVectorSignalsCher() {return VectorSignalsCher;}
	    std::vector<G4double>& GetVecTowerE() {return VecTowerE;}
	    std::vector<G4double>& GetVecSPMT() {return VecSPMT;}
	    std::vector<G4double>& GetVecCPMT() {return VecCPMT;}


        //Fill vector of scintillating fibers with energy deposition
        //
        void AddVectorScin(G4double de, G4int fiber); 
        //Fill vector of cherenkov fibers with chernekov photoelectrons
        //
        void AddVectorCher(G4int fiber, G4int n);
        //Fill vector of energy in each tower
	    //
	    void AddVecTowerE(G4int TowerID, G4double de);
        //Fill vector of signals in scintillating PMTs
	    //
    	void AddVecSPMT(G4int PMTID, G4double de);
    	//Fill vector of signals in Cherenkov PMTs
        //
	    void AddVecCPMT(G4int PMTID, G4double de);
        //
        //
        // insert funtions per hits
        // Vector of photoelectrons in each hit (S fibers)
        std::vector<double>& GetHitPheSvector() {return fHitPheSvector;}
        std::vector<double>& GetHitZcoordSvector() {return fHitZcoordSvector;}
        std::vector<int>& GetHitSiPMIDSvector() {return fHitSiPMIDSvector;}
        std::vector<double>& GetHitPheCvector() {return fHitPheCvector;}
        std::vector<double>& GetHitZcoordCvector() {return fHitZcoordCvector;}
        std::vector<int>& GetHitSiPMIDCvector() {return fHitSiPMIDCvector;}


        void SaveHitPheSvector(std::vector<double> HitPheSvector);
        void AddNewHit();
        void InsertHitPheS(G4double phe);

    private:
        G4double  EnergyScin; //Energy in scintillating fibers
        G4double  EnergyCher; //Energy in Cherenkov fibers
        G4int     NofCherDet; //Number of Cherenkov p.e. detected 
    	G4int     NofScinDet; //Number of Scintillating p.e. detected
        G4double  EnergyTot;  //Total energy deposited (does not count invisibile energy)
        G4int     PrimaryPDGID; //PDGID of primary particle
        G4double  PrimaryParticleEnergy; //Primary particle energy
        G4double  PrimaryX; //Primary particle energy
        G4double  PrimaryY; //Primary particle energy
        G4double  EscapedEnergy; //Energy deposited in leakage absorber
        G4double  EscapedEnergyl; //Energy deposited in leakage absorber
        G4double  EscapedEnergyd; //Energy deposited in leakage absorber
	    G4double  PSEnergy;
        //Vector of SiPMs filled with scintillating signals
    	//
        std::vector<G4double> VectorSignals;
        //Vector of SiPMs filled with Cherenkov signals
	    //
        std::vector<G4double> VectorSignalsCher;
	    //Vector of PMTs filled with scintillating signals
    	//
    	std::vector<G4double> VecSPMT;
    	//Vector of PMTs filled with Cherenkov signals
    	//
	    std::vector<G4double> VecCPMT;
    	//Vector of energy deposited in towers
	    //
    	std::vector<G4double> VecTowerE;
        //
        // insert members per hits
        std::vector<double> fHitPheSvector;
        std::vector<double> fHitZcoordSvector;
        std::vector<int> fHitSiPMIDSvector; 
        //G4int NofHitSInEvt;
        std::vector<double> fHitPheCvector;
        std::vector<double> fHitZcoordCvector;
        std::vector<int> fHitSiPMIDCvector; 
        //G4int NofHitSInEvt;


        HidraSimCalorimeterHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;

        void PrintEventStatistics(
                              G4double SEdep, G4double SZcoord, G4int SSiPMID,
                              G4double CEdep, G4double CZcoord, G4int CSiPMID) const;

        G4int fSfiberHCID;
        G4int fCfiberHCID;

};

//Inline functions definition
//
inline void HidraSimEventAction::AddEscapedEnergy(G4double escapedenergy){
    EscapedEnergy += escapedenergy;
}
inline void HidraSimEventAction::AddEscapedEnergyl(G4double escapedenergy){
    EscapedEnergyl += escapedenergy;
}
inline void HidraSimEventAction::AddEscapedEnergyd(G4double escapedenergy){
    EscapedEnergyd += escapedenergy;
}

inline void HidraSimEventAction::SavePrimaryPDGID(G4int pdgid){
    PrimaryPDGID = pdgid;
}
inline void HidraSimEventAction::SavePrimaryXY(G4double x, G4double y){
    PrimaryX = x;
    PrimaryY = y;
}


inline void HidraSimEventAction::SavePrimaryEnergy(G4double primaryparticleenergy){
    PrimaryParticleEnergy = primaryparticleenergy;
}

inline void HidraSimEventAction::AddVectorScin(G4double de, G4int fiber) {
    VectorSignals.at(fiber) += de;
}

inline void HidraSimEventAction::AddVectorCher(G4int fiber, G4int n) {
    VectorSignalsCher.at(fiber) = VectorSignalsCher.at(fiber) + n;
}

inline void HidraSimEventAction::AddVecTowerE(G4int TowerID, G4double de) {
    VecTowerE.at(TowerID) += de;
}

inline void HidraSimEventAction::AddVecSPMT(G4int PMTID, G4double de) {
    VecSPMT.at(PMTID) += de;
}

inline void HidraSimEventAction::AddVecCPMT(G4int PMTID, G4double de) {
    VecCPMT.at(PMTID) += de;
}

inline void HidraSimEventAction::AddScin(G4double de){
    EnergyScin += de;
}

inline void HidraSimEventAction::AddCher(G4double de){
    EnergyCher += de;
}

inline void HidraSimEventAction::Addenergy(G4double de){
    EnergyTot += de;
}

inline void HidraSimEventAction::AddPSEnergy(G4double de){
    PSEnergy += de;
}



#endif

//**************************************************
