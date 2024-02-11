//**************************************************
// \file HidraSimDetectorConstruction.hh
// \brief: Definition of HidraSimDetectorConstruction class
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim) @lopezzot
// \start date: 7 July 2021
//**************************************************

//Prevent including header multplie times
//
#ifndef HidraSimDetectorConstruction_h
#define HidraSimDetectorConstruction_h 1

//Includers from Geant4
//
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"

//Include geometrical parameters
#include "HidraSimGeoPar.hh"
#include "HidraSimGeoMessenger.hh"
//Forward declaration
//
class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class HidraSimDetectorConstruction : public G4VUserDetectorConstruction {
  
    public:
        //Constructor
        //
        HidraSimDetectorConstruction();
        //De-constructor
        //
        virtual ~HidraSimDetectorConstruction();

    public:
        virtual G4VPhysicalVolume* Construct(); 			// returns DefineVolumes()  

        G4LogicalVolume* constructscinfiber(double tolerance,
                                            G4double tuberadius,
                                            G4double fiberZ,
                                            G4Material* absorberMaterial,
                                            G4double coreradius,
                                            G4double coreZ, 
                                            G4Material* ScinMaterial, 
                                            G4double claddingradiusmin,
                                            G4double claddingradiusmax,
                                            G4double claddingZ,
                                            G4Material* CherMaterial,
                                            G4String logicname,
                                            G4LogicalVolume* logic_Core_S_fiber,
                                            G4LogicalVolume* logic_Clad_S_fiber);
    
        G4LogicalVolume* constructcherfiber(double tolerance, 
                                            G4double tuberadius,
                                            G4double fiberZ,
                                            G4Material* absorberMaterial,
                                            G4double coreradius,
                                            G4double coreZ, 
                                            G4Material* CherMaterial, 
                                            G4double claddingradiusmin,
                                            G4double claddingradiusmax,
                                            G4double claddingZ,
                                            G4Material* CladCherMaterial,
                                            G4String logicname,
                                            G4LogicalVolume* logic_Core_C_fiber,
                                            G4LogicalVolume* logic_Clad_C_fiber);

        //Getters
    	//
	const G4VPhysicalVolume* GetLeakCntPV() const;
	const G4VPhysicalVolume* GetLeakCntlPV() const;
	const G4VPhysicalVolume* GetLeakCntdPV() const;    
    const G4VPhysicalVolume* GetWorldPV() const;

    G4double GetXshift() const {return fXshift;};
    G4double GetYshift() const {return fYshift;};
    G4double GetOrzrot() const {return fOrzrot;};
    G4double GetVerrot() const {return fVerrot;};

    //Setters
    //
    void SetXshift(const G4double& val) {fXshift=val;};
    void SetYshift(const G4double& val) {fYshift=val;};
    void SetOrzrot(const G4double& val) {fOrzrot=val;};
    void SetVerrot(const G4double& val) {fVerrot=val;};    

        //Other methods
	//
	G4int GetTowerID( const G4int& cpno ) const;
    G4int GetSiPMID(const G4int& cpno ) const; 
	G4int GetSiPMTower(const G4int& town ) const;
       
        //
	//  Build contour in x-y plane of a module as 
	//  an hexcell shape
	//
	std::vector<G4TwoVector> calcmod(double radius, int nrow, int ncol); 

    // make fibers sensitive detector to save hits
    virtual void ConstructSDandField();


    private:
        
        //Mandatory method for Geant4
        //
        G4VPhysicalVolume* DefineVolumes();

        // Members
        G4bool  fCheckOverlaps; // option for checking volumes overlaps
				
        G4VPhysicalVolume* fLeakCntPV; //PV: lekage counter
        G4VPhysicalVolume* fLeakCntlPV; //PV: lekage counter
        G4VPhysicalVolume* fLeakCntdPV; //PV: lekage counter
        G4VPhysicalVolume* fWorldPV;   //PV: wourld volume
 
        //Pointer to messenger for UI
        //
        HidraSimGeoMessenger* fGeoMessenger;

        //Parameters selectable via UI
        //
        G4double fXshift{0.}, fYshift{0.}, fVerrot{0.}, fOrzrot{0.};

};

inline G4int HidraSimDetectorConstruction::GetTowerID( const G4int& cpno ) const {
// remap as for 2021 hardware numbering from front face
//    const G4int idmap[9]={1,2,3,4,0,5,6,7,8};
// test:remap as for output of old simulation
//    const G4int idmap[9]={3,2,1,5,0,4,8,7,6};
//    return idmap[cpno-1];
    return cpno;		
}

inline G4int HidraSimDetectorConstruction::GetSiPMID( const G4int& cpno ) const {
// kept for compatibility with old simulation. Dummy for now
    return cpno;		
}

inline G4int HidraSimDetectorConstruction::GetSiPMTower( const G4int& town ) const {
    G4int SiPMTower=-1;
    for(int i=0;i<NofModulesSiPM;i++){
      if(town==SiPMMod[i])SiPMTower=i;
    }
    return SiPMTower;		
}

inline const G4VPhysicalVolume* HidraSimDetectorConstruction::GetLeakCntPV() const {
    return fLeakCntPV;
}
inline const G4VPhysicalVolume* HidraSimDetectorConstruction::GetLeakCntlPV() const {
    return fLeakCntlPV;
}
inline const G4VPhysicalVolume* HidraSimDetectorConstruction::GetLeakCntdPV() const {
    return fLeakCntdPV;
}

inline const G4VPhysicalVolume* HidraSimDetectorConstruction::GetWorldPV() const {
    return fWorldPV;
}

#endif

//**************************************************
