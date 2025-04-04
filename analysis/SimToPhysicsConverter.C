//**************************************************
// \file SimToPhysicsConverter.C
// \brief:  Converter from HidraSim (TB24 geo) to 2024 test beam formats
// \author: Andrea Pareti (Pavia Uni) 
//          based on previous converter scripts from G. Polesello
//          and L. Pezzotti      
// \start date: Nov. 13, 2024
//**************************************************
//
////usage: root -l -b -q 'SimToPhysicsConverter.C("filename")'
///  filename is the name of the sim output ntuple
//
//   It produces an histogram file 
//   physics+"filename"+.root
//
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <array>
#include <stdint.h>
#include <string>
#include <fstream>
#include <string>
#include <numeric>
#include <cstring>
// include file with geometry of module
#include "HidraConfig.h"

// test


struct PMTCalibration
{
  std::array<double, 1> PheGeVPS, PheGeVPC;
  PMTCalibration()
  {
    PheGeVPS[0] = SciPheGeV_Steel;
    PheGeVPC[0] = CerPheGeV_Steel;
  }
};


class EventOut
{
  public:
    EventOut(){};
    ~EventOut(){};
    uint32_t EventID;

    double TS00, TS10, TS11, TS12, TS13, TS14, TS15, TS16, TS17, TS20, TS21, TS22, TS23, TS24, TS25, TS30, TS31, TS32, TS33, TS34, TS35, TS40, TS41, TS42, TS43, TS44, TS45, TS50, TS51, TS52, TS53, TS54, TS55, TS60, TS61, TS62;
    double TC00, TC10, TC11, TC12, TC13, TC14, TC15, TC16, TC17, TC20, TC21, TC22, TC23, TC24, TC25, TC30, TC31, TC32, TC33, TC34, TC35, TC40, TC41, TC42, TC43, TC44, TC45, TC50, TC51, TC52, TC53, TC54, TC55, TC60, TC61, TC62;
    double L02, L03, L04, L05, L07, L08, L09, L10, L11, L12, L13, L14, L15, L16, L20;    
    //float SiPMPheC[160] = {0};
    //float SiPMPheS[160] = {0};
    float VecTowerE[36] = {0};
    float EnergyTot;
    //float totSiPMCene = 0.;
    //float totSiPMSene = 0.;
    //int NSiPMZero = 0.;
    float totPMTSene = 0.;
    float totPMTCene = 0.;
    float XDWC1, XDWC2, YDWC1, YDWC2;
    double PShower, TailCatcher, C1, C2;
    //int PShowerall;
    double ene_R0_S, ene_R1_S, ene_R2_S, ene_R3_S, ene_R4_S, ene_R5_S, ene_R6_S;
    double ene_R0_C, ene_R1_C, ene_R2_C, ene_R3_C, ene_R4_C, ene_R5_C, ene_R6_C;
    double totLeakage;

    void CompLeakage()
    {
	totLeakage = L02+L03+L04+L05+L07+L08+L09+L10+L11+L12+L13+L14+L15+L16+L20;
    }

    void CompSPMTene()
    {
      totPMTSene = TS00+TS10+TS11+TS12+TS13+TS14+TS15+TS16+TS17+TS20+TS21+TS22+TS23+TS24+TS25+TS30+TS31+TS32+TS33+TS34+TS35+TS40+TS41+TS42+TS43+TS44+TS45+TS50+TS51+TS52+TS53+TS54+TS55+TS60+TS61+TS62;
      ene_R0_S = TS00;
      ene_R1_S = TS16+TS15+TS14+TS17+TS13+TS10+TS11+TS12;
      ene_R2_S = TS20+TS21+TS22+TS23+TS24+TS25;
      ene_R3_S = TS30+TS31+TS32+TS33+TS34+TS35;
      ene_R4_S = TS40+TS41+TS42+TS43+TS44+TS45;
      ene_R5_S = TS50+TS51+TS52+TS53+TS54+TS55;
      ene_R6_S = TS60+TS61+TS62;
    }
    void CompCPMTene()
    {
      totPMTCene = TC00+TC10+TC11+TC12+TC13+TC14+TC15+TC16+TC17+TC20+TC21+TC22+TC23+TC24+TC25+TC30+TC31+TC32+TC33+TC34+TC35+TC40+TC41+TC42+TC43+TC44+TC45+TC50+TC51+TC52+TC53+TC54+TC55+TC60+TC61+TC62;
      ene_R0_C = TC00;
      ene_R1_C = TC16+TC15+TC14+TC17+TC13+TC10+TC11+TC12;
      ene_R2_C = TC20+TC21+TC22+TC23+TC24+TC25;
      ene_R3_C = TC30+TC31+TC32+TC33+TC34+TC35;
      ene_R4_C = TC40+TC41+TC42+TC43+TC44+TC45;
      ene_R5_C = TC50+TC51+TC52+TC53+TC54+TC55;
      ene_R6_C = TC60+TC61+TC62;

    }

    /*int SiPMCol(int index) { return index % 16; }
    int SiPMRow(int index) { return index / 16; }
    pair<double, double> SiPMSpos(int index)
    {
      int row = index / 16;
      int column = index % 16;
      double x = (column - 7) * 2 - 1.5;
      double y = 2. * sq3 * (4 - row) + sq3 / 2;
      return pair<double, double>(x, y);
    }
    pair<double, double> SiPMCpos(int index)
    {
      int row = index / 16;
      int column = index % 16;
      double x = (column - 7) * 2 - 0.5;
      double y = 2. * sq3 * (4 - row) + 1.5 * sq3;
      return pair<double, double>(x, y);
    }*/

};

class Event
{
  public:
  Event(){};
  ~Event(){};

  // Data Members
  double TS00, TS10, TS11, TS12, TS13, TS14, TS15, TS16, TS17, TS20, TS21, TS22, TS23, TS24, TS25, TS30, TS31, TS32, TS33, TS34, TS35, TS40, TS41, TS42, TS43, TS44, TS45, TS50, TS51, TS52, TS53, TS54, TS55, TS60, TS61, TS62;
  double TC00, TC10, TC11, TC12, TC13, TC14, TC15, TC16, TC17, TC20, TC21, TC22, TC23, TC24, TC25, TC30, TC31, TC32, TC33, TC34, TC35, TC40, TC41, TC42, TC43, TC44, TC45, TC50, TC51, TC52, TC53, TC54, TC55, TC60, TC61, TC62;
  double L02, L03, L04, L05, L07, L08, L09, L10, L11, L12, L13, L14, L15, L16, L20;
  double beamX, beamY;

  void calibratePMT(const PMTCalibration&, EventOut*);

};




void SimToPhysicsConverter(const string run){
  //Open ntuples
  string infile = infolder+run;
  //string infile = run;

  std::cout<<"Using file: "<<infile<<std::endl;
  char cinfile[infile.size() + 1];
  strcpy(cinfile, infile.c_str());

  string outfile = "physics_" + run;
  char coutfile[outfile.size() + 1];
  strcpy(coutfile, outfile.c_str());




  auto simfile = new TFile(cinfile, "READ");
  //auto *simtree = (TTree*)simfile->Get( "HidraSimout" );
  auto *simtree = (TTree*)simfile->Get( "HidraSimout" );

  // Create new tree and Event objects
  auto Outfile = new TFile(coutfile, "RECREATE");
  auto ftree = new TTree("Ftree", "Ftree");
  ftree->SetDirectory(Outfile);

  auto ev = new Event();
  auto evout = new EventOut();
  ftree->Branch("Events", evout);


  // Check entries in trees
  //
  std::cout << "Entries in PMT " << simtree->GetEntries() << std::endl;

// Allocate branch pointers
  int pdg;
  simtree->SetBranchAddress("PrimaryPDGID", &pdg);
  double venergy;
  simtree->SetBranchAddress("PrimaryParticleEnergy", &venergy);
  //double lenergy;
  //simtree->SetBranchAddress("EscapedEnergy", &lenergy);
  double edep;
  simtree->SetBranchAddress("EnergyTot", &edep);
  double Stot;
  simtree->SetBranchAddress("NofScinDet", &Stot);
  double Ctot;
  simtree->SetBranchAddress("NofCherDet", &Ctot);
  double PSdep;
  simtree->SetBranchAddress("PSEnergy", &PSdep);
  double beamX;
  simtree->SetBranchAddress("PrimaryX", &beamX);
  double beamY;
  simtree->SetBranchAddress("PrimaryY", &beamY);
  vector<double>* TowerE = NULL;
  simtree->SetBranchAddress("VecTowerE", &TowerE);
  vector<double>* SPMT = NULL;
  simtree->SetBranchAddress("VecSPMT", &SPMT);
  vector<double>* CPMT = NULL;
  simtree->SetBranchAddress("VecCPMT", &CPMT);
  vector<double>* SSiPM = NULL;
  simtree->SetBranchAddress("VectorSignals", &SSiPM);
  vector<double>* CSiPM = NULL;
  simtree->SetBranchAddress("VectorSignalsCher", &CSiPM);  
  vector<double>* LeakCounter = NULL;
  simtree->SetBranchAddress("VecLeakCounter", &LeakCounter);

  //double SciPheGeV = SciPheGeV_Steel;
  //double CerPheGeV = CerPheGeV_Steel;
   double SciPheGeV = 100;
   double CerPheGeV = 26;


for (unsigned int i = 0; i < simtree->GetEntries(); i++) {
    simtree->GetEntry(i);
    evout->EventID = i;
    //std::cout << SPMT->at(0) << std::endl;

    // Map S towers
    evout->TS60 = SPMT->at(0)/SciPheGeV;
    evout->TS61 = SPMT->at(1)/SciPheGeV;
    evout->TS62 = SPMT->at(2)/SciPheGeV;
    evout->TS50 = SPMT->at(3)/SciPheGeV;
    evout->TS51 = SPMT->at(4)/SciPheGeV;
    evout->TS52 = SPMT->at(5)/SciPheGeV;
    evout->TS40 = SPMT->at(6)/SciPheGeV;
    evout->TS41 = SPMT->at(7)/SciPheGeV;
    evout->TS42 = SPMT->at(8)/SciPheGeV;
    evout->TS30 = SPMT->at(9)/SciPheGeV;
    evout->TS31 = SPMT->at(10)/SciPheGeV;
    evout->TS32 = SPMT->at(11)/SciPheGeV;
    evout->TS20 = SPMT->at(12)/SciPheGeV;
    evout->TS21 = SPMT->at(13)/SciPheGeV;
    evout->TS22 = SPMT->at(14)/SciPheGeV;
    evout->TS10 = SPMT->at(15)/SciPheGeV;
    evout->TS11 = SPMT->at(16)/SciPheGeV;
    evout->TS12 = SPMT->at(17)/SciPheGeV;
    evout->TS17 = SPMT->at(18)/SciPheGeV;
    evout->TS00 = SPMT->at(19)/SciPheGeV;
    evout->TS13 = SPMT->at(20)/SciPheGeV;
    evout->TS16 = SPMT->at(21)/SciPheGeV;
    evout->TS15 = SPMT->at(22)/SciPheGeV;
    evout->TS14 = SPMT->at(23)/SciPheGeV;
    evout->TS25 = SPMT->at(24)/SciPheGeV;
    evout->TS24 = SPMT->at(25)/SciPheGeV;
    evout->TS23 = SPMT->at(26)/SciPheGeV;
    evout->TS35 = SPMT->at(27)/SciPheGeV;
    evout->TS34 = SPMT->at(28)/SciPheGeV;
    evout->TS33 = SPMT->at(29)/SciPheGeV;
    evout->TS45 = SPMT->at(30)/SciPheGeV;
    evout->TS44 = SPMT->at(31)/SciPheGeV;
    evout->TS43 = SPMT->at(32)/SciPheGeV;
    evout->TS55 = SPMT->at(33)/SciPheGeV;
    evout->TS54 = SPMT->at(34)/SciPheGeV;
    evout->TS53 = SPMT->at(35)/SciPheGeV;
    evout->CompSPMTene();

    // Map C towers
    evout->TC60 = CPMT->at(0)/CerPheGeV;
    evout->TC61 = CPMT->at(1)/CerPheGeV;
    evout->TC62 = CPMT->at(2)/CerPheGeV;
    evout->TC50 = CPMT->at(3)/CerPheGeV;
    evout->TC51 = CPMT->at(4)/CerPheGeV;
    evout->TC52 = CPMT->at(5)/CerPheGeV;
    evout->TC40 = CPMT->at(6)/CerPheGeV;
    evout->TC41 = CPMT->at(7)/CerPheGeV;
    evout->TC42 = CPMT->at(8)/CerPheGeV;
    evout->TC30 = CPMT->at(9)/CerPheGeV;
    evout->TC31 = CPMT->at(10)/CerPheGeV;
    evout->TC32 = CPMT->at(11)/CerPheGeV;
    evout->TC20 = CPMT->at(12)/CerPheGeV;
    evout->TC21 = CPMT->at(13)/CerPheGeV;
    evout->TC22 = CPMT->at(14)/CerPheGeV;
    evout->TC10 = CPMT->at(15)/CerPheGeV;
    evout->TC11 = CPMT->at(16)/CerPheGeV;
    evout->TC12 = CPMT->at(17)/CerPheGeV;
    evout->TC17 = CPMT->at(18)/CerPheGeV;
    evout->TC00 = CPMT->at(19)/CerPheGeV;
    evout->TC13 = CPMT->at(20)/CerPheGeV;
    evout->TC16 = CPMT->at(21)/CerPheGeV;
    evout->TC15 = CPMT->at(22)/CerPheGeV;
    evout->TC14 = CPMT->at(23)/CerPheGeV;
    evout->TC25 = CPMT->at(24)/CerPheGeV;
    evout->TC24 = CPMT->at(25)/CerPheGeV;
    evout->TC23 = CPMT->at(26)/CerPheGeV;
    evout->TC35 = CPMT->at(27)/CerPheGeV;
    evout->TC34 = CPMT->at(28)/CerPheGeV;
    evout->TC33 = CPMT->at(29)/CerPheGeV;
    evout->TC45 = CPMT->at(30)/CerPheGeV;
    evout->TC44 = CPMT->at(31)/CerPheGeV;
    evout->TC43 = CPMT->at(32)/CerPheGeV;
    evout->TC55 = CPMT->at(33)/CerPheGeV;
    evout->TC54 = CPMT->at(34)/CerPheGeV;
    evout->TC53 = CPMT->at(35)/CerPheGeV;
    evout->CompCPMTene();

    // Map ancillaries
    evout->XDWC1 = beamX; evout->XDWC2 = beamX;
    evout->YDWC1 = beamY;  evout->YDWC2 = beamY;
    evout->PShower = PSdep; 
    // Leakage counters
    // sim contains one extra volume (above calo, first layer) -> Not mapped
                                       evout->L07 = LeakCounter->at(4);   evout->L11 = LeakCounter->at(8);  evout->L15 = LeakCounter->at(12);  // Above calo 
    evout->L03 = LeakCounter->at(1);   evout->L08 = LeakCounter->at(5);   evout->L12 = LeakCounter->at(9);  evout->L16 = LeakCounter->at(13);  // calo right (as seen from beam)
    evout->L04 = LeakCounter->at(2);   evout->L09 = LeakCounter->at(6);   evout->L13 = LeakCounter->at(10); evout->L20 = LeakCounter->at(14);  // Below calo
    evout->L02 = LeakCounter->at(3);   evout->L05 = LeakCounter->at(7);   evout->L10 = LeakCounter->at(11); evout->L14 = LeakCounter->at(15);  // calo left (as seen from beam)

    // tail catcher
    evout->TailCatcher = LeakCounter->at(16);

    evout->CompLeakage();

    // Truth total deposited energy
    evout->EnergyTot = edep;

    ftree->Fill();
}




  ftree->Write();
  Outfile->Close();

}

//**************************************************
