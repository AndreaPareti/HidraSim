//**************************************************
// \file HidraAna.C
// \brief:  analysis skeleton for HidraSim ntuples
// \author: Giacomo Polesello (INFN Pavia) 
//          Edoardo Proserpio (Uni Insubria)
// \start date: May 3, 2022
//**************************************************
//
////usage: root -l -b -q 'HidraAna.C(energy,"filename")'
///  where energy is the energy of the beam, and filename
//   the name of the data ntuple
//
//   It produces an histogram file 
//   Hidra+"energy"+.root
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
#include <algorithm>
// include file with geometry of module
#include "HidraConfig.h"



void HidraAna(double energy, const string intup){
//Open ntuples
  //string infile = intup;
  string infile = infolder+intup;

  std::cout<<"Using file: "<<infile<<std::endl;
  char cinfile[infile.size() + 1];
  strcpy(cinfile, infile.c_str());
  auto simfile = new TFile(cinfile, "READ");
  auto *simtree = (TTree*)simfile->Get( "HidraSimout" );
  std::ostringstream os;
  os << energy;
  std::string enstr = os.str();
// 
  double sciPheGeV=0.;
  double cerPheGeV=0.;  
  double num_phe_sci=0;
  double num_phe_cev=0;
  double chi=0.;   
  double elcont=1.;
  double picont=1.;
  double sigma=0.4*sqrt(energy);

  std::string material;
  std::string filename = infile;
  std::string brass = "brass";
  std::string steel = "steel";

  bool isElectromagnetic = false;
  //if (filename.find("ele") != std::string::npos) {isElectromagnetic = true; std::cout << "\tElectron File\t" << std::endl;}

  if(filename.find(brass) != std::string::npos)
  {  
    std::cout << "\n USING BRASS ABSORBER" << std::endl; material=brass;
    sciPheGeV=SciPheGeV_Brass;	
    cerPheGeV=CerPheGeV_Brass;
    elcont=brass_elcontainment;
    chi=chi_brass;
    picont = EnergyContainmentConverter(energy, picont, brass_picontainment);   
  } 
  else  // Default material is Steel
  {  
    std::cout << "\nUSING STEEL ABSORBER" << std::endl;  material=steel;
    elcont=steel_elcontainment;
    sciPheGeV=SciPheGeV_Steel;
    cerPheGeV=CerPheGeV_Steel;
    chi=chi_steel;
    picont = EnergyContainmentConverter(energy, picont, steel_picontainment);    
  } 

  string outfile="hidra"+enstr+"_"+intup;
  TFile f(outfile.c_str(), "RECREATE");

//  build vectors with row and column position of each of the MiniModules
  int modcol[NofActiveModules];
  int modrow[NofActiveModules];
  //int nPMT = (sizeof(modflag)/sizeof(*modflag)); 
  for(int i=0;i<nPMT;i++){
    int row=i/NofmodulesX;
    int col=i%NofmodulesX;
    int imod=modflag[i];
    if(imod>=0){
      modcol[imod]=col;
      modrow[imod]=row;
    }
  }


// vector with row and column of each SiPM
  int NofFibersX = NofFibersrow*NofSiPMTowersX;
  int NofFibersY = NofFiberscolumn;
  int SiPMrow[NofSciSiPM];
  int SiPMcol[NofSciSiPM];

  for(int id=0; id<NofSciSiPM;id++)
  {
    int rowId=id/NofFibersrow;
    int colId=id%NofFibersrow;
  }


  // book histograms  
  double bmin=energy-sqrt(energy)*10.;
  //double bmax=energy+0.4*sqrt(energy)*10.;
  double bmax=energy+0.4*sqrt(energy)*10.;
  double bmin_val = (0 > bmin) ? 0 : bmin; // Get the maximum value

  int blow = int(bmin)+100;
  auto sciene = new TH1F("sciene", "Reco E_{S}; E [GeV]; Counts",100, bmin_val,bmax); 
  auto cerene = new TH1F("cerene", "Reco E_{C}; E [GeV]; Counts",100, bmin_val,bmax); 
  auto totene = new TH1F("totene", "Reco E (E_{S}+E_{C})/2; E [GeV]; Counts",100,bmin_val,bmax); 
  auto drene = new TH1F("drene", "Reco E_{DR} = S-#chiC / 1-#chi; E [GeV]; Counts",100,bmin_val,bmax);  
  auto totdep = new TH1F("totdep", "totdep",100,bmin_val,bmax);
  // Albedo represents beam energy - deposited E in the calo - truth leaked energy from side and back of calo  
  auto leakene = new TH1F("leakene", "leakene",100,0.,1.); 
  auto leakene_lat = new TH1F("leakene_lat", "leakene_lat", 100, 0., 1.);
  auto leakene_dep = new TH1F("leakene_dep", "leakene_dep", 100, 0., 1.);
  auto albedo = new TH1F("albedo", "albedo", 100, 0., 1.);
  auto chidist = new TH1F("chidist", "chidist",100,0.,1.);  
  auto histofem = new TH2F("histofem", "histofem", 100, 0., 1.5, 100, 0., 1.5);
  auto beamene = new TH1F("beamene", "beamene", 100, bmin_val, bmax);
  auto trueX = new TH1F("TrueX", "TrueX; X[mm]", 100, -15., 15.);
  auto trueY = new TH1F("TrueY", "TrueX; Y[mm]", 100, -15., 15.);
  auto S_recoX = new TH1F("recoX, S fiber", "recoX, S fiber; X[mm]", 100, -15., 15.);
  auto S_recoY = new TH1F("recoY, S fiber", "recoY, S fiber; Y[mm]", 100, -15., 15.);
  auto C_recoX = new TH1F("recoX, C fiber", "recoX, C fiber; X[mm]", 100, -15., 15.);
  auto C_recoY = new TH1F("recoY, C fiber", "recoY, C fiber; Y[mm]", 100, -15., 15.);
  auto scienec = new TH1F("scienec", "S energy (Corrected for pion containment)",100,bmin_val,bmax); 
  auto cerenec = new TH1F("cerenec", "C energy (Corrected for pion containment)",100,bmin_val,bmax);
  auto drenec = new TH1F("drenec", "Reco E_{DR} (corrected for pion containment); E [GeV]; Counts",100,bmin_val,bmax);  
  auto SCplot = new TH2F("SCplot", "S over C signals; S/E; C/E", 100, 0., 1.5, 100, 0., 1.5);
  auto TowerMapS = new TH2F("TowerMapS", "TowerMapS; Col; Row", NofmodulesX, 0, NofmodulesX, NofmodulesY, 0, NofmodulesY);
  auto TowerMapC = new TH2F("TowerMapC", "TowerMapC; Col; Row", NofmodulesX, 0, NofmodulesX, NofmodulesY, 0, NofmodulesY);
  auto SipmMapS = new TH2F("SipmMapS", "SipmS; Col; Row", NofSiPMTowersX*NofFiberscolumn, 0, NofSiPMTowersX*NofFiberscolumn, NofSiPMTowersY*NofFibersrow/2, 0, NofSiPMTowersY*NofFibersrow);
  auto SipmMapC = new TH2F("SipmMapC", "SipmC; Col; Row", NofSiPMTowersX*NofFiberscolumn, 0, NofSiPMTowersX*NofFiberscolumn, NofSiPMTowersY*NofFibersrow/2, 0, NofSiPMTowersY*NofFibersrow);

  auto LeakCounterSum = new TH1F("LeakCounterSum", "Energy deposit in lateral leakage counters; E [GeV]; ", 1000, 0., energy);
  auto TailCatcher = new TH1F("TailCatcherSum", "Energy deposit in tail catcher; E [GeV]", 1000, 0., energy);

  auto LeakProfile = new TProfile("LeakEdepProfile", "Signal in leakage counters Vs Truth deposited energy in Calo; Edep [GeV]; LeakageCounter sum [GeV]", 100, 0., energy, 0., energy);


  int nentries=simtree->GetEntries();
  std::cout<<"Entries "<<nentries<<std::endl;
  //Allocate branch pointers
  int pdg; simtree->SetBranchAddress( "PrimaryPDGID", &pdg ); 
  double venergy; simtree->SetBranchAddress( "PrimaryParticleEnergy", &venergy );
  double lenergy; simtree->SetBranchAddress( "EscapedEnergyl", &lenergy );
  double denergy; simtree->SetBranchAddress( "EscapedEnergyd", &denergy );
  double edep; simtree->SetBranchAddress( "EnergyTot", &edep );
  double Stot; simtree->SetBranchAddress( "NofScinDet", &Stot );
  double Ctot; simtree->SetBranchAddress( "NofCherDet", &Ctot );
  double PSdep; simtree->SetBranchAddress( "PSEnergy", &PSdep );
  double beamX; simtree->SetBranchAddress( "PrimaryX", &beamX );
  double beamY; simtree->SetBranchAddress( "PrimaryY", &beamY );
  vector<double>* TowerE = NULL;	simtree->SetBranchAddress( "VecTowerE", &TowerE );
  vector<double>* SPMT = NULL;	simtree->SetBranchAddress( "VecSPMT", &SPMT );
  vector<double>* CPMT = NULL;	simtree->SetBranchAddress( "VecCPMT", &CPMT );
  vector<double>* SSiPM = NULL;	simtree->SetBranchAddress( "VectorSignals", &SSiPM );
  vector<double>* CSiPM = NULL;	simtree->SetBranchAddress( "VectorSignalsCher", &CSiPM );
  vector<double>* LeakCounter = NULL;	simtree->SetBranchAddress( "VecLeakCounter", &LeakCounter );

  // Note that  SPMT, CPMT, SSiPM and CSiPM are given in photoelectrons
  // Other variables are given in MeV

  // Variables for the Event Loop
  double sum_towerE = 0;
  double sum_scin = 0;
  double sum_chev = 0; 
  double total_E = 0.;
  double beam_energy = 0.;
  double totleak=0.;
  double alb_energy=0.;
  double calib_sci=0.;
  double calib_cer=0.;
  



  // Loop on events 
  for( unsigned int i=0; i<simtree->GetEntries(); i++){
    //Check one event
    //int evtID = 105;  if(i<evtID) continue; if(i>evtID) break;    

    simtree->GetEntry(i);
    double totsci=0.;
    double totcer=0.;
    double tottow=0.;
    double partial_scin=0.;
    double partial_chev=0.;
    totleak = lenergy+denergy;
    double ecalo=energy-totleak/1000;  

    // Sum energy over all MiniModules
    for(unsigned int j=0; j<SPMT->size(); j++)
    {
      totsci+=SPMT->at(j)/sciPheGeV;
      totcer+=CPMT->at(j)/cerPheGeV;
      tottow+=TowerE->at(j);
      partial_scin+=SPMT->at(j);
      partial_chev+=CPMT->at(j);
      TowerMapS->Fill(modcol[j], modrow[j], SPMT->at(j)/sciPheGeV);
      TowerMapC->Fill(modcol[j], modrow[j], CPMT->at(j)/cerPheGeV);
    }
    //if(isElectromagnetic == false){
      //GetAttCorrection(energy);
    //  totsci = totsci/S_attenuation_correction;
    //  totcer = totcer/C_attenuation_correction;
    //}  

    for(unsigned int N=0; N<SSiPM->size(); N++){        // Loop over SiPMs - S Fibers
      double content = SSiPM->at(N)/sciPheGeV;
      unsigned int towID = static_cast<unsigned int>( N/(NofFiberscolumn*NofFibersrow/2) );
      unsigned int SiPMID = N%(NofFiberscolumn*NofFibersrow/2);
      unsigned int colID = static_cast<unsigned int>(SiPMID/(NofFibersrow/2));
      unsigned int rowID = 2*static_cast<unsigned int>(SiPMID%(NofFibersrow/2)); 
      SipmMapS->Fill( modcol[towID]*NofFiberscolumn + colID, modrow[towID]*NofFibersrow+rowID, content); 
    }


   for(unsigned int N=0; N<CSiPM->size(); N++){        // Loop over SiPMs - C Fibers
      double content = CSiPM->at(N)/cerPheGeV;
      unsigned int towID = static_cast<unsigned int>( N/(NofFiberscolumn*NofFibersrow/2) );
      unsigned int SiPMID = N%(NofFiberscolumn*NofFibersrow/2);
      unsigned int colID = static_cast<unsigned int>(SiPMID/(NofFibersrow/2));
      unsigned int rowID = 2*static_cast<unsigned int>(SiPMID%(NofFibersrow/2)) + 1; 
      SipmMapC->Fill( modcol[towID]*NofFiberscolumn + colID, modrow[towID]*NofFibersrow+rowID, content); 
    }    

    double tot_leakCount=0;
    // sum of leakage counters
    for(unsigned int N = 0; N<LeakCounter->size()-1; N++){
      tot_leakCount+=LeakCounter->at(N);

    }  
    LeakCounterSum->Fill(tot_leakCount/1000);
    TailCatcher->Fill(LeakCounter->at(16)/1000);

    alb_energy = (venergy-edep-totleak)/energy/1000;    
    sciene->Fill(totsci);    
    cerene->Fill(totcer); 
    double eneS = picont*totsci; double eneC = picont*totcer;
    scienec->Fill(eneS); cerenec->Fill(eneC);
    drene->Fill((totsci-chi*totcer)/(1-chi));


    totene->Fill(elcont*0.5*(totsci+totcer));   
    drenec->Fill(picont*(totsci-chi*totcer)/(1-chi));   
    totdep->Fill(tottow/1000.);   
    leakene_lat->Fill(lenergy/1000/energy);
    leakene_dep->Fill(denergy/1000/energy);
    leakene->Fill(totleak/1000/energy);  
    albedo->Fill(alb_energy);
    chidist->Fill((totsci-ecalo)/(totcer-ecalo)); 
    histofem->Fill(totsci/ecalo, totcer/ecalo, 1);
    sum_towerE+=tottow;
    sum_scin+=partial_scin;
    sum_chev+=partial_chev;
    total_E+=sum_towerE;
    double trueE = venergy/1000;
    beamene->Fill(trueE);
    trueX->Fill(beamX); trueY->Fill(beamY);

    SCplot->Fill( (eneS/trueE), (eneC/trueE), 1);


    int TowID, n, colID, rowID;
    double sciSiPM_totene=0;
    double cerSiPM_totene=0;
    double sci_Xweighted=0;
    double sci_Yweighted=0;
    double cer_Xweighted=0;
    double cer_Yweighted=0;
    double sciBarZ=0;
    double cerBarZ=0;
    double SiPM_X=1., SiPM_Y=0.;
    double SciBarX=0.;
    double SciBarY=0.;
    double CerBarX=0.;
    double CerBarY=0;

    LeakProfile->Fill(edep/1000, (tot_leakCount+LeakCounter->at(16))/1000 );


  //break;

  }  // Next event 


  
  std::cout << "\n\tLOOP ENDED \n" << std::endl;

  

  double leakedE = leakene->GetMean()*energy;
  double err_leakage = (leakene->GetRMS())/sqrt(nentries); 
  double mean_containment = energy - leakedE;
  double phe_GeV_sci = sum_scin/(energy-leakedE)/nentries;
  double phe_GeV_cer = sum_chev/(energy-leakedE)/nentries;
  std::cout << "\nphe/GeV (S): " << Stot/(energy-leakedE)/nentries << "\nphe/GeV (C): " << Ctot/(energy-leakedE)/nentries << std::endl;
  double chidist_mean = 0.;


  double min = energy - 1.5*sigma;
  /*
  //TF1 *fit1;
  if(isElectromagnetic){totene->Fit("gaus","Q","", bmin, 100*energy); fit1 = totene->GetFunction("gaus");}
  else{drenec->Fit("gaus","Q","", bmin, 10*energy); fit1 = drenec->GetFunction("gaus");}
 
  std::cout << "\n\tFIT OK\n " << std::endl;
  double peak1=fit1->GetParameter(1);
  double epeak1=fit1->GetParError(1);
  double rms1=fit1->GetParameter(2);
  double erms1=fit1->GetParError(2);*/


  // Fit total S energy
  TF1 *fit_sciene; double peak_sciene, epeak_sciene, rms_sciene, erms_sciene;
  TFitResultPtr is_sciene_fit = sciene->Fit("gaus","QS","", bmin_val, bmax); 
  if(is_sciene_fit->IsValid()){
    fit_sciene = sciene->GetFunction("gaus");
    peak_sciene = fit_sciene->GetParameter(1);
    epeak_sciene = fit_sciene->GetParError(1);
    rms_sciene = fit_sciene->GetParameter(2);
    erms_sciene = fit_sciene->GetParError(2);
  }

  // Fit total C energy
  TF1 *fit_cerene; double peak_cerene, epeak_cerene, rms_cerene, erms_cerene;
  TFitResultPtr is_cerene_fit = cerene->Fit("gaus","QS","", bmin_val, bmax); 
  if(is_cerene_fit->IsValid()){
    fit_cerene = cerene->GetFunction("gaus");
    peak_cerene = fit_cerene->GetParameter(1);
    epeak_cerene = fit_cerene->GetParError(1);
    rms_cerene = fit_cerene->GetParameter(2);
    erms_cerene = fit_cerene->GetParError(2);
  }
 
  // Fit combined energy (S+C)/2
  TF1 *fit_totene; double peak_totene, epeak_totene, rms_totene, erms_totene;
  TFitResultPtr is_totene_fit = totene->Fit("gaus","QS","", bmin_val, bmax); 
  if(is_totene_fit->IsValid()){
    fit_totene = totene->GetFunction("gaus");
    peak_totene = fit_totene->GetParameter(1);
    epeak_totene = fit_totene->GetParError(1);
    rms_totene = fit_totene->GetParameter(2);
    erms_totene = fit_totene->GetParError(2);
  }
 
  // Fit dual-readout 
  TF1 *fit_drene; double peak_drene, epeak_drene, rms_drene, erms_drene;
  TFitResultPtr is_drene_fit = drene->Fit("gaus","QS","", bmin_val, bmax); 
  if(is_drene_fit->IsValid()){
    fit_drene = drene->GetFunction("gaus");
    peak_drene = fit_drene->GetParameter(1);
    epeak_drene = fit_drene->GetParError(1);
    rms_drene = fit_drene->GetParameter(2);
    erms_drene = fit_drene->GetParError(2);
  }
 
  
  // Fit S, C and DR but corrected for estimated containment for pions
  TF1 *fit_scienec; double peak_scienec, epeak_scienec, rms_scienec, erms_scienec;
  TFitResultPtr is_scienec_fit = scienec->Fit("gaus","Q S","", bmin_val, bmax); 
  if(is_scienec_fit->IsValid()){
    fit_scienec = scienec->GetFunction("gaus");
    peak_scienec = fit_scienec->GetParameter(1);
    epeak_scienec = fit_scienec->GetParError(1);
    rms_scienec = fit_scienec->GetParameter(2);
    erms_scienec = fit_scienec->GetParError(2);
  }
 

  TF1 *fit_C; double peak_C, epeak_C, rms_C, erms_C;
  TFitResultPtr is_C_fit = cerenec->Fit("gaus","Q S","", bmin_val, bmax); 
  if(is_C_fit->IsValid()){
  fit_C = cerenec->GetFunction("gaus");
  peak_C = fit_C->GetParameter(1);
  epeak_C = fit_C->GetParError(1);
  rms_C = fit_C->GetParameter(2);
  erms_C = fit_C->GetParError(2);
  }


  TF1 *fit_drenec; double peak_drenec, epeak_drenec, rms_drenec, erms_drenec;
  TFitResultPtr is_drenec_fit = drenec->Fit("gaus","Q S","", bmin_val, bmax); 
  if(is_drenec_fit->IsValid()){
  fit_drenec = drenec->GetFunction("gaus");
  peak_drenec = fit_drenec->GetParameter(1);
  epeak_drenec = fit_drenec->GetParError(1);
  rms_drenec = fit_drenec->GetParameter(2);
  erms_drenec = fit_drenec->GetParError(2);
  }
 

  // PRINTOUT

  // Printout for Calibration: phe/GeV ratio for electrons, chi factor and containment
  //std::cout << " ? " << energy << "\t" << phe_GeV_sci << "\t" << phe_GeV_cer << "\t" << chidist_mean << "\t"  << energy/mean_containment << "\t" << leakene->GetMean() << "\t"  << leakene->GetRMS() << std::endl;
  // Print calorimeter parameters (energy)
  //cout << " # " << energy << " " << peak1 << " " << epeak1 << " " << rms1 << " " << erms1  << " "  << chi << " " << mean_containment/energy << " " << err_leakage/energy << endl;
 
  //std::cout << " % " << energy << " " << peak_drenec << " " << epeak_drenec << " " << rms_drenec << " " << erms_drenec << 
  //" " << peak_S << " " << epeak_S << " " << rms_S << " " << erms_S << " " << peak_C << " " << epeak_C << " " << rms_C << " " << erms_C << std::endl;


  f.Write();


}

//**************************************************
