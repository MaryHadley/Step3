#include "step2.cc"

void testStep2(){
  
  //TString inputFile="root://cmseos.fnal.gov//store/user/lpcljm/2018/LJMet94X_1lepTT_091518/nominal/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8_99.root";
  //TString inputFile="root://cmseos.fnal.gov//store/user/lpcljm/Mary_4Tops_Third_Real_Run_23Jan19_22h00Cern/nominal/4T_17_MC_Signal_List/4T_17_MC_Signal_List_2.root";
//TString inputFile="root://cmseos.fnal.gov//store/user/lpcljm/2018/LJMet94X_1lepTT_091518/nominal/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/TTTT_TuneCP5_13TeV-amcatnlo-pythia8_12.root";
//
  TString inputFile="/mnt/hadoop/users/mhadley/TTTT/LJMet94X_1lep_040819_step2haddsL1Prefiring/nominal/SingleMuon_31Mar18_hadd.root";
// TString inputFile="/mnt/hadoop/users/mhadley/TTTT/LJMet94X_1lep_040819_step2haddsL1Prefiring/nominal/SingleElectron_31Mar18_hadd.root";
  TString outputFile="test.root";
  
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  
  step2 t(inputFile,outputFile);
  t.Loop();
}


