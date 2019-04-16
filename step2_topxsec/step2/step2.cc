#define step2_cxx
#include "step2.h"
//#include "GeneralFunctions.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <iostream>	// std::cout..
#include <algorithm>	// std::sort
#include <TRandom.h>
#include <TRandom3.h>
#include <sstream>
#include <string>
#include <vector>
#include "TMath.h"
#include <cmath>
#include <TMatrixD.h>
#include <TVectorT.h>
//#include "Davismt2.h"

using namespace std;

const double MTOP  = 173.5;
const double MW    = 80.4;
TLorentzVector Wlv;

TMatrixD SpheAplaTensor(const vector<TLorentzVector> allJets){
      TMatrixD MomentumTensor(3,3); 
      Double_t p2_sum=0.0;
      for(unsigned int njet = 0; njet < allJets.size(); ++njet){
        MomentumTensor(0, 0) += allJets[njet].Px()*allJets[njet].Px();
        MomentumTensor(0, 1) += allJets[njet].Px()*allJets[njet].Py();
        MomentumTensor(0, 2) += allJets[njet].Px()*allJets[njet].Pz();
        MomentumTensor(1, 0) += allJets[njet].Py()*allJets[njet].Px();
        MomentumTensor(1, 1) += allJets[njet].Py()*allJets[njet].Py();
        MomentumTensor(1, 2) += allJets[njet].Py()*allJets[njet].Pz();
        MomentumTensor(2, 0) += allJets[njet].Pz()*allJets[njet].Px();
        MomentumTensor(2, 1) += allJets[njet].Pz()*allJets[njet].Py();
        MomentumTensor(2, 2) += allJets[njet].Pz()*allJets[njet].Pz();
        p2_sum += allJets[njet].Px()*allJets[njet].Px()+allJets[njet].Py()*allJets[njet].Py()+allJets[njet].Pz()*allJets[njet].Pz();
      }
      if (p2_sum != 0){
          MomentumTensor(0, 0) /= p2_sum;
          MomentumTensor(0, 1) /= p2_sum;
          MomentumTensor(0, 2) /= p2_sum;
          MomentumTensor(1, 0) /= p2_sum;
          MomentumTensor(1, 1) /= p2_sum;
          MomentumTensor(1, 2) /= p2_sum;
          MomentumTensor(2, 0) /= p2_sum;
          MomentumTensor(2, 1) /= p2_sum;
          MomentumTensor(2, 2) /= p2_sum;
      }
      TVectorD *_pv = new TVectorD(3);
      return MomentumTensor;
}

bool compareFunc(const pair<TLorentzVector, double> &a, const pair<TLorentzVector, double> &b){
    return (a.second > b.second);
}

bool comparepair( const std::pair<double,int> a, const std::pair<double,int> b) { return a.first > b.first; }




void step2::Loop()
{   
   if (inputTree == 0) return;
   outputFile->cd();
   TTree *outputTree = inputTree->CloneTree(); //Copy of Input Tree
//    TTree *outputTree = new TTree("ljmet","ljmet"); //No Copy of Input Tree   
                              
   TBranch *b_atLeast3BJets   = outputTree->Branch("atLeast3BJets", &atLeast3BJets, "atLeast3BJets/I");
   TBranch *b_atLeast2BJets   = outputTree->Branch("atLeast2BJets", &atLeast2BJets, "atLeast2BJets/I");
   TBranch *b_atLeast1BJet    = outputTree->Branch("atLeast1BJet", &atLeast1BJet,   "atLeast1BJet/I");
   TBranch *b_atMost5Jets     = outputTree->Branch("atMost5Jets", &atMost5Jets,   "atMost5Jets/I");
   TBranch *b_v_DCSV_allJets_DCSV_Ordered = outputTree->Branch("v_DCSV_allJets_DCSV_Ordered", &v_DCSV_allJets_DCSV_Ordered);


   //Define histograms //
   TH2D *h_den_eta_vs_pt_atMost5jtot_atLeast2bj = new TH2D("h_den_eta_vs_pt_atMost5jtot_atLeast2bj", "Eta vs. pT of add. jets (min 2 b tags in sample); Eta; pT [GeV]", 60, -3., 3., 25, 0, 500);
   TH2D *h_num_eta_vs_pt_atMost5jtot_atLeast2bj = new TH2D("h_num_eta_vs_pt_atMost5jtot_atLeast2bj", "Eta vs. pT of add. b tagged jets (min 2 b tags in sample); Eta; pT [GeV]", 60, -3., 3., 25, 0, 500);
   TH2D *h_den_eta_vs_pt_atMost5jtot_atLeast3bj = new TH2D("h_den_eta_vs_pt_atMost5jtot_atLeast3bj", "Eta vs. pT of add. jets (min 3 b tags in sample); Eta; pT [GeV]", 60, -3., 3., 25, 0, 500);
   TH2D *h_num_eta_vs_pt_atMost5jtot_atLeast3bj = new TH2D("h_num_eta_vs_pt_atMost5jtot_atLeast3bj", "Eta vs. pT of add. b tagged jets (min 3 b tags in sample); Eta; pT [GeV] ", 60, -3., 3., 25, 0, 500);
   
   ///// end definition of histograms /////
   
   Long64_t nentries = inputTree->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
   int tot_atLeast1BJet = 0; //leaving these branches so I do not have to change the header file, ok I guess I could leave them anyway, but I am just changing as few things as possible and I might as well fill them, it is info that does not take any time to get, in the scheme of things
   int tot_atLeast2BJets = 0; //same comment as above
   int tot_atLeast3BJets =0; //same comment as above
   int tot_atMost5Jets = 0; //same comment as above
   

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //std::cout<<"jentry : "<<jentry<<std::endl;
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = inputTree->GetEntry(jentry);   nbytes += nb;
     if (Cut(ientry) != 1) continue;
     if(jentry % 1000 ==0) std::cout<<"Completed "<<jentry<<" out of "<<nentries<<" events"<<std::endl; 
//     if(jentry > 10) break;
     atLeast1BJet = 0; //0 aka false
     atLeast2BJets = 0; //0 aka False
     atLeast3BJets = 0; //0 aka False
     atMost5Jets = 0; //0 aka False
if (NJets_JetSubCalc <= 5){    
     tot_atMost5Jets +=1;
     atMost5Jets = 1; //aka True
     
     if (NJetsCSV_JetSubCalc >= 1){ 
         atLeast1BJet = 1; //aka True
         tot_atLeast1BJet +=1;
         } 
         
     if (NJetsCSV_JetSubCalc >= 2){ 
         atLeast2BJets = 1; //aka True
         tot_atLeast2BJets +=1;
         vector<pair<double,int>> jetcsvindpair;
          for(unsigned int ijet=0; ijet < theJetPt_JetSubCalc_PtOrdered->size(); ijet++){
	jetcsvindpair.push_back(std::make_pair(theJetDeepCSVb_JetSubCalc_PtOrdered->at(ijet)+theJetDeepCSVbb_JetSubCalc_PtOrdered->at(ijet),ijet));
      }
    std::sort(jetcsvindpair.begin(), jetcsvindpair.end(), comparepair);
   // h_den_eta_vs_pt_atMost5jtot_atLeast2bj->Fill(theJetEta_JetSubCalc_PtOrdered->at(jetcsvindpair[2].second), theJetPt_JetSubCalc_PtOrdered->at(jetcsvindpair[2].second)); //this takes us to the third entry in the jetcsvindpair vector (aka [2]) and then gets us the second thing in the pair, aka the index //example from sinan, just leaving it here to remember the idea
  // std::cout<<"jetcsvindpair.size() is: "<<jetcsvindpair.size()<<std::endl;
     for (unsigned int pjet = 2; pjet < jetcsvindpair.size(); pjet++){ //2 because I want to exclude the first two highest CSV scored jets and start at the third
    //  std::cout<<jetcsvindpair[pjet].first<<std::endl; 
        h_den_eta_vs_pt_atMost5jtot_atLeast2bj->Fill(theJetEta_JetSubCalc_PtOrdered->at(jetcsvindpair[pjet].second), theJetPt_JetSubCalc_PtOrdered->at(jetcsvindpair[pjet].second));
        if (jetcsvindpair[pjet].first > 0.4941){
         h_num_eta_vs_pt_atMost5jtot_atLeast2bj->Fill(theJetEta_JetSubCalc_PtOrdered->at(jetcsvindpair[pjet].second), theJetPt_JetSubCalc_PtOrdered->at(jetcsvindpair[pjet].second));
} //close the if the jet is b tagged loop
}//close the loop over pjet
        }//closing the NJetsCSV_JetSubCalc >= 2 loop
         
         
             

     if (NJetsCSV_JetSubCalc >=3){
         atLeast3BJets =1; //aka True
         tot_atLeast3BJets += 1;
         vector<pair<double,int>> jetcsvindpair_3; //3 for min 3 b jets
     for(unsigned int ijet=0; ijet < theJetPt_JetSubCalc_PtOrdered->size(); ijet++){ //iterators defined in for loops will get forgotten outside the for loop
	jetcsvindpair_3.push_back(std::make_pair(theJetDeepCSVb_JetSubCalc_PtOrdered->at(ijet)+theJetDeepCSVbb_JetSubCalc_PtOrdered->at(ijet),ijet));
      }
      std::sort(jetcsvindpair_3.begin(), jetcsvindpair_3.end(), comparepair);
      for (unsigned int pjet = 3; pjet < jetcsvindpair_3.size(); pjet++){ //start at pjet =3 because I want to exclude the first two highest CSV scored jets and start at the fourth
      h_den_eta_vs_pt_atMost5jtot_atLeast3bj->Fill(theJetEta_JetSubCalc_PtOrdered->at(jetcsvindpair_3[pjet].second), theJetPt_JetSubCalc_PtOrdered->at(jetcsvindpair_3[pjet].second));
      if (jetcsvindpair_3[pjet].first > 0.4941){
         h_num_eta_vs_pt_atMost5jtot_atLeast3bj->Fill(theJetEta_JetSubCalc_PtOrdered->at(jetcsvindpair_3[pjet].second), theJetPt_JetSubCalc_PtOrdered->at(jetcsvindpair_3[pjet].second));
      } // close the if it met  b tag requirement loop
      }//close the pjet loop
      
         }//closing the NJetsCSV_JetSubCalc >=3 loop 
          
}


    b_atLeast1BJet->Fill();
    b_atLeast2BJets->Fill(); 
    b_atLeast3BJets->Fill(); 
    b_atMost5Jets->Fill();
    b_v_DCSV_allJets_DCSV_Ordered->Fill(); 		  
   }



std::cout<<"DONE "<<nentries<<std::endl;   


h_den_eta_vs_pt_atMost5jtot_atLeast2bj->Write();
h_num_eta_vs_pt_atMost5jtot_atLeast2bj->Write();
h_den_eta_vs_pt_atMost5jtot_atLeast3bj->Write();
h_num_eta_vs_pt_atMost5jtot_atLeast3bj->Write();
outputFile->Write();
delete outputFile;
delete inputFile;
}
