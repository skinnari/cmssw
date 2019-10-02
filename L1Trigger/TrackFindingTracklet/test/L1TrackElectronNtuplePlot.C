// ----------------------------------------------------------------------------------------------------------------
// Basic example ROOT script for making tracking performance plots using the ntuples produced by L1TrackNtupleMaker.cc
//
//    e.g.  .x  L1TrackNtuplePlot.C("TTbar_PU200_hybrid")
// 
// By Louise Skinnari, June 2013 <==I was in in high school >_<
//
// I essentially stealed Halil's function and inserted it into "L1TrackNtuplePlot.C"
//
// e.g.  .x L1TrackElectronNtuplePlot.C("D35_single_electron")
//
// it also takes several bool arguments
//
// e.g.  .x L1TrackElectronNtuplePlot.C("D35_single_electron",true,true,false,0)
//
// 1st argument: turn on tracks correction (default yes)
// 2nd argument: turn on elliptical matching (default yes)
// 3rd argument: turn on tracks isolation (default no)
// 4th argument: barrel = 0, endcap = 1, barrel+endcap = else (default 0)
//
//
// ----------------------------------------------------------------------------------------------------------------

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include <TError.h>
#include "TControlBar.h"

#include <ROOT/RVec.hxx>
#include "TMath.h"
#include "TLorentzVector.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text);

float phiCorr (float &pt,float &eta,float &phi,int &charge,int &hgc,float phicorrection); // phi correction
float etaCorr (float &eta,float &z0,int &hgc); // eta correction
vector<float> dZAtDRMatchPtVMaker ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<int> &mainlooseTkID, vector<int> &mainhgc, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, vector<float> &loopz0, float dr, bool phietacorr, float correction=0, float trackZ0multiplier=1.0 );
vector<float> dLAtDRMatchPtVMaker ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<int> &mainlooseTkID, vector<int> &mainhgc, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, vector<float> &loopz0, float dr, bool phietacorr, float correction=0, float trackZ0multiplier=1.0 );
vector<float> dPhiAtDRMatchPtVMaker ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<int> &mainlooseTkID, vector<int> &mainhgc, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, vector<float> &loopz0, vector<float> &loopd0, float dr,bool phietacorr, bool dphimatch, bool iso);
vector<float> dEtaAtDRMatchPtVMaker ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<int> &mainlooseTkID, vector<int> &mainhgc, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, vector<float> &loopz0, float dr, bool phietacorr, bool dphimatch, bool iso);
vector<float> dPhiTrack ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, float dr);
vector<float> dEgamma ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi,vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &looplooseTkID, vector<int> &loophgc, vector<float> &loopdphi, vector<float> &loopdeta, vector<float> &loopdz, vector<float> &loopdl, float dr, bool tk_eg_match, bool gen_eg_delta, int delta);
vector<float> dR ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi);
float rel_isolation (float &mainpt, float &maineta, float &mainphi, float &mainz0, float &maind0, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<float> &loopz0, vector<float> &loopd0 );

// ----------------------------------------------------------------------------------------------------------------
// Main script
// ----------------------------------------------------------------------------------------------------------------


void L1TrackElectronNtuplePlot(TString type, bool ifcorr=true, bool dphimatching=false, bool ifiso=false, int isbarrel=0,TString treeName="",float TP_maxEta=2.4,float TP_minPt=2.0) {

 

  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;
  
  SetPlotStyle();
  
  
  // ----------------------------------------------------------------------------------------------------------------
  // define input options

  


  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  TChain* tree = new TChain("L1TkElNtuple"+treeName+"/eventTree");
  tree->Add(type+".root");
  
  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }
  

  // ----------------------------------------------------------------------------------------------------------------
  // define leafs & branches
  
  // all L1 tracks
  vector<float>* gen_pt;
  vector<float>* gen_eta;
  vector<float>* gen_phi;
    
  // *L1 cluster
  vector<float>* elec_et;
  vector<float>* elec_phi;
  vector<float>* elec_eta;
  vector<int>* elec_looseTkID;
  vector<int>* elec_hgc;


  // all L1 tracks
  vector<float>* trk_pt;
  vector<float>* trk_eta;
  vector<float>* trk_phi;
  vector<int>* trk_charge;
  vector<float>* trk_z0;
  vector<float>* trk_d0;
 
  TBranch* b_gen_pt;
  TBranch* b_gen_eta;
  TBranch* b_gen_phi;

  TBranch* b_elec_et;
  TBranch* b_elec_phi;
  TBranch* b_elec_eta;
  TBranch* b_elec_looseTkID;
  TBranch* b_elec_hgc;


  TBranch* b_trk_pt;
  TBranch* b_trk_eta;
  TBranch* b_trk_phi;
  TBranch* b_trk_z0;
  TBranch* b_trk_charge;
  TBranch* b_trk_d0;
    
  gen_pt=0;
  gen_eta=0;
  gen_phi=0;
  
  elec_et=0;
  elec_phi=0;
  elec_eta=0;
  elec_looseTkID=0;
  elec_hgc=0;
    
    
  trk_pt=0;
  trk_eta=0;
  trk_phi=0;
  trk_charge=0;
  trk_z0=0;
  trk_d0=0;
    
    
  tree->SetBranchAddress("gen_pt",   &gen_pt,   &b_gen_pt);
  tree->SetBranchAddress("gen_phi",  &gen_phi,  &b_gen_phi);
  tree->SetBranchAddress("gen_eta",  &gen_eta,  &b_gen_eta);

  tree->SetBranchAddress("elec_et",     &elec_et,     &b_elec_et);
  tree->SetBranchAddress("elec_phi",    &elec_phi,    &b_elec_phi);
  tree->SetBranchAddress("elec_eta",    &elec_eta,    &b_elec_eta);
  tree->SetBranchAddress("elec_looseTkID", &elec_looseTkID, &b_elec_looseTkID);
  tree->SetBranchAddress("elec_hgc",    &elec_hgc,    &b_elec_hgc);


  tree->SetBranchAddress("trk_pt",   &trk_pt,   &b_trk_pt);
  tree->SetBranchAddress("trk_phi",  &trk_phi,  &b_trk_phi);
  tree->SetBranchAddress("trk_eta",  &trk_eta,  &b_trk_eta);
  tree->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
  tree->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
  tree->SetBranchAddress("trk_charge",   &trk_charge,   &b_trk_charge);
    
    
  vector<float> GenP_pt;
  vector<float> GenP_eta;
  vector<float> GenP_dR;
  vector<float> eg_dphitrk;
  vector<float> eg_dphieg;
  vector<float> eg_dphitkeg;
  vector<float> EG_dPhi;
  vector<float> EG_dEta;
  vector<float> EG_dZ;
  vector<float> EG_dL;
  vector<float> EG_et;
  vector<float> EG_eta;
  vector<float> EG_hgc;
    
  int count_all=0;
  int count_all_H=0;
  int count_all_L=0;
    
  int count_trk=0;
  int count_trk_H=0;
  int count_trk_L=0;
    
  int count_eg=0;
  int count_eg_H=0;
  int count_eg_L=0;
    
  int count_egtrk=0;
  int count_egtrk_H=0;
  int count_egtrk_L=0;
  

  // ----------------------------------------------------------------------------------------------------------------
  // histograms
  // ----------------------------------------------------------------------------------------------------------------

  /////////////////////////////////////////////////
  // NOTATION:                                   //
  // 'C' - Central eta range, |eta|<0.8          //
  // 'I' - Intermediate eta range, 0.8<|eta|<1.6 //
  // 'F' - Forward eta range, |eta|>1.6          //
  //                                             //
  // 'L' - Low pt range,  pt<8 GeV               //
  // 'H' - High pt range, pt>8 GeV               //
  /////////////////////////////////////////////////



  // ----------------------------------------------------------------------------------------------------------------
  // for efficiencies

  TH1F* h_gen_pt    = new TH1F("gen_pt",   ";Gen Particle p_{T} [GeV]; Gen Particle / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_gen_eta   = new TH1F("gen_eta",  ";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  TH1F* h_gen_eta_L = new TH1F("gen_eta_L",";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  TH1F* h_gen_eta_H = new TH1F("gen_eta_H",";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
    
  TH1F* h_match_trk_pt    = new TH1F("match_trk_pt",   ";Gen Particle p_{T} [GeV]; Gen Particle / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_match_trk_eta   = new TH1F("match_trk_eta",  ";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  TH1F* h_match_trk_eta_L = new TH1F("match_trk_eta_L",";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  TH1F* h_match_trk_eta_H = new TH1F("match_trk_eta_H",";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  
  TH1F* h_match_eg_pt    = new TH1F("match_eg_pt",   ";Gen Particle p_{T} [GeV]; Gen Particle / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_match_eg_eta   = new TH1F("match_eg_eta",  ";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  TH1F* h_match_eg_eta_L = new TH1F("match_eg_eta_L",";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  TH1F* h_match_eg_eta_H = new TH1F("match_eg_eta_H",";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  
  TH1F* h_EG_et    = new TH1F("EG_et",   ";L1EG et [GeV]; L1EG", 120,  0,   120.0);
  TH1F* h_EG_eta_BB    = new TH1F("EG_eta_BB",   ";L1EG #eta; L1EG", 50,  -1*TP_maxEta,   TP_maxEta);
  TH1F* h_EG_et_BB    = new TH1F("EG_et_BB",   ";L1EG et [GeV]; L1EG", 100,  0,   100.0);
  TH1F* h_matched_EG_et_BB  = new TH1F("matched_EG_et_BB",   ";L1EG et [GeV]; L1EG", 100,  0,   100.0);
  TH1F* h_matched_EG_eta_BB    = new TH1F("matched_EG_eta_BB",   ";L1EG #eta; L1EG", 50,  -1*TP_maxEta,   TP_maxEta);
  TH1F* h_EG_et_BE    = new TH1F("EG_et_BE",   ";L1EG et [GeV]; L1EG", 100,  0,   100.0);
  TH1F* h_EG_eta   = new TH1F("EG_eta",  ";L1EG #eta; L1EG / 0.1",             50, -2.5,   2.5);
  
    
  TH1F* h_match_egtrk_pt    = new TH1F("match_egtrk_pt",   ";Gen Particle p_{T} [GeV]; Gen Particle / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_match_egtrk_eta   = new TH1F("match_egtrk_eta",  ";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  TH1F* h_match_egtrk_eta_L = new TH1F("match_egtrk_eta_L",";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
  TH1F* h_match_egtrk_eta_H = new TH1F("match_egtrk_eta_H",";Gen Particle #eta; Gen Particle / 0.1",             50, -2.5,   2.5);
    
  TH1F* h_res_dR_1 = new TH1F("res_dR_1", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_2 = new TH1F("res_dR_2", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_3 = new TH1F("res_dR_3", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_4 = new TH1F("res_dR_4", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_5 = new TH1F("res_dR_5", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_6 = new TH1F("res_dR_6", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_7 = new TH1F("res_dR_7", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_8 = new TH1F("res_dR_8", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_9 = new TH1F("res_dR_9", ";Gen Particle p_{T} [GeV]; Events", 100,  -0.1,   0.1);
  TH1F* h_res_dR_10 = new TH1F("res_dR_10", ";Gen Particle p_{T} [GeV]s; Events",100,  -0.1,   0.1);
    
  TH1F* h_res_dR = new TH1F("res_dR", ";dR; # of GenParticles",100,  -0.1,   0.1);
    
  TH1F* h_EG_dPhi = new TH1F("EG_dPhi", ";dPhi; # of L1EG",100,  -0.1,   0.1);
  TH1F* h_EG_dEta = new TH1F("EG_dEta", ";dEta; # of L1EG",100,  -0.1,   0.1);
  TH1F* h_EG_dZ = new TH1F("EG_dZ", ";dZ [cm]; # of L1EG",100,  -15,   15);
  TH1F* h_EG_dL = new TH1F("EG_dL", ";dL [cm]; # of L1EG",100,  -5,   5);

  TH2F * h_EG_dPhivsdEta = new TH2F("h_EG_dPhivsdEta","dPhi vs dEta;#Delta#phi(L1EG,L1Tk);#Delta#eta(L1EG,L1Tk)",120,-0.12,0.12,120,-0.12,0.12);
    
  TH2F * h_EG_dPhivsPt = new TH2F("h_EG_dPhivsP","dPhi vs Pt;EGamma pT [Gev];EGamma dPhi",120,0,120,100,-0.1,0.1);
  TH2F * h_EG_dPhivsEta = new TH2F("h_EG_dPhivsEta","dPhi vs Eta;L1EG #eta;#Delta#phi(L1EG,L1Tk)",100,-2.5,2.5,100,-0.1,0.1);
    
  TH2F * h_EG_dEtavsPt = new TH2F("h_EG_dEtavsP","dEta vs Pt;EGamma pT [Gev];EGamma dEta",120,0,120,100,-0.1,0.1);
  TH2F * h_EG_dEtavsEta = new TH2F("h_EG_dEtavsEta","dEta vs Eta;L1EG #eta;L1EG #Delta#eta",100,-2.5,2.5,100,-0.1,0.1);
    
  TH2F * h_EG_dZvsPt = new TH2F("h_EG_dZvsP","dZ vs Pt;L1EG pT [Gev];L1EG #Delta Z",120,0,120,100,-15,15);
  TH2F * h_EG_dZvsEta = new TH2F("h_EG_dZvsEta","dZ vs Eta;L1EG #eta;L1EG #Delta Z",100,-2.5,2.5,100,-15,15);
    
  TH2F * h_EG_dLvsPt = new TH2F("h_EG_dLvsP","dL vs Pt;L1EG pT [Gev];L1EG #Delta L",120,0,120,100,-5,5);
  TH2F * h_EG_dLvsEta = new TH2F("h_EG_dLvsEta","dL vs Eta;L1EG #eta;L1EG #Delta L",100,-2.5,2.5,100,-5,5);
  
  


  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
  
  int nevt = tree->GetEntries();
  std::cout << "number of events = " << nevt << endl;


  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  for (int i=0; i<nevt; i++) {
    int barWidth = 50;
    std::cout << "[";
    float pos = barWidth * (i+1) /nevt;
    for (int k = 0; k < barWidth; ++k) {
      if(k<pos) std::cout << "=";
      else if (k == pos) std::cout << ">";
      else std::cout << " ";
              }
    if(i+1<nevt){
    std::cout << "] " << int((i+1) * 100.0/nevt) << " % Looping Over Event\r";}
    else{
    std::cout << "] " << int((i+1) * 100.0/nevt) << " % Done with Event Loop\r";
      }
    std::cout.flush();
    

    tree->GetEntry(i,0);
    
    auto dr=dR(*gen_pt,*gen_eta,*gen_phi,*trk_pt,*trk_eta,*trk_phi); //return delta R for further study
    auto dphitrack=dPhiTrack(*gen_pt,*gen_eta,*gen_phi,*trk_pt,*trk_eta,*trk_phi,*trk_charge,0.1);// gen particle match to track
    auto elec_dphi=dPhiAtDRMatchPtVMaker(*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,*trk_pt,*trk_eta,*trk_phi,*trk_charge,*trk_z0,*trk_d0,0.08,ifcorr,dphimatching,ifiso);//L1EG match to track, dphi=9999 means no match
    auto elec_deta=dEtaAtDRMatchPtVMaker(*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,*trk_pt,*trk_eta,*trk_phi,*trk_charge,*trk_z0,0.1,ifcorr,dphimatching,ifiso); // return deta for further study
    auto elec_dz=dZAtDRMatchPtVMaker(*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,*trk_pt,*trk_eta,*trk_phi,*trk_charge,*trk_z0,0.1,ifcorr);//return dz for further study
    auto elec_dl=dLAtDRMatchPtVMaker(*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,*trk_pt,*trk_eta,*trk_phi,*trk_charge,*trk_z0,0.1,ifcorr);//return dl for further study
    auto dphiegamma=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,false,false,1); //gen particle match to l1eg, dphi=9999 means no match
    auto dphitkeg=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,false,1);//gen particle match to l1electron, dphi=9999 means no match
    auto dphi_pt=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,true,1);// return all kinds of parameters for further study
    auto deta_pt=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,true,2);
    auto dz_pt=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,true,3);
    auto dl_pt=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,true,4);
    auto et_pt=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,true,5);
    auto eta_pt=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,true,6);
    auto phi_pt=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,true,7);
    auto hgc_pt=dEgamma(*gen_pt,*gen_eta,*gen_phi,*elec_et,*elec_eta,*elec_phi,*elec_looseTkID,*elec_hgc,elec_dphi,elec_deta,elec_dz,elec_dl,0.3,true,true,8);
      
      for (int j=0; j<gen_pt->size();j++){
          GenP_pt.push_back(gen_pt->at(j));
          GenP_eta.push_back(gen_eta->at(j));
          eg_dphitrk.push_back(dphitrack.at(j));// gen matched to l1tk?
          eg_dphieg.push_back(dphiegamma.at(j));// gen matched to l1eg?
          eg_dphitkeg.push_back(dphitkeg.at(j));// gen mathced to l1electron?
          GenP_dR.push_back(dr.at(j));
          
          EG_et.push_back(et_pt.at(j));
          EG_eta.push_back(eta_pt.at(j));
          EG_hgc.push_back(hgc_pt.at(j));
          EG_dPhi.push_back(dphi_pt.at(j));
          EG_dEta.push_back(deta_pt.at(j));
          EG_dZ.push_back(dz_pt.at(j));
          EG_dL.push_back(dl_pt.at(j));
      }
   }
    std::cout << std::endl;

    for (int it=0;it<GenP_pt.size();it++){
        
        int barWidth = 50;
        std::cout << "[";
        float pos = barWidth * (it+1) /GenP_pt.size();
        for (int k = 0; k < barWidth; ++k) {
            if(k<pos) std::cout << "=";
            else if (k == pos) std::cout << ">";
            else std::cout << " ";
        }
        if(it+1<GenP_pt.size()){
            std::cout << "] " << int((it+1) * 100.0/GenP_pt.size()) << " % Looping Over Gen Particles\r";}
        else{
            std::cout << "] " << int((it+1) * 100.0/GenP_pt.size()) << " % Done with Gen Particle Loop\r";
        }
        std::cout.flush();
        
        if (fabs(GenP_eta.at(it)) > 2.4) continue;
        if (isbarrel==0&&fabs(GenP_eta.at(it)) > 1.479) continue;
        if (isbarrel==1&&fabs(GenP_eta.at(it)) < 1.479) continue;
        
        count_all++;// all the gen particles
        h_gen_pt->Fill(GenP_pt.at(it));
        h_gen_eta->Fill(GenP_eta.at(it));
        if (GenP_pt.at(it)<40){
            count_all_L++;
            h_gen_eta_L->Fill(GenP_eta.at(it));
        }
        else{
            count_all_H++;
            h_gen_eta_H->Fill(GenP_eta.at(it));
        }
        if (eg_dphitrk.at(it)<9998){ // eff for l1tk
            count_trk++;
        h_match_trk_pt->Fill(GenP_pt.at(it));
        h_match_trk_eta->Fill(GenP_eta.at(it));
        if (GenP_pt.at(it)<40){
            count_trk_L++;
            h_match_trk_eta_L->Fill(GenP_eta.at(it));
        }
        else{
            count_trk_H++;
            h_match_trk_eta_H->Fill(GenP_eta.at(it));
        }
        }
        
        if (eg_dphieg.at(it)<9998){// eff for l1eg
            count_eg++;
            h_match_eg_pt->Fill(GenP_pt.at(it));
            h_match_eg_eta->Fill(GenP_eta.at(it));
            if (GenP_pt.at(it)<40){
                count_eg_L++;
                h_match_eg_eta_L->Fill(GenP_eta.at(it));
            }
            else{
                count_eg_H++;
                h_match_eg_eta_H->Fill(GenP_eta.at(it));
            }
        }
        
        if (eg_dphitkeg.at(it)<9998){ // eff for l1electron, this is the one I used to make efficiency plots
            count_egtrk++;
            h_match_egtrk_pt->Fill(GenP_pt.at(it));
            h_match_egtrk_eta->Fill(GenP_eta.at(it));
            if (GenP_pt.at(it)<40){
                count_egtrk_L++;
                h_match_egtrk_eta_L->Fill(GenP_eta.at(it));
            }
            else{
                count_egtrk_H++;
                h_match_egtrk_eta_H->Fill(GenP_eta.at(it));
            }
        }
        
        h_res_dR->Fill(GenP_dR.at(it));// dR resoltution
        
        if (GenP_pt.at(it)>=0&&GenP_pt.at(it)<10){
            h_res_dR_1->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=10&&GenP_pt.at(it)<20){
            h_res_dR_2->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=20&&GenP_pt.at(it)<30){
            h_res_dR_3->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=30&&GenP_pt.at(it)<40){
            h_res_dR_4->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=40&&GenP_pt.at(it)<50){
            h_res_dR_5->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=50&&GenP_pt.at(it)<60){
            h_res_dR_6->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=60&&GenP_pt.at(it)<70){
            h_res_dR_7->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=70&&GenP_pt.at(it)<80){
            h_res_dR_8->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=80&&GenP_pt.at(it)<90){
            h_res_dR_9->Fill(GenP_dR.at(it));
        }
        if (GenP_pt.at(it)>=90&&GenP_pt.at(it)<100){
            h_res_dR_10->Fill(GenP_dR.at(it));
        }
        
        
    }
    std::cout << std::endl;
    
    float dphi_matching=9998;
    float deta_matching=9998;
    
    for (int i=0;i<EG_et.size();i++){
        
        int barWidth = 50;
        std::cout << "[";
        float pos = barWidth * (i+1) /EG_et.size();
        for (int k = 0; k < barWidth; ++k) {
            if(k<pos) std::cout << "=";
            else if (k == pos) std::cout << ">";
            else std::cout << " ";
        }
        if(i+1<GenP_pt.size()){
            std::cout << "] " << int((i+1) * 100.0/EG_et.size()) << " % Looping Over L1EGs\r";}
        else{
            std::cout << "] " << int((i+1) * 100.0/EG_et.size()) << " % Done with L1EG Loop\r";
        }
        std::cout.flush();
        
        if (fabs(EG_eta.at(i))>TP_maxEta) continue;
        //if (EG_et.at(i)>20) continue;
        //if (EG_et.at(i)>50) continue;
        
        if (isbarrel==0&&fabs(EG_eta.at(i)) > 1.479) continue;
        if (isbarrel==1&&fabs(EG_eta.at(i)) < 1.479) continue;
        
        if (EG_dPhi.at(i)>=9998) continue;
        
        h_EG_dPhivsdEta->Fill(EG_dPhi.at(i),EG_dEta.at(i));
        
        h_EG_dPhivsPt->Fill(EG_et.at(i),EG_dPhi.at(i));
        h_EG_dPhivsEta->Fill(EG_eta.at(i),EG_dPhi.at(i));
        
        h_EG_dEtavsPt->Fill(EG_et.at(i),EG_dEta.at(i));
        h_EG_dEtavsEta->Fill(EG_eta.at(i),EG_dEta.at(i));
        
        h_EG_dZvsPt->Fill(EG_et.at(i),EG_dZ.at(i));
        h_EG_dZvsEta->Fill(EG_eta.at(i),EG_dZ.at(i));
        h_EG_dZ->Fill(EG_dZ.at(i));
        
        h_EG_dLvsPt->Fill(EG_et.at(i),EG_dL.at(i));
        h_EG_dLvsEta->Fill(EG_eta.at(i),EG_dL.at(i));
        h_EG_dL->Fill(EG_dL.at(i));
        
        h_EG_dPhi->Fill(EG_dPhi.at(i));
        h_EG_dEta->Fill(EG_dEta.at(i));
        h_EG_et->Fill(EG_et.at(i));
        h_EG_eta->Fill(EG_eta.at(i));
        
        if (EG_hgc.at(i)==1){
            h_EG_eta_BB->Fill(EG_eta.at(i));
            h_EG_et_BB->Fill(EG_et.at(i));
            if (fabs(EG_dPhi.at(i))<dphi_matching&&fabs(EG_dEta.at(i))<deta_matching){
            h_matched_EG_eta_BB->Fill(EG_eta.at(i));
            h_matched_EG_et_BB->Fill(EG_et.at(i));
            }
        }
        if (EG_hgc.at(i)){
            h_EG_et_BE->Fill(EG_et.at(i));
        }
    }
    std::cout << std::endl;
    
  // ----------------------------------------------------------------------------------------------------------------
//
//      TF1* fit = new TF1("fit", "gaus", -100.0,100.0);
//      std::vector<double> sigma1;
//      std::vector<double> error1;
//      float bins[]={0,10,20,30,40,50,60,70,80,90,100};
//      int binnum = 10;
//      TH1F *resdR = new TH1F("resdR","dR Resolution",binnum,bins);
//      resdR->GetXaxis()->SetTitle("Gen Particle p_{T} [GeV]");
//      resdR->GetYaxis()->SetTitle("#sigma(#DeltaR)");
//      resdR->SetMinimum(0.0);
//      resdR->SetStats(false);
//
//      h_res_dR_1->Fit(fit,"R");
//      h_res_dR_2->Fit(fit,"R");
//      h_res_dR_3->Fit(fit,"R");
//      h_res_dR_4->Fit(fit,"R");
//      h_res_dR_5->Fit(fit,"R");
//      h_res_dR_6->Fit(fit,"R");
//      h_res_dR_7->Fit(fit,"R");
//      h_res_dR_8->Fit(fit,"R");
//      h_res_dR_9->Fit(fit,"R");
//      h_res_dR_10->Fit(fit,"R");
//
//      //h_res_dR->Fit(fit,"R");
//
//      sigma1.push_back(h_res_dR_1->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_2->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_3->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_4->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_5->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_6->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_7->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_8->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_9->GetFunction("fit")->GetParameter(2));
//      sigma1.push_back(h_res_dR_10->GetFunction("fit")->GetParameter(2));
//
//      error1.push_back(h_res_dR_1->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_2->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_3->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_4->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_5->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_6->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_7->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_8->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_9->GetFunction("fit")->GetParError(2));
//      error1.push_back(h_res_dR_10->GetFunction("fit")->GetParError(2));
//
//      for(int n=0;n<10;n++) {
//            resdR->SetBinContent(n+1,sigma1[n]);
//            resdR->SetBinError(n+1,error1[n]);
//            //cout<<sigma1[n]<<endl;
//        }

bool barrel=isbarrel==0;
bool endcap=isbarrel==1;
  
  TFile* fout;
  fout = new TFile("output_"+type+treeName+string(ifcorr?"_withcorr":"_withoutcorr")+string(dphimatching?"_withdphi":"_withoutdphi")+string(ifiso?"_withiso":"_withoutiso")+string(barrel?"_BB":endcap?"_BE":"_all")+".root","recreate");

    
  char ctxt[500];
  TCanvas c;
  TString DIR = "TrkPlots/";

  // ----------------------------------------------------------------------------------------------------------------
  // efficiency plots  
  // ----------------------------------------------------------------------------------------------------------------

//  // rebin pt/phi plots
  h_gen_pt->Rebin(4);
  h_match_eg_pt->Rebin(4);
  h_match_trk_pt->Rebin(4);
  h_match_egtrk_pt->Rebin(4);

  h_gen_eta->Rebin(2);
  h_match_eg_eta->Rebin(2);
  h_match_trk_eta->Rebin(2);
  h_match_egtrk_eta->Rebin(2);
  h_gen_eta_L->Rebin(2);
  h_match_eg_eta_L->Rebin(2);
  h_match_egtrk_eta_L->Rebin(2);
  h_match_trk_eta_L->Rebin(2);
  h_gen_eta_H->Rebin(2);
  h_match_eg_eta_H->Rebin(2);
  h_match_trk_eta_H->Rebin(2);
  h_match_egtrk_eta_H->Rebin(2);

  // calculate the effeciency
  h_match_eg_pt->Sumw2();
  h_gen_pt->Sumw2();
  TH1F* h_eff_pt_eg = (TH1F*) h_match_eg_pt->Clone();
  h_eff_pt_eg->SetName("eff_pt_l1eg");
  h_eff_pt_eg->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_eg->Divide(h_match_eg_pt, h_gen_pt, 1.0, 1.0, "B");

  h_match_eg_eta->Sumw2();
  h_gen_eta->Sumw2();
  TH1F* h_eff_eta_eg = (TH1F*) h_match_eg_eta->Clone();
  h_eff_eta_eg->SetName("eff_eta_l1eg");
  h_eff_eta_eg->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_eg->Divide(h_match_eg_eta, h_gen_eta, 1.0, 1.0, "B");

  h_match_eg_eta_L->Sumw2();
  h_gen_eta_L->Sumw2();
  TH1F* h_eff_eta_eg_L = (TH1F*) h_match_eg_eta_L->Clone();
  h_eff_eta_eg_L->SetName("eff_eta_l1eg_L");
  h_eff_eta_eg_L->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_eg_L->Divide(h_match_eg_eta_L, h_gen_eta_L, 1.0, 1.0, "B");

  h_match_eg_eta_H->Sumw2();
  h_gen_eta_H->Sumw2();
  TH1F* h_eff_eta_eg_H = (TH1F*) h_match_eg_eta_H->Clone();
  h_eff_eta_eg_H->SetName("eff_eta_l1eg_H");
  h_eff_eta_eg_H->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_eg_H->Divide(h_match_eg_eta_H, h_gen_eta_H, 1.0, 1.0, "B");



h_match_trk_pt->Sumw2();
TH1F* h_eff_pt_trk = (TH1F*) h_match_trk_pt->Clone();
h_eff_pt_trk->SetName("eff_pt_l1tk");
h_eff_pt_trk->GetYaxis()->SetTitle("Efficiency");
h_eff_pt_trk->Divide(h_match_trk_pt, h_gen_pt, 1.0, 1.0, "B");

h_match_trk_eta->Sumw2();
TH1F* h_eff_eta_trk = (TH1F*) h_match_trk_eta->Clone();
h_eff_eta_trk->SetName("eff_eta_l1tk");
h_eff_eta_trk->GetYaxis()->SetTitle("Efficiency");
h_eff_eta_trk->Divide(h_match_trk_eta, h_gen_eta, 1.0, 1.0, "B");

h_match_trk_eta_L->Sumw2();
TH1F* h_eff_eta_trk_L = (TH1F*) h_match_trk_eta_L->Clone();
h_eff_eta_trk_L->SetName("eff_eta_l1tk_L");
h_eff_eta_trk_L->GetYaxis()->SetTitle("Efficiency");
h_eff_eta_trk_L->Divide(h_match_trk_eta_L, h_gen_eta_L, 1.0, 1.0, "B");

h_match_trk_eta_H->Sumw2();
TH1F* h_eff_eta_trk_H = (TH1F*) h_match_trk_eta_H->Clone();
h_eff_eta_trk_H->SetName("eff_eta_l1tk_H");
h_eff_eta_trk_H->GetYaxis()->SetTitle("Efficiency");
h_eff_eta_trk_H->Divide(h_match_trk_eta_H, h_gen_eta_H, 1.0, 1.0, "B");





h_match_egtrk_pt->Sumw2();
TH1F* h_eff_pt_egtrk = (TH1F*) h_match_egtrk_pt->Clone();
h_eff_pt_egtrk->SetName("eff_pt_l1ele");
h_eff_pt_egtrk->GetYaxis()->SetTitle("Efficiency");
h_eff_pt_egtrk->Divide(h_match_egtrk_pt, h_gen_pt, 1.0, 1.0, "B");

h_match_egtrk_eta->Sumw2();
TH1F* h_eff_eta_egtrk = (TH1F*) h_match_egtrk_eta->Clone();
h_eff_eta_egtrk->SetName("eff_eta_l1ele");
h_eff_eta_egtrk->GetYaxis()->SetTitle("Efficiency");
h_eff_eta_egtrk->Divide(h_match_egtrk_eta, h_gen_eta, 1.0, 1.0, "B");

h_match_egtrk_eta_L->Sumw2();
TH1F* h_eff_eta_egtrk_L = (TH1F*) h_match_egtrk_eta_L->Clone();
h_eff_eta_egtrk_L->SetName("eff_eta_l1ele_L");
h_eff_eta_egtrk_L->GetYaxis()->SetTitle("Efficiency");
h_eff_eta_egtrk_L->Divide(h_match_egtrk_eta_L, h_gen_eta_L, 1.0, 1.0, "B");

h_match_egtrk_eta_H->Sumw2();
TH1F* h_eff_eta_egtrk_H = (TH1F*) h_match_egtrk_eta_H->Clone();
h_eff_eta_egtrk_H->SetName("eff_eta_l1ele_H");
h_eff_eta_egtrk_H->GetYaxis()->SetTitle("Efficiency");
h_eff_eta_egtrk_H->Divide(h_match_egtrk_eta_H, h_gen_eta_H, 1.0, 1.0, "B");



  // set the axis range
  h_eff_pt_trk  ->SetAxisRange(0,1.1,"Y");
  h_eff_eta_trk ->SetAxisRange(0,1.1,"Y");
  h_eff_eta_trk_L ->SetAxisRange(0,1.1,"Y");
  h_eff_eta_trk_H ->SetAxisRange(0,1.1,"Y");

h_eff_pt_eg  ->SetAxisRange(0,1.1,"Y");
h_eff_eta_eg ->SetAxisRange(0,1.1,"Y");
h_eff_eta_eg_L ->SetAxisRange(0,1.1,"Y");
h_eff_eta_eg_H ->SetAxisRange(0,1.1,"Y");

h_eff_pt_egtrk  ->SetAxisRange(0,1.1,"Y");
h_eff_eta_egtrk ->SetAxisRange(0,1.1,"Y");
h_eff_eta_egtrk_L ->SetAxisRange(0,1.1,"Y");
h_eff_eta_egtrk_H ->SetAxisRange(0,1.1,"Y");

  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetLogz();

   //draw and save plots
  h_eff_pt_trk->Draw();
  h_eff_pt_trk->Write();
  c.SaveAs(DIR+type+"_eff_pt_l1tk.pdf");

  h_eff_eta_trk->Draw();
  h_eff_eta_trk->Write();
  c.SaveAs(DIR+type+"_eff_eta_l1tk.pdf");


  h_eff_eta_trk_L->Draw();
  h_eff_eta_trk_L->Write();
  sprintf(ctxt,"p_{T} < 40 GeV");
  mySmallText(0.45,0.5,1,ctxt);
  c.SaveAs(DIR+type+"_eff_eta_l1tk_L.pdf");

  h_eff_eta_trk_H->Draw();
  h_eff_eta_trk_H->Write();
  sprintf(ctxt,"p_{T} > 40 GeV");
  mySmallText(0.45,0.5,1,ctxt);
  c.SaveAs(DIR+type+"_eff_eta_l1tk_H.pdf");


h_eff_pt_eg->Draw();
h_eff_pt_eg->Write();
c.SaveAs(DIR+type+"_eff_pt_l1eg.pdf");

h_eff_eta_eg->Draw();
h_eff_eta_eg->Write();
c.SaveAs(DIR+type+"_eff_eta_l1eg.pdf");


h_eff_eta_eg_L->Draw();
h_eff_eta_eg_L->Write();
sprintf(ctxt,"p_{T} < 40 GeV");
mySmallText(0.45,0.5,1,ctxt);
c.SaveAs(DIR+type+"_eff_eta_l1eg_L.pdf");

h_eff_eta_eg_H->Draw();
h_eff_eta_eg_H->Write();
sprintf(ctxt,"p_{T} > 40 GeV");
mySmallText(0.45,0.5,1,ctxt);
c.SaveAs(DIR+type+"_eff_eta_l1eg_H.pdf");




h_eff_pt_egtrk->Draw();
h_eff_pt_egtrk->Write();
c.SaveAs(DIR+type+"_eff_pt_l1ele.pdf");

h_eff_eta_egtrk->Draw();
h_eff_eta_egtrk->Write();
c.SaveAs(DIR+type+"_eff_eta_l1ele.pdf");


h_eff_eta_egtrk_L->Draw();
h_eff_eta_egtrk_L->Write();
sprintf(ctxt,"p_{T} < 40 GeV");
mySmallText(0.45,0.5,1,ctxt);
c.SaveAs(DIR+type+"_eff_eta_l1ele_L.pdf");

h_eff_eta_egtrk_H->Draw();
h_eff_eta_egtrk_H->Write();
sprintf(ctxt,"p_{T} > 40 GeV");
mySmallText(0.45,0.5,1,ctxt);
c.SaveAs(DIR+type+"_eff_eta_l1ele_H.pdf");
    
cout<<"number of gen particles ="<<GenP_pt.size()<<endl;
cout<<"number of EGamma ="<<EG_et.size()<<endl;
    
  float k = (float)count_trk;
  float N = (float)count_all;
  if (fabs(N)>0) cout << endl << "efficiency for L1Tk = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)count_trk_L;
  N = (float)count_all_L;
  if (fabs(N)>0) cout << "efficiency for L1Tk with pT < 40 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)count_trk_H;
  N = (float)count_all_H;
  if (fabs(N)>0) cout << "efficiency for L1Tk with pT > 40 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;

  k = (float)count_eg;
  N = (float)count_all;
  if (fabs(N)>0) cout << endl << "efficiency for L1EG = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)count_eg_L;
  N = (float)count_all_L;
  if (fabs(N)>0) cout << "efficiency for L1EG with pT < 40 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)count_eg_H;
  N = (float)count_all_H;
  if (fabs(N)>0) cout << "efficiency for L1EG with pT > 40 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
    
  k = (float)count_egtrk;
  N = (float)count_all;
  if (fabs(N)>0) cout << endl << "efficiency for L1Eletron = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)count_egtrk_L;
  N = (float)count_all_L;
  if (fabs(N)>0) cout << "efficiency for L1Eletron with pT < 40 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)count_egtrk_H;
  N = (float)count_all_H;
  if (fabs(N)>0) cout << "efficiency for L1Eletron with pT > 40 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  
//  resdR->Draw();
//  resdR->Write();
//  c.SaveAs(DIR+type+"_res_dR.pdf");
//
//  h_res_dR->Draw();
//  h_res_dR->Write();
//  c.SaveAs(DIR+type+"_dis_dR.pdf");
    
    
    
  h_EG_et_BE->Draw();
  h_EG_et_BE->Write();
  c.SaveAs(DIR+type+"_dis_et_BE.pdf");
    
  h_EG_et_BB->Draw();
  h_EG_et_BB->Write();
  c.SaveAs(DIR+type+"_dis_et_BB.pdf");
    
  h_EG_dPhi->Draw();
  h_EG_dPhi->Write();
  c.SaveAs(DIR+type+"_dis_dPhi.pdf");
    
  h_EG_dEta->Draw();
  h_EG_dEta->Write();
  c.SaveAs(DIR+type+"_dis_dEta.pdf");
    
  h_EG_dZ->Draw();
  h_EG_dZ->Write();
  c.SaveAs(DIR+type+"_dis_dZ.pdf");
    
  h_EG_dL->Draw();
  h_EG_dL->Write();
  c.SaveAs(DIR+type+"_dis_dL.pdf");
    

    
  h_EG_dPhivsdEta->Draw("COLZ");
  h_EG_dPhivsdEta->Write();
  c.SaveAs(DIR+type+"_dis_dPhivsdEta.pdf");
    

  h_EG_dPhivsPt->Draw("COLZ");
  h_EG_dPhivsPt->Write();
  c.SaveAs(DIR+type+"_dis_dPhivsPt.pdf");
    
  h_EG_dPhivsEta->Draw("COLZ");
  h_EG_dPhivsEta->Write();
  c.SaveAs(DIR+type+"_dis_dPhivsEta.pdf");
    
  h_EG_dEtavsPt->Draw("COLZ");
  h_EG_dEtavsPt->Write();
  c.SaveAs(DIR+type+"_dis_dEtavsPt.pdf");
    
  h_EG_dEtavsEta->Draw("COLZ");
  h_EG_dPhivsEta->Write();
  c.SaveAs(DIR+type+"_dis_dEtavsEta.pdf");
  
  h_EG_dZvsPt->Draw("COLZ");
  h_EG_dZvsPt->Write();
  h_EG_dZvsPt->SetAxisRange(-15,15,"Y");
  c.SaveAs(DIR+type+"_dis_dZvsPt.pdf");

  h_EG_dZvsEta->Draw("COLZ");
  h_EG_dZvsEta->Write();
  h_EG_dZvsEta->SetAxisRange(-15,15,"Y");
  c.SaveAs(DIR+type+"_dis_dZvsEta.pdf");
    
  h_EG_dLvsPt->Draw("COLZ");
  h_EG_dLvsPt->Write();
  h_EG_dLvsPt->SetAxisRange(-5,5,"Y");
  c.SaveAs(DIR+type+"_dis_dLvsPt.pdf");
    
  h_EG_dLvsEta->Draw("COLZ");
  h_EG_dLvsEta->Write();
  h_EG_dLvsEta->SetAxisRange(-5,5,"Y");
  c.SaveAs(DIR+type+"_dis_dLvsEta.pdf");
    
    

  h_EG_et_BB->Rebin(4);
  h_matched_EG_et_BB->Rebin(4);
    
  h_matched_EG_et_BB->Sumw2();
  h_EG_et_BB->Sumw2();
  TH1F* h_eff_pt_matched_eg_BB = (TH1F*) h_matched_EG_et_BB->Clone();
  h_eff_pt_matched_eg_BB->SetName("eff_pt_matched_eg_BB");
  h_eff_pt_matched_eg_BB->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_matched_eg_BB->Divide(h_eff_pt_matched_eg_BB, h_EG_et_BB, 1.0, 1.0, "B");
    
  h_eff_pt_matched_eg_BB->Draw();
  h_eff_pt_matched_eg_BB->Write();
  h_eff_pt_matched_eg_BB  ->SetAxisRange(0,1.1,"Y");
  c.SaveAs(DIR+type+"_eff_pt_matched_eg_BB.pdf");
    
    
  h_EG_eta_BB->Rebin(2);
  h_matched_EG_eta_BB->Rebin(2);
    
  h_matched_EG_eta_BB->Sumw2();
  h_EG_eta_BB->Sumw2();
  TH1F* h_eff_eta_matched_eg_BB = (TH1F*) h_matched_EG_eta_BB->Clone();
  h_eff_eta_matched_eg_BB->SetName("eff_eta_matched_eg_BB");
  h_eff_eta_matched_eg_BB->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_matched_eg_BB->Divide(h_eff_eta_matched_eg_BB, h_EG_eta_BB, 1.0, 1.0, "B");
    
  h_eff_eta_matched_eg_BB->Draw();
  h_eff_eta_matched_eg_BB->Write();
  h_eff_eta_matched_eg_BB  ->SetAxisRange(0,1.1,"Y");
  c.SaveAs(DIR+type+"_eff_eta_matched_eg_BB.pdf");
    
  
    
    
    
  fout->Close();

//    gROOT->Reset("a");
//    TControlBar bar("vertical");
//    bar.AddButton("Help to run demos",".x demoshelp.C",
//                  "Explains how to run the demos");
//    bar.AddButton("framework",        ".x framework.C",
//                  "An Example of Object Oriented User Interface");
//    bar.AddButton("hsimple",          ".x hsimple.C",
//                  "An Example Creating Histograms/Ntuples on File");
//    bar.AddButton("hsum",             ".x hsum.C",
//                  "Filling histograms and some graphics options");
//    bar.AddButton("canvas",           ".x canvas.C",
//                  "Canvas and Pad Management");
//    bar.AddButton("formula1",         ".x formula1.C",
//                  "Simple Formula and Functions");
//    bar.AddButton("fillrandom",       ".x fillrandom.C",
//                  "Histograms with Random Numbers from a Function");
//    bar.AddButton("fit1",             ".x fit1.C",
//                  "A Simple Fitting Example");
//    bar.AddButton("h1draw",           ".x h1draw.C",
//                  "Drawing Options for 1D Histograms");
//    bar.AddButton("graph",            ".x graph.C",
//                  "Examples of a simple graph");
//    bar.AddButton("tornado",          ".x tornado.C",
//                  "Examples of 3-D PolyMarkers");
//    bar.AddButton("shapes",           ".x shapes.C",
//                  "The Geometry Shapes");
//    bar.AddButton("atlasna49",        ".x atlasna49.C",
//                  "Creating and Viewing Geometries");
//    bar.AddButton("file_layout",      ".x file.C",
//                  "The ROOT file format");
//    bar.AddButton("tree_layout",      ".x tree.C",
//                  "The Tree Data Structure");
//    bar.AddButton("ntuple1",          ".x ntuple1.C",
//                  "Ntuples and Selections");
//    bar.AddButton("run benchmarks",   ".x benchmarks.C",
//                  "Runs all the ROOT benchmarks");
//    bar.AddButton("rootmarks",        ".x rootmarks.C",
//                  "Prints an estimated ROOTMARKS for your machine");
//    bar.AddButton("edit_hsimple",     ".!ved hsimple.C &",
//                  "Invokes the text editor on file hsimple.C");
//    bar.AddButton("Close Bar",        "gROOT.Reset(\"a\")",
//                  "Close ControlBar");
//    bar.Show();
//    gROOT->SaveContext();

}


void SetPlotStyle() {

  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(kBird);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}

void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
    Double_t tsize=0.044;
    TLatex l;
    l.SetTextSize(tsize);
    l.SetNDC();
    l.SetTextColor(color);
    l.DrawLatex(x,y,text);
}

// This is one of the most important functions in this program, it takes l1eg objects and try to match it with l1 tracks, corrections/isolation/matching criteria can be specified. This function returns a vector of dphi values. dphi = 9999 means l1eg has no matched l1tk.
vector<float> dPhiAtDRMatchPtVMaker ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<int> &mainlooseTkID, vector<int> &mainhgc, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, vector<float> &loopz0, vector<float> &loopd0, float dr,bool phietacorr, bool dphimatch, bool iso){
    vector<float> perobjectmindphi;
    TLorentzVector r, g, g1, g2;
    //
    for (std::size_t i = 0; i != maineta.size(); ++i) {
        float dphimin = 9999.;
        float ptmatch   = 0.;
        float detacut = 0;
        float dphicut = 0;
        float dphideta = 9999. ;
        float max_iso = 0;
        if (mainhgc.at(i)==1){
            max_iso=0.4;// isolation cut in hgcal
            dphicut=0.07;// dimensions of ellipse
            detacut=0.0075;
        }
        else{
            max_iso= 0.12;
            dphicut=0.03;
            detacut=0.025;
            if (fabs(maineta.at(i))>0.9){
                detacut=0.015;
            }
            if (mainpt.at(i)<50){
                dphicut=0.045;
            }
        }
        
        r.SetPtEtaPhiM( mainpt.at(i), maineta.at(i), mainphi.at(i), 0. );//assume this is cluster (doesnt have to be)
        //
        if (mainlooseTkID.at(i)==1){
            vector<float> candpt;// vectors to store electron candidate
            vector<float> candeta;
            vector<float> candphi;
            vector<float> candz0;
            vector<float> candd0;
            //
            for (std::size_t j = 0; j != loopeta.size(); ++j) { // assume this is loop over tracks (doesnt have to be)
                if (looppt.at(j)<5) continue; // add pt cut on tracks
                if (std::isnan(loopeta.at(j))) continue;
//                if (mainhgc.at(i)==1) {
//                    float theta = 2*TMath::ATan(TMath::Exp(-(loopeta.at(j))));
//                    if(fabs(2.47*0.3*fabs(TMath::Tan(theta)*320. /130. )*loopcharge.at(j)/mainpt.at(i))>1) continue;
//                }
                auto loopphicorr=phiCorr(mainpt.at(i), loopeta.at(j), loopphi.at(j),  loopcharge.at(j),  mainhgc.at(i), 2.47*0.3); // make phi corrections
                auto loopetacorr=etaCorr(loopeta.at(j), loopz0.at(j), mainhgc.at(i));// make eta corrections
                g.SetPtEtaPhiM( looppt.at(j), loopeta.at(j), loopphi.at(j), 0. );// assume no correcions
                if(phietacorr){
                    g.SetPtEtaPhiM( looppt.at(j), loopetacorr, loopphicorr, 0. );// turn on corrections
                }
                if(!iso){
                    //&&fabs(r.DeltaPhi(g))<0.07
                if(!dphimatch&& fabs(r.DeltaR(g)) < dr && ptmatch < g.Pt() &&fabs(r.DeltaPhi(g))<0.07){// standard matching
                    ptmatch = g.Pt();
                    dphimin = r.DeltaPhi(g);
                }
                dphideta=r.DeltaPhi(g)*r.DeltaPhi(g)/(dphicut*dphicut)+(r.Eta()-g.Eta())*(r.Eta()-g.Eta())/(detacut*detacut);// elliptical matching
                if(dphimatch&& dphideta<1 && ptmatch < g.Pt()){
                    ptmatch = g.Pt();
                    dphimin = r.DeltaPhi(g);
                }
                }
                if(iso){
                    if(!dphimatch&& fabs(r.DeltaR(g)) < dr &&fabs(r.DeltaPhi(g))<0.07){
                        candpt.push_back(looppt.at(j));// only consider electron candidate
                        candeta.push_back(loopeta.at(j));
                        candphi.push_back(loopphi.at(j));
                        candz0.push_back(loopz0.at(j));
                        candd0.push_back(loopd0.at(j));

                    }
                    dphideta=r.DeltaPhi(g)*r.DeltaPhi(g)/(dphicut*dphicut)+(r.Eta()-g.Eta())*(r.Eta()-g.Eta())/(detacut*detacut);
                    if(dphimatch&& dphideta<1){
                        candpt.push_back(looppt.at(j));
                        candeta.push_back(loopeta.at(j));
                        candphi.push_back(loopphi.at(j));
                        candz0.push_back(loopz0.at(j));
                        candd0.push_back(loopd0.at(j));
                        
                    }
                    
                }
                
            }
            if(iso){// add isolation cut
                for (std::size_t k = 0; k != candeta.size(); ++k){
                    g1.SetPtEtaPhiM( candpt.at(k), candeta.at(k), candphi.at(k), 0. );
                    auto rel_iso=rel_isolation(candpt.at(k),candeta.at(k),candphi.at(k),candz0.at(k),candd0.at(k),looppt,loopeta,loopphi,loopz0,loopd0);
                    if (rel_iso<max_iso&&ptmatch<g1.Pt()){
                        ptmatch = g1.Pt();
                        dphimin=r.DeltaPhi(g1);
                    }
                }
            }
        }
        perobjectmindphi.emplace_back( dphimin );
    }
    //
    return perobjectmindphi;
}
// This should be straightfoward, it returns a vector of deta values.
vector<float> dEtaAtDRMatchPtVMaker ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<int> &mainlooseTkID, vector<int> &mainhgc, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, vector<float> &loopz0, float dr, bool phietacorr, bool dphimatch, bool iso){
    vector<float> perobjectmindeta;
    TLorentzVector r, g, g1, g2;
    //
    for (std::size_t i = 0; i != maineta.size(); ++i) {
        float detamin = 9999.;
        float ptmatch   = 0.;
        float detacut = 0;
        float dphicut = 0;
        float dphideta = 9999. ;
        if (mainhgc.at(i)==1){
            dphicut=0.07;
            detacut=0.007;
        }
        else{
            dphicut=0.03;
            detacut=0.015;
            if (fabs(maineta.at(i))>0.9){
                detacut=0.015;
            }
            if (mainpt.at(i)<50){
                dphicut=0.045;
            }
        }
        
        r.SetPtEtaPhiM( mainpt.at(i), maineta.at(i), mainphi.at(i), 0. );//assume this is cluster (doesnt have to be)
        
        if (mainlooseTkID.at(i)==1){
            //
            for (std::size_t j = 0; j != loopeta.size(); ++j) { // assume this is loop over tracks (doesnt have to be)
                if (looppt.at(j)<5) continue;
                if (std::isnan(loopeta.at(j))) continue;
//                if (mainhgc.at(i)==1) {
//                    float theta = 2*TMath::ATan(TMath::Exp(-(loopeta.at(j))));
//                    if(fabs(2.47*0.3*fabs(TMath::Tan(theta)*320. /130. )*loopcharge.at(j)/mainpt.at(i))>1) continue;
//                }
                auto loopphicorr=phiCorr(mainpt.at(i), loopeta.at(j), loopphi.at(j),  loopcharge.at(j),  mainhgc.at(i), 2.47*0.3);
                auto loopetacorr=etaCorr(loopeta.at(j), loopz0.at(j), mainhgc.at(i));
                g.SetPtEtaPhiM( looppt.at(j), loopeta.at(j), loopphi.at(j), 0. );
                if(phietacorr){
                    g.SetPtEtaPhiM( looppt.at(j), loopetacorr, loopphicorr, 0. );
                }
                if(!iso&&!dphimatch&& fabs(r.DeltaR(g)) < dr && ptmatch < g.Pt() ){
                    ptmatch = g.Pt();
                    detamin = r.Eta()-g.Eta();// use cluster (aka main) pt to correct for phi!
                }
                dphideta=r.DeltaPhi(g)*r.DeltaPhi(g)/(dphicut*dphicut)+(r.Eta()-g.Eta())*(r.Eta()-g.Eta())/(detacut*detacut);
                if(!iso&&dphimatch&& dphideta<1 && ptmatch < g.Pt() ){
                    ptmatch = g.Pt();
                    detamin = r.Eta()-g.Eta();
                }
            }
        }
        perobjectmindeta.emplace_back( detamin );
    }
    //
    return perobjectmindeta;
}

vector<float> dZAtDRMatchPtVMaker ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<int> &mainlooseTkID, vector<int> &mainhgc, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, vector<float> &loopz0, float dr, bool phietacorr, float correction=0, float trackZ0multiplier=1.0 ){
    vector<float> perobjectmindz;
    TLorentzVector r, g, g1, g2;
    //
    for (std::size_t i = 0; i != maineta.size(); ++i) {
        float dzmin   = 9999.;
        float ptmatch = 0.;
        r.SetPtEtaPhiM( mainpt.at(i), maineta.at(i), mainphi.at(i), 0. );//assume this is cluster (doesnt have to be)
        if (mainlooseTkID.at(i)==1){
        for (std::size_t j = 0; j != loopeta.size(); ++j) { // assume this is loop over tracks (doesnt have to be)
            if (std::isnan(loopeta.at(j))) continue;
            if (mainhgc.at(i)==1) continue;
            auto loopphicorr=phiCorr(mainpt.at(i), loopeta.at(j), loopphi.at(j), loopcharge.at(j), mainhgc.at(i), 2.47*0.3);
            auto loopetacorr=etaCorr(loopeta.at(j), loopz0.at(j),mainhgc.at(i));
            g.SetPtEtaPhiM( looppt.at(j), loopeta.at(j), loopphi.at(j), 0. );
            if(phietacorr){
                g.SetPtEtaPhiM( looppt.at(j), loopetacorr, loopphicorr, 0. );
            }
            if( fabs(r.DeltaR(g)) < dr && ptmatch < g.Pt() ){
                ptmatch = g.Pt();
                float thetamain = 2*TMath::ATan(TMath::Exp(-(maineta.at(i))));
                float zmain     = 130.*(1.+correction)/TMath::Tan(thetamain);
                float thetaloop = 2*TMath::ATan(TMath::Exp(-(loopeta.at(j))));
                float zloop     = 130.*(1.+correction)/TMath::Tan(thetaloop);
                if(phietacorr){
                    zloop = zloop+loopz0.at(j);
                }
                //float zloop     = trackZ0multiplier*loopz0.at(j) + 130.*(1.+correction)/TMath::Tan(thetaloop);
                if (fabs(zmain)>320||fabs(zloop)>320) continue;
                dzmin = zmain - zloop;
            }
        }
        // ECAL crystal geometry correction in PLUS and MINUS sides
//        if( maineta.at(i)<0.0) dzmin = dzmin + 0.72;
//        if( maineta.at(i)>0.0) dzmin = dzmin - 0.72;
        }
        perobjectmindz.emplace_back( dzmin );
    }
    //
    return perobjectmindz;
}

vector<float> dLAtDRMatchPtVMaker ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<int> &mainlooseTkID, vector<int> &mainhgc, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, vector<float> &loopz0, float dr, bool phietacorr, float correction=0, float trackZ0multiplier=1.0 ){
    vector<float> perobjectmindz;
    TLorentzVector r, g, g1, g2;
    //
    for (std::size_t i = 0; i != maineta.size(); ++i) {
        float dzmin   = 9999.;
        float ptmatch = 0.;
        r.SetPtEtaPhiM( mainpt.at(i), maineta.at(i), mainphi.at(i), 0. );//assume this is cluster (doesnt have to be)
        if (mainlooseTkID.at(i)==1){
            for (std::size_t j = 0; j != loopeta.size(); ++j) { // assume this is loop over tracks (doesnt have to be)
                if (std::isnan(loopeta.at(j))) continue;
                if (mainhgc.at(i)==0) continue;
                auto loopphicorr=phiCorr(mainpt.at(i), loopeta.at(j), loopphi.at(j), loopcharge.at(j), mainhgc.at(i), 2.47*0.3);
                auto loopetacorr=etaCorr(loopeta.at(j), loopz0.at(j),mainhgc.at(i));
                g.SetPtEtaPhiM( looppt.at(j), loopeta.at(j), loopphi.at(j), 0. );
                if(phietacorr){
                    g.SetPtEtaPhiM( looppt.at(j), loopetacorr, loopphicorr, 0. );
                }
                if( fabs(r.DeltaR(g)) < dr && ptmatch < g.Pt() ){
                    ptmatch = g.Pt();
                    float thetamain = 2*TMath::ATan(TMath::Exp(-(maineta.at(i))));
                    float lmain     = 320.*TMath::Tan(thetamain);
                    float thetaloop = 2*TMath::ATan(TMath::Exp(-(loopeta.at(j))));
                    float lloop     = 320.*TMath::Tan(thetaloop);
                    if(phietacorr){
                        lloop = fabs(lloop/fabs(lloop)*320.-loopz0.at(j))*TMath::Tan(thetaloop);
                    }
                    //float zloop     = trackZ0multiplier*loopz0.at(j) + 130.*(1.+correction)/TMath::Tan(thetaloop);
                    if (fabs(lmain)>130||fabs(lloop)>130) continue;
                    dzmin = lmain - lloop;
                }
            }
            // ECAL crystal geometry correction in PLUS and MINUS sides
            //        if( maineta.at(i)<0.0) dzmin = dzmin + 0.72;
            //        if( maineta.at(i)>0.0) dzmin = dzmin - 0.72;
        }
        perobjectmindz.emplace_back( dzmin );
    }
    //
    return perobjectmindz;
}
// This function is also import, it returns corrected phi angle.
float phiCorr (float &pt,float &eta,float &phi,int &charge,int &hgc,float phicorrection){
        float pi = (float)(TMath::Pi());
        float effectiveradiusratio = 1.;
        float theta            = 2*TMath::ATan(TMath::Exp(-(eta)));
        float effectiveradius      = fabs(TMath::Tan(theta)*320.);
        if(hgc==1){ effectiveradiusratio = effectiveradius / 130.; }
        float correctedphi = phi;
        float dsin=phicorrection*effectiveradiusratio*charge/pt;
        //if(std::isnan(theta)){cout<<"Warning:NAN"<<eta<<endl;}
        //if (dsin>1){cout<<dsin << pt<<endl;dsin=1;}
        //if (dsin<-1){dsin=-1;}
    
        correctedphi = phi - TMath::ASin(dsin);//Jack's correction
        //correctedphi = phi - 0.85/phicorrection*dsin;//Halil's correction
    
    if(correctedphi>pi){correctedphi = correctedphi - 2*pi;}
    if(correctedphi<-1*pi){correctedphi = correctedphi + 2*pi;}
    if(hgc==0&&eta>0){ //manually shift the distribution
        correctedphi=correctedphi-0.0025;
    }
    if(hgc==0&&eta<0){
        correctedphi=correctedphi+0.0025;
    }
    if(hgc==1&&eta>0){
        correctedphi=correctedphi+0.0015;
    }
    if(hgc==1&&eta<0){
        correctedphi=correctedphi+0.0015;
    }
    return correctedphi;
}

// corrections for eta
float etaCorr (float &eta,float &z0,int &hgc){
    float thetacorr;
    float perp=130;
    float pi = (float)(TMath::Pi());
    float theta=2*TMath::ATan(TMath::Exp(-(eta)));
    float z=perp/TMath::Tan(theta);
    float zcorr=z+z0;
    if (hgc==1&&fabs(eta)>1.4){
        zcorr=320*z/fabs(z);
        perp=fabs((zcorr-z0)*TMath::Tan(theta));
    }
    if (zcorr>0){
        thetacorr=TMath::ATan(perp/zcorr);
        }
    else {
        thetacorr=pi-TMath::ATan(perp/zcorr);
    }
    if( thetacorr >pi ) {thetacorr = 2*pi - thetacorr;}
    if( thetacorr <0 && thetacorr > -1 * pi) {thetacorr = -1 * thetacorr;}
    if( thetacorr < -1* pi ) {thetacorr = 2*pi + thetacorr;}
    float correctedeta=-1.0 * log(TMath::Tan(thetacorr/2.0));
    if(hgc==0&&eta>0){
        correctedeta=correctedeta-0.0025;
    }
    if(hgc==0&&eta<0){
        correctedeta=correctedeta+0.0025;
    }
    if(hgc==1&&eta>0){
        correctedeta=correctedeta+0.0015;
    }
    if(hgc==1&&eta<0){
        correctedeta=correctedeta-0.0015;
    }
    return correctedeta;
}

// This is simple function: match gen particle with l1tk, no correcions
vector<float> dPhiTrack ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &loopcharge, float dr){
    vector<float> perobjectmindphi;
    TLorentzVector r, g, g1, g2;
    //
    for (std::size_t i = 0; i != maineta.size(); ++i) {
        float dphimin = 9999.;
        float ptmatch   = 0.;
        if (maineta.at(i)>4) continue;
        r.SetPtEtaPhiM( mainpt.at(i), maineta.at(i), mainphi.at(i), 0. );//assume this is cluster (doesnt have to be)

        for (std::size_t j = 0; j != loopeta.size(); ++j) { // assume this is loop over tracks (doesnt have to be)
            //if (looppt.at(j)<40||fabs(loopeta.at(j))>1.5){dr=0.8;}
//            float factor=5;
//            if (mainpt.at(i)>0&&mainpt.at(i)<=10)
//            {dr=0.0287*factor;}
//            else if (mainpt.at(i)>10&&mainpt.at(i)<=20)
//            {dr=0.0250*factor;}
//            else if (mainpt.at(i)>20&&mainpt.at(i)<=30)
//            {dr=0.0244*factor;}
//            else if (mainpt.at(i)>30&&mainpt.at(i)<=40)
//            {dr=0.0212*factor;}
//            else if (mainpt.at(i)>40&&mainpt.at(i)<=50)
//            {dr=0.0180*factor;}
//            else if (mainpt.at(i)>50&&mainpt.at(i)<=60)
//            {dr=0.0169*factor;}
//            else if (mainpt.at(i)>60&&mainpt.at(i)<=70)
//            {dr=0.0157*factor;}
//            else if (mainpt.at(i)>70&&mainpt.at(i)<=80)
//            {dr=0.0157*factor;}
//            else if (mainpt.at(i)>80&&mainpt.at(i)<=90)
//            {dr=0.0144*factor;}
//            else if (mainpt.at(i)>90&&mainpt.at(i)<=100)
//            {dr=0.0138*factor;}
//            dr=0.8;
            g.SetPtEtaPhiM( looppt.at(j), loopeta.at(j), loopphi.at(j), 0. );
            if(fabs(r.DeltaR(g)) < dr && ptmatch < g.Pt() ){ // find the max pt in cone of interest     // &&
                ptmatch = g.Pt();
                dphimin = r.DeltaPhi(g);
            }
        }
        perobjectmindphi.emplace_back( dphimin );
    }
    //
    return perobjectmindphi;
}

// This one is rather complecated, it takes two bool arguments. use the first argument "tk_eg_match" to tell this function we are interested in l1eg or l1electron. The second argument is somewhat redundant, I was lasy to modify it. Anyway it helps to return all kinds of parameters.
vector<float> dEgamma ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<int> &looplooseTkID, vector<int> &loophgc, vector<float> &loopdphi, vector<float> &loopdeta, vector<float> &loopdz, vector<float> &loopdl, float dr, bool tk_eg_match, bool gen_eg_delta, int delta){
    vector<float> perobjectmindphi;
    TLorentzVector r, g, g1, g2;
    //
    for (std::size_t i = 0; i != maineta.size(); ++i) {
        float dphimin = 9999.;
        float ptmatch   = 0.;
        //float dphicut =dphimatch?0.03:9998;
        //float detacut =dphimatch?0.025:9998;
            //
            for (std::size_t j = 0; j != loopeta.size(); ++j) { // assume this is loop over tracks (doesnt have to be)
                if (looplooseTkID.at(j)==1){
                    if (looppt.at(j)<20) continue;// add pt cut on l1eg
                    if (tk_eg_match){// we require l1eg to match with l1tk?
                        if(loopdphi.at(j)>=9998) continue;
                    }
//                    if (loophgc.at(j)==1) {
//                        float theta = 2*TMath::ATan(TMath::Exp(-(maineta.at(i))));
//                        if(fabs(2.47*0.3*fabs(TMath::Tan(theta)*320./130)*maincharge.at(i)/looppt.at(j))>1) continue;
//                    }
//
//                auto loopphicorr=phiCorr(looppt.at(j), maineta.at(i), mainphi.at(i),  maincharge.at(i), loophgc.at(j), 2.47*0.3);
//                auto loopetacorr=etaCorr(maineta.at(i), mainz0.at(i), loophgc.at(j));
                r.SetPtEtaPhiM( mainpt.at(i), maineta.at(i), mainphi.at(i), 0. );
////                if(phietacorr){
//                r.SetPtEtaPhiM( mainpt.at(i), loopetacorr, loopphicorr, 0. );//assume this is cluster (doesnt have to be)
////                    }
                g.SetPtEtaPhiM( looppt.at(j), loopeta.at(j), loopphi.at(j), 0. );
                if( fabs(r.DeltaR(g)) < dr && ptmatch < g.Pt() ){
                    ptmatch = g.Pt();
                    if (!gen_eg_delta){
                        dphimin = r.DeltaPhi(g);}
                        if(gen_eg_delta&&delta==1){
                            dphimin=loopdphi.at(j);
                        }
                    if(gen_eg_delta&&delta==2){
                        dphimin=loopdeta.at(j);
                    }
                    if(gen_eg_delta&&delta==3){
                        dphimin=loopdz.at(j);
                    }
                    if(gen_eg_delta&&delta==4){
                        dphimin=loopdl.at(j);
                    }
                    if(gen_eg_delta&&delta==5){
                        dphimin=looppt.at(j);
                    }
                    if(gen_eg_delta&&delta==6){
                        dphimin=loopeta.at(j);
                    }
                    if(gen_eg_delta&&delta==7){
                        dphimin=loopphi.at(j);
                    }
                    if(gen_eg_delta&&delta==8){
                        dphimin=loophgc.at(j);
                    }
                }
            }
        }
        perobjectmindphi.emplace_back( dphimin );
    }
    //
    return perobjectmindphi;
}

vector<float> dR ( vector<float> &mainpt, vector<float> &maineta, vector<float> &mainphi, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi){
    vector<float> dr;
    TLorentzVector r, g;
    //
    for (std::size_t i = 0; i != maineta.size(); ++i) {
        float drmin = 0;
        float ptmatch   = 0.;
        r.SetPtEtaPhiM( mainpt.at(i), maineta.at(i), mainphi.at(i), 0. );//assume this is cluster (doesnt have to be)
        for (std::size_t j = 0; j != loopeta.size(); ++j) { // assume this is loop over tracks (doesnt have to be)
            g.SetPtEtaPhiM( looppt.at(j), loopeta.at(j), loopphi.at(j), 0. );
            if(fabs(r.DeltaR(g)) < 0.1 && ptmatch < g.Pt() ){ // find the max pt in cone of interest     // && fabs(g.Pt()-r.Pt())/(g.Pt()+r.Pt())<99999.0 ){
                ptmatch = g.Pt();
                drmin = r.DeltaR(g);
            }
        }
        dr.emplace_back( drmin );
    }
    //
    return dr;
}

// track isolations
float rel_isolation (float &mainpt, float &maineta, float &mainphi, float &mainz0, float &maind0, vector<float> &looppt, vector<float> &loopeta, vector<float> &loopphi, vector<float> &loopz0, vector<float> &loopd0 ){
    TLorentzVector r, g;
    float pt_sum=0;
    float detacut=0.01;
    float dphicut=0.15;
    float dphideta;
    r.SetPtEtaPhiM( mainpt, maineta, mainphi, 0. );
    for (std::size_t j = 0; j != loopeta.size(); ++j) {
        g.SetPtEtaPhiM( looppt.at(j), loopeta.at(j), loopphi.at(j), 0. );
        dphideta=r.DeltaPhi(g)*r.DeltaPhi(g)/(dphicut*dphicut)+(r.Eta()-g.Eta())*(r.Eta()-g.Eta())/(detacut*detacut);
        if (fabs(mainz0-loopz0.at(j))<0.6&&fabs(maind0-loopd0.at(j))<0.5&&r.DeltaR(g)>0.03&&r.DeltaR(g)<0.2){
            pt_sum+=g.Pt();
        }
    }
    return pt_sum/r.Pt();
}
