#include "TMath.h"
#include "TRint.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TGaxis.h"
#include <fstream>
#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "plotHist.cc"
#include <string>
#include <map>



void tc_seedpairs(){


  gROOT->Reset();

  gROOT->SetStyle("Plain");

  gStyle->SetCanvasColor(kWhite);

  gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(1111);
  gStyle->SetOptTitle(1);

  // For publishing:
  gStyle->SetLineWidth(1);
  gStyle->SetTextSize(1.1);
  gStyle->SetLabelSize(0.06,"xy");
  gStyle->SetTitleSize(0.06,"xy");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.0,"y");
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.12);


  TCanvas* c1 = new TCanvas("c1","Track performance",200,10,700,800);
  c1->Divide(2,2);
  c1->SetFillColor(0);
  c1->SetGrid();
  
  
  double max=512.0;
  double nbin=100;
  
  TH1 *hist_L1L2 = new TH1F("L1L2","Pairs Tried",nbin,-0.5,max+0.5);

  ifstream in("tc_seedpairs.txt");

  int count=0;

  std::map<TString, std::pair<int,int> > spairs;

  while (in.good()){

    TString name1;
    TString name2;
    int pass;
    
    in >>name1>>name2>>pass;

    TSubString sub1=name1(6,name1.Sizeof()-9);
    TSubString sub2=name2(6,name2.Sizeof()-9);
    //cout << "sub1 "<<sub1<<endl;

    if (!in.good()) continue;

    TString name="SP_"+sub1+"_"+sub2;
    
    std::map<TString, std::pair<int, int> >::iterator it=spairs.find(name);

    if (it==spairs.end()) {
      std::pair<int,int> tmp(1,pass);
      spairs[name]=tmp;
    } else {
      spairs[name].first++;
      spairs[name].second+=pass;
    }

    count++;

 }

  cout << "count = "<<count<<endl;


  std::map<TString, std::pair<int, int> >::iterator it=spairs.begin();

  int countall=0;
  int countunused=0;
  
  while(it!=spairs.end()) {
    countall++;
    if (!it->second.second) {
      cout << it->first << endl;
    }
    //cout << it->first << " " << it->second.first << " " << it->second.second << endl;
    if (it->second.second==0) countunused++;
    it++;
 }

  cout << "Total seedpairs: "<<countall<<" not used: "<<countunused<<endl;

}
