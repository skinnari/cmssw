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
#include "plotHist2.cc"
#include <string>
#include <map>



void trackletprocessor(){


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

  int ncut1=72*4;
  int ncut2=108*4;


  TCanvas* c1 = new TCanvas("c1","Track performance",200,10,700,800);
  c1->Divide(2,2);
  c1->SetFillColor(0);
  c1->SetGrid();
  
  TCanvas* c2 = new TCanvas("c2","Track performance",200,10,700,800);
  c2->Divide(2,2);
  c2->SetFillColor(0);
  c2->SetGrid();
  
  
  double max=512.0;
  double nbin=100;
  
  TH1 *hist_L1L2 = new TH1F("L1L2","Pairs Tried",nbin,-0.5,max+0.5);
  TH1 *hist_L2L3 = new TH1F("L2L3","Pairs Tried",nbin,-0.5,max+0.5);
  TH1 *hist_L3L4 = new TH1F("L3L4","Pairs Tried",nbin,-0.5,max+0.5);
  TH1 *hist_L5L6 = new TH1F("L5L6","Pairs Tried",nbin,-0.5,max+0.5);

  TH1 *hist_D1D2 = new TH1F("D1D2","Pairs Tried",nbin,-0.5,max+0.5);
  TH1 *hist_D3D4 = new TH1F("D3D4","Pairs Tried",nbin,-0.5,max+0.5);
  TH1 *hist_L1D1 = new TH1F("L1D1","Pairs Tried",nbin,-0.5,max+0.5);
  TH1 *hist_L2D1 = new TH1F("L2D1","Pairs Tried",nbin,-0.5,max+0.5);

  TH1 *hist_L1L2pass = new TH1F("L1L2Pass","Pairs Passed",nbin,-0.5,max+0.5);
  TH1 *hist_L2L3pass = new TH1F("L2L3Pass","Pairs Passed",nbin,-0.5,max+0.5);
  TH1 *hist_L3L4pass = new TH1F("L3L4Pass","Pairs Passed",nbin,-0.5,max+0.5);
  TH1 *hist_L5L6pass = new TH1F("L5L6Pass","Pairs Passed",nbin,-0.5,max+0.5);

  TH1 *hist_D1D2pass = new TH1F("D1D2Pass","Pairs Passed",nbin,-0.5,max+0.5);
  TH1 *hist_D3D4pass = new TH1F("D3D4Pass","Pairs Passed",nbin,-0.5,max+0.5);
  TH1 *hist_L1D1pass = new TH1F("L1D1Pass","Pairs Passed",nbin,-0.5,max+0.5);
  TH1 *hist_L2D1pass = new TH1F("L2D1Pass","Pairs Passed",nbin,-0.5,max+0.5);


  ifstream in("trackletprocessor.txt");

  int count=0;

  std::map<TString, TH1*> hists;

  TString nameold="";
  
  while (in.good()){

    TString name;
    int matchesall;
    int matchespass;
    
    in >>name>>matchesall>>matchespass;

    TSubString sub=name(0,9);
    TString subname(sub);
    cout << "subname "<<subname<<endl;

    if (nameold=="") {
      nameold=subname;
    }

    if (!in.good()) continue;
    
    if (matchesall>max) matchesall=max;
    if (matchespass>max) matchespass=max;
    
    if (name.Contains("_L1L2")) {
      hist_L1L2->Fill(matchesall);
      hist_L1L2pass->Fill(matchespass);
    }
    if (name.Contains("_L2L3")) {
      hist_L2L3->Fill(matchesall);
      hist_L2L3pass->Fill(matchespass);
    }
    if (name.Contains("_L3L4")) {
      hist_L3L4->Fill(matchesall);
      hist_L3L4pass->Fill(matchespass);
    }
    if (name.Contains("_L5L6")) {
      hist_L5L6->Fill(matchesall);
      hist_L5L6pass->Fill(matchespass);
    }

    if (name.Contains("_D1D2")) {
      hist_D1D2->Fill(matchesall);
      hist_D1D2pass->Fill(matchespass);
    }
    if (name.Contains("_D3D4")) {
      hist_D3D4->Fill(matchesall);
      hist_D3D4pass->Fill(matchespass);
    }
    if (name.Contains("_L1D1")) {
      hist_L1D1->Fill(matchesall);
      hist_L1D1pass->Fill(matchespass);
    }
    if (name.Contains("_L2D1")) {
      hist_L2D1->Fill(matchesall);
      hist_L2D1pass->Fill(matchespass);
    }


    std::map<TString, TH1*>::iterator it=hists.find(name);

    if (it==hists.end()) {
      TH1 *hist = new TH1F(name,name,nbin,-0.5,max+0.5);
      hist->Fill(matchesall);
      hists[name]=hist;
      cout << "hists size = "<<hists.size()<<" "<<name<<endl;
    } else {
      hists[name]->Fill(matchesall);
    }

    count++;

 }

  cout << "count = "<<count<<endl;

  bool logscale=true;

  c1->cd(1);
  plotHist2(hist_L1L2,hist_L1L2pass,0.05,ncut1,ncut2,logscale);
  
  c1->cd(2);
  plotHist2(hist_L2L3,hist_L2L3pass,0.05,ncut1,ncut2,logscale);

  c1->cd(3);
  plotHist2(hist_L3L4,hist_L3L4pass,0.05,ncut1,ncut2,logscale);

  c1->cd(4);
  plotHist2(hist_L5L6,hist_L5L6pass,0.05,ncut1,ncut2,logscale);

  c1->Print("trackletprocessor.pdf(","pdf");

  c2->cd(1);
  plotHist2(hist_D1D2,hist_D1D2pass,0.05,ncut1,ncut2,logscale);

  c2->cd(2);
  plotHist2(hist_D3D4,hist_D3D4pass,0.05,ncut1,ncut2,logscale);
  
  c2->cd(3);
  plotHist2(hist_L1D1,hist_L1D1pass,0.05,ncut1,ncut2,logscale);

  c2->cd(4);
  plotHist2(hist_L2D1,hist_L2D1pass,0.05,ncut1,ncut2,logscale);
  

  
  c2->Print("trackletprocessor.pdf","pdf");
  
  int pages=0;

  std::map<TString, TH1*>::iterator it=hists.begin();

  TCanvas* c=0;
  
  
  while(it!=hists.end()) {
    
    if (pages%4==0) {
      
     c = new TCanvas(it->first,"Track performance",200,50,600,700);
     c->Divide(2,2);
     c->SetFillColor(0);
     c->SetGrid();

   }

   c->cd(pages%4+1);
   //gPad->SetLogy();
   plotHist(it->second,0.05,ncut1,ncut2);

   pages++;

   ++it;

   if (it==hists.end()) {
     c->Print("trackletprocessor.pdf)","pdf");
   }
   else {
     if (pages%4==0) {
       c->Print("trackletprocessor.pdf","pdf");
     }
   }

   


 }

  //c->Print("matchengine.pdf)","pdf");

}
