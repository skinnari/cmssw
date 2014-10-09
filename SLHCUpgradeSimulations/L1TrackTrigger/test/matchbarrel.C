#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;


void matchbarrel(int layer) {

  bool L1L2 = false;
  bool L2L3 = false;
  bool L3L4 = false;
  bool L4L5 = false;
  bool L5L6 = false;

  if (layer==1) L1L2 = true;
  else if (layer==2) L2L3 = true;
  else if (layer==3) L3L4 = true;
  else if (layer==4) L4L5 = true;
  else if (layer==5) L5L6 = true;
  else return;

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1); 
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(1);
  //gStyle->SetOptStat(0);
  

  // ----------------------------------------------------------------------------------  
  // canvas
  // ----------------------------------------------------------------------------------  

  TCanvas* c1 = new TCanvas("c1","Matching to tracklets in barrel",200,10,400,600);
  c1->Divide(2,4);
  c1->SetFillColor(0);
  c1->SetGrid();
  

  // ----------------------------------------------------------------------------------  
  // histograms
  // ----------------------------------------------------------------------------------  

  TH1F *hist11, *hist12, *hist13, *hist14, *hist15, *hist16, *hist17, *hist18;

  if (L1L2) {
    hist11 = new TH1F("h11","Tracklet in L1+2 stub in L3",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in L1+2 stub in L3",100,-10.0,10.0);
    hist13 = new TH1F("h13","Tracklet in L1+2 stub in L4",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in L1+2 stub in L4",100,-10.0,10.0);
    hist15 = new TH1F("h15","Tracklet in L1+2 stub in L5",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in L1+2 stub in L5",100,-10.0,10.0);
    hist17 = new TH1F("h17","Tracklet in L1+2 stub in L6",100,-0.5,0.5);
    hist18 = new TH1F("h18","Tracklet in L1+2 stub in L6",100,-10.0,10.0);
  }
  else if (L2L3) {
    hist11 = new TH1F("h11","Tracklet in L2+3 stub in L1",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in L2+3 stub in L1",100,-10.0,10.0);
    hist13 = new TH1F("h13","Tracklet in L2+3 stub in L4",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in L2+3 stub in L4",100,-10.0,10.0);
    hist15 = new TH1F("h15","Tracklet in L2+3 stub in L5",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in L2+3 stub in L5",100,-10.0,10.0);
    hist17 = new TH1F("h17","Tracklet in L2+3 stub in L6",100,-0.5,0.5);
    hist18 = new TH1F("h18","Tracklet in L2+3 stub in L6",100,-10.0,10.0);
  }
  else if (L3L4) {
    hist11 = new TH1F("h11","Tracklet in L3+4 stub in L1",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in L3+4 stub in L1",100,-10.0,10.0);
    hist13 = new TH1F("h13","Tracklet in L3+4 stub in L2",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in L3+4 stub in L2",100,-10.0,10.0);
    hist15 = new TH1F("h15","Tracklet in L3+4 stub in L5",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in L3+4 stub in L5",100,-10.0,10.0);
    hist17 = new TH1F("h17","Tracklet in L3+4 stub in L6",100,-0.5,0.5);
    hist18 = new TH1F("h18","Tracklet in L3+4 stub in L6",100,-10.0,10.0);
  }
  else if (L4L5) {
    hist11 = new TH1F("h11","Tracklet in L4+5 stub in L1",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in L4+5 stub in L1",100,-10.0,10.0);
    hist13 = new TH1F("h13","Tracklet in L4+5 stub in L2",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in L4+5 stub in L2",100,-10.0,10.0);
    hist15 = new TH1F("h15","Tracklet in L4+5 stub in L3",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in L4+5 stub in L3",100,-10.0,10.0);
    hist17 = new TH1F("h17","Tracklet in L4+5 stub in L6",100,-0.5,0.5);
    hist18 = new TH1F("h18","Tracklet in L4+5 stub in L6",100,-10.0,10.0);
  }
  else if (L5L6) {
    hist11 = new TH1F("h11","Tracklet in L5+6 stub in L1",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in L5+6 stub in L1",100,-10.0,10.0);
    hist13 = new TH1F("h13","Tracklet in L5+6 stub in L2",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in L5+6 stub in L2",100,-10.0,10.0);
    hist15 = new TH1F("h15","Tracklet in L5+6 stub in L3",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in L5+6 stub in L3",100,-10.0,10.0);
    hist17 = new TH1F("h17","Tracklet in L5+6 stub in L4",100,-0.5,0.5);
    hist18 = new TH1F("h18","Tracklet in L5+6 stub in L4",100,-10.0,10.0);
  }

  hist11->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist12->GetXaxis()->SetTitle("#Delta z (cm)");
  hist13->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist14->GetXaxis()->SetTitle("#Delta z (cm)");
  hist15->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist16->GetXaxis()->SetTitle("#Delta z (cm)");
  hist17->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist18->GetXaxis()->SetTitle("#Delta z (cm)");


  // ----------------------------------------------------------------------------------  
  // input file
  // ----------------------------------------------------------------------------------  

  ifstream in("barrelmatch.txt");
  
  int count=0;
  

  // ----------------------------------------------------------------------------------  
  // loop 
  // ----------------------------------------------------------------------------------  

  while (in.good()){
    
    double rtracklet,rstub,deltarphi,deltaz;
    in >>rtracklet>>rstub>>deltarphi>>deltaz;
    
    count++;
    if (count%1000==0) cout << "count="<<count<<" "<<endl;
    

    // ----------------------------------------------------------------------------------  
    // R(L1) = 23 cm
    // R(L2) = 36 cm
    // R(L3) = 51 cm
    // R(L4) = 69 cm
    // R(L5) = 89 cm
    // R(L6) = 109 cm
    // ----------------------------------------------------------------------------------  

    if (L1L2) {
      if (rtracklet<30) { 
	if (fabs(rstub-51)<5) {
	  if (fabs(deltaz)<0.5) hist11->Fill(deltarphi);
	  if (fabs(deltarphi)<0.075) hist12->Fill(deltaz);
	}
	if (fabs(rstub-69)<10) {
	  if (fabs(deltaz)<3.0) hist13->Fill(deltarphi);
	  if (fabs(deltarphi)<0.12) hist14->Fill(deltaz);       
	}
	if (fabs(rstub-89)<10) {
	  if (fabs(deltaz)<3.0) hist15->Fill(deltarphi);
	  if (fabs(deltarphi)<0.15) hist16->Fill(deltaz);       
	}
	if (fabs(rstub-109)<10) {
	  if (fabs(deltaz)<3.0) hist17->Fill(deltarphi);
	  if (fabs(deltarphi)<0.2) hist18->Fill(deltaz);       
	}
      }
    }
    else if (L2L3) {
      if (fabs(rtracklet-36)<5) { 
	if (fabs(rstub-23)<5) {
	  if (fabs(deltaz)<0.5) hist11->Fill(deltarphi);
	  if (fabs(deltarphi)<0.04) hist12->Fill(deltaz);
	}
	if (fabs(rstub-69)<10) {
	  if (fabs(deltaz)<3.0) hist13->Fill(deltarphi);
	  if (fabs(deltarphi)<0.075) hist14->Fill(deltaz);       
	}
	if (fabs(rstub-89)<10) {
	  if (fabs(deltaz)<3.0) hist15->Fill(deltarphi);
	  if (fabs(deltarphi)<0.1) hist16->Fill(deltaz);       
	}
	if (fabs(rstub-109)<10) {
	  if (fabs(deltaz)<3.0) hist17->Fill(deltarphi);
	  if (fabs(deltarphi)<0.15) hist18->Fill(deltaz);       
	}
      }
    }
    else if (L3L4) {
      if (fabs(rtracklet-51)<10) { 
	if (fabs(rstub-23)<5) {
	  if (fabs(deltaz)<5.0) hist11->Fill(deltarphi);
	  if (fabs(deltarphi)<0.025) hist12->Fill(deltaz);
	}
	if (fabs(rstub-36)<5) {
	  if (fabs(deltaz)<3.0) hist13->Fill(deltarphi);
	  if (fabs(deltarphi)<0.025) hist14->Fill(deltaz);       
	}
	if (fabs(rstub-89)<10) {
	  if (fabs(deltaz)<8.0) hist15->Fill(deltarphi);
	  if (fabs(deltarphi)<0.03) hist16->Fill(deltaz);       
	}
	if (fabs(rstub-109)<10) {
	  if (fabs(deltaz)<10.0) hist17->Fill(deltarphi);
	  if (fabs(deltarphi)<0.05) hist18->Fill(deltaz);       
	}
      }
    }
    else if (L4L5) {
      if (fabs(rtracklet-69)<10) { 
	if (fabs(rstub-23)<5) {
	  if (fabs(deltaz)<5.0) hist11->Fill(deltarphi);
	  if (fabs(deltarphi)<0.025) hist12->Fill(deltaz);
	}
	if (fabs(rstub-36)<5) {
	  if (fabs(deltaz)<3.0) hist13->Fill(deltarphi);
	  if (fabs(deltarphi)<0.025) hist14->Fill(deltaz);       
	}
	if (fabs(rstub-51)<5) {
	  if (fabs(deltaz)<8.0) hist15->Fill(deltarphi);
	  if (fabs(deltarphi)<0.03) hist16->Fill(deltaz);       
	}
	if (fabs(rstub-109)<10) {
	  if (fabs(deltaz)<10.0) hist17->Fill(deltarphi);
	  if (fabs(deltarphi)<0.05) hist18->Fill(deltaz);       
	}
      }
    }
    else if (L5L6) {
      if (fabs(rtracklet-89)<10) { 
	if (fabs(rstub-23)<5) {
	  if (fabs(deltaz)<5.0) hist11->Fill(deltarphi);
	  if (fabs(deltarphi)<0.025) hist12->Fill(deltaz);
	}
	if (fabs(rstub-36)<5) {
	  if (fabs(deltaz)<5.0) hist13->Fill(deltarphi);
	  if (fabs(deltarphi)<0.025) hist14->Fill(deltaz);       
	}
	if (fabs(rstub-51)<5) {
	  if (fabs(deltaz)<8.0) hist15->Fill(deltarphi);
	  if (fabs(deltarphi)<0.03) hist16->Fill(deltaz);       
	}
	if (fabs(rstub-69)<10) {
	  if (fabs(deltaz)<10.0) hist17->Fill(deltarphi);
	  if (fabs(deltarphi)<0.05) hist18->Fill(deltaz);       
	}
      }
    }

    
  }
  
  
  // draw
  c1->cd(1);
  h11->Draw();
  c1->cd(2);
  h12->Draw();
  c1->cd(3);
  h13->Draw();
  c1->cd(4);
  h14->Draw();
  c1->cd(5);
  h15->Draw();
  c1->cd(6);
  h16->Draw();
  c1->cd(7);
  h17->Draw();
  c1->cd(8);
  h18->Draw();
  
  if (L1L2) {
    c1->Print("matchbarrel_L1L2.pdf");
    c1->Print("matchbarrel_L1L2.png");
    c1->Print("matchbarrel_L1L2.eps");
  }
  else if (L2L3) {
    c1->Print("matchbarrel_L2L3.pdf");
    c1->Print("matchbarrel_L2L3.png");
    c1->Print("matchbarrel_L2L3.eps");
  }
  else if (L3L4) {
    c1->Print("matchbarrel_L3L4.pdf");
    c1->Print("matchbarrel_L3L4.png");
    c1->Print("matchbarrel_L3L4.eps");
  }
  else if (L4L5) {
    c1->Print("matchbarrel_L4L5.pdf");
    c1->Print("matchbarrel_L4L5.png");
    c1->Print("matchbarrel_L4L5.eps");
  }
  else if (L5L6) {
    c1->Print("matchbarrel_L5L6.pdf");
    c1->Print("matchbarrel_L5L6.png");
    c1->Print("matchbarrel_L5L6.eps");
  }

}

