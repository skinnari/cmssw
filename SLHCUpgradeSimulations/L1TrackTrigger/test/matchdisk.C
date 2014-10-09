#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
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


void matchdisk(int layer) {

  bool L1L2 = false;
  bool L2L3 = false;
  bool L3L4 = false;
  bool L4L5 = false;

  if (layer==1) L1L2 = true;
  else if (layer==2) L2L3 = true;
  else if (layer==3) L3L4 = true;
  else if (layer==4) L4L5 = true;
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

  TCanvas* c1 = new TCanvas("c1","Matching to tracklets disks",200,10,400,440);
  c1->Divide(2,3);
  c1->SetFillColor(0);
  c1->SetGrid();

  TCanvas* c2 = new TCanvas("c2","Matching to tracklets disks",200,10,400,440);
  c2->Divide(2,3);
  c2->SetFillColor(0);
  c2->SetGrid();

  TCanvas* c3 = new TCanvas("c3","Matching to tracklets disks",200,10,400,440);
  c3->Divide(2,3);
  c3->SetFillColor(0);
  c3->SetGrid();

  TCanvas* c4 = new TCanvas("c4","Matching to tracklets disks",200,10,400,440);
  c4->Divide(2,3);
  c4->SetFillColor(0);
  c4->SetGrid();
  

  // ----------------------------------------------------------------------------------  
  // histograms
  // ----------------------------------------------------------------------------------  

  TH1F *hist11, *hist12, *hist13, *hist14, *hist15, *hist16;
  TH1F *hist21, *hist22, *hist23, *hist24, *hist25, *hist26;
  TH1F *hist31, *hist32, *hist33, *hist34, *hist35, *hist36;
  TH1F *hist41, *hist42, *hist43, *hist44, *hist45, *hist46;

  if (L1L2) {
    hist11 = new TH1F("h11","Tracklet in D1+2 stub in D3 (PS)",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in D1+2 stub in D3 (PS)",100,-5.0,5.0);
    hist13 = new TH1F("h13","Tracklet in D1+2 stub in D4 (PS)",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in D1+2 stub in D4 (PS)",100,-5.0,5.0);
    hist15 = new TH1F("h15","Tracklet in D1+2 stub in D5 (PS)",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in D1+2 stub in D5 (PS)",100,-5.0,5.0);

    hist31 = new TH1F("h31","Tracklet in D1+2 stub in D3 (2S)",100,-0.5,0.5);
    hist32 = new TH1F("h32","Tracklet in D1+2 stub in D3 (2S)",100,-5.0,5.0);
    hist33 = new TH1F("h33","Tracklet in D1+2 stub in D4 (2S)",100,-0.5,0.5);
    hist34 = new TH1F("h34","Tracklet in D1+2 stub in D4 (2S)",100,-5.0,5.0);
    hist35 = new TH1F("h35","Tracklet in D1+2 stub in D5 (2S)",100,-0.5,0.5);
    hist36 = new TH1F("h36","Tracklet in D1+2 stub in D5 (2S)",100,-5.0,5.0);

    hist21 = new TH1F("h21","Tracklet in outer D1+2 stub in D3 (2S)",100,-0.5,0.5);
    hist22 = new TH1F("h22","Tracklet in outer D1+2 stub in D3 (2S)",100,-5.0,5.0);
    hist23 = new TH1F("h23","Tracklet in outer D1+2 stub in D4 (2S)",100,-0.5,0.5);
    hist24 = new TH1F("h24","Tracklet in outer D1+2 stub in D4 (2S)",100,-5.0,5.0);
    hist25 = new TH1F("h25","Tracklet in outer D1+2 stub in D5 (2S)",100,-0.5,0.5);
    hist26 = new TH1F("h26","Tracklet in outer D1+2 stub in D5 (2S)",100,-5.0,5.0);

    hist41 = new TH1F("h41","Tracklet in outer D1+2 stub in D3 (PS)",100,-0.5,0.5);
    hist42 = new TH1F("h42","Tracklet in outer D1+2 stub in D3 (PS)",100,-5.0,5.0);
    hist43 = new TH1F("h43","Tracklet in outer D1+2 stub in D4 (PS)",100,-0.5,0.5);
    hist44 = new TH1F("h44","Tracklet in outer D1+2 stub in D4 (PS)",100,-5.0,5.0);
    hist45 = new TH1F("h45","Tracklet in outer D1+2 stub in D5 (PS)",100,-0.5,0.5);
    hist46 = new TH1F("h46","Tracklet in outer D1+2 stub in D5 (PS)",100,-5.0,5.0);
  }
  else if (L2L3) {
    hist11 = new TH1F("h11","Tracklet in D2+3 stub in D1 (PS)",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in D2+3 stub in D1 (PS)",100,-5.0,5.0);
    hist13 = new TH1F("h13","Tracklet in D2+3 stub in D4 (PS)",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in D2+3 stub in D4 (PS)",100,-5.0,5.0);
    hist15 = new TH1F("h15","Tracklet in D2+3 stub in D5 (PS)",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in D2+3 stub in D5 (PS)",100,-5.0,5.0);

    hist31 = new TH1F("h31","Tracklet in D2+3 stub in D1 (2S)",100,-0.5,0.5);
    hist32 = new TH1F("h32","Tracklet in D2+3 stub in D1 (2S)",100,-5.0,5.0);
    hist33 = new TH1F("h33","Tracklet in D2+3 stub in D4 (2S)",100,-0.5,0.5);
    hist34 = new TH1F("h34","Tracklet in D2+3 stub in D4 (2S)",100,-5.0,5.0);
    hist35 = new TH1F("h35","Tracklet in D2+3 stub in D5 (2S)",100,-0.5,0.5);
    hist36 = new TH1F("h36","Tracklet in D2+3 stub in D5 (2S)",100,-5.0,5.0);

    hist21 = new TH1F("h21","Tracklet in outer D2+3 stub in D1 (2S)",100,-0.5,0.5);
    hist22 = new TH1F("h22","Tracklet in outer D2+3 stub in D1 (2S)",100,-5.0,5.0);
    hist23 = new TH1F("h23","Tracklet in outer D2+3 stub in D4 (2S)",100,-0.5,0.5);
    hist24 = new TH1F("h24","Tracklet in outer D2+3 stub in D4 (2S)",100,-5.0,5.0);
    hist25 = new TH1F("h25","Tracklet in outer D2+3 stub in D5 (2S)",100,-0.5,0.5);
    hist26 = new TH1F("h26","Tracklet in outer D2+3 stub in D5 (2S)",100,-5.0,5.0);

    hist41 = new TH1F("h41","Tracklet in outer D2+3 stub in D1 (PS)",100,-0.5,0.5);
    hist42 = new TH1F("h42","Tracklet in outer D2+3 stub in D1 (PS)",100,-5.0,5.0);
    hist43 = new TH1F("h43","Tracklet in outer D2+3 stub in D4 (PS)",100,-0.5,0.5);
    hist44 = new TH1F("h44","Tracklet in outer D2+3 stub in D4 (PS)",100,-5.0,5.0);
    hist45 = new TH1F("h45","Tracklet in outer D2+3 stub in D5 (PS)",100,-0.5,0.5);
    hist46 = new TH1F("h46","Tracklet in outer D2+3 stub in D5 (PS)",100,-5.0,5.0);
  }
  else if (L3L4) {
    hist11 = new TH1F("h11","Tracklet in D3+4 stub in D1 (PS)",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in D3+4 stub in D1 (PS)",100,-5.0,5.0);
    hist13 = new TH1F("h13","Tracklet in D3+4 stub in D2 (PS)",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in D3+4 stub in D2 (PS)",100,-5.0,5.0);
    hist15 = new TH1F("h15","Tracklet in D3+4 stub in D5 (PS)",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in D3+4 stub in D5 (PS)",100,-5.0,5.0);

    hist31 = new TH1F("h31","Tracklet in D3+4 stub in D1 (2S)",100,-0.5,0.5);
    hist32 = new TH1F("h32","Tracklet in D3+4 stub in D1 (2S)",100,-5.0,5.0);
    hist33 = new TH1F("h33","Tracklet in D3+4 stub in D2 (2S)",100,-0.5,0.5);
    hist34 = new TH1F("h34","Tracklet in D3+4 stub in D2 (2S)",100,-5.0,5.0);
    hist35 = new TH1F("h35","Tracklet in D3+4 stub in D5 (2S)",100,-0.5,0.5);
    hist36 = new TH1F("h36","Tracklet in D3+4 stub in D5 (2S)",100,-5.0,5.0);

    hist21 = new TH1F("h21","Tracklet in outer D3+4 stub in D1 (2S)",100,-0.5,0.5);
    hist22 = new TH1F("h22","Tracklet in outer D3+4 stub in D1 (2S)",100,-5.0,5.0);
    hist23 = new TH1F("h23","Tracklet in outer D3+4 stub in D2 (2S)",100,-0.5,0.5);
    hist24 = new TH1F("h24","Tracklet in outer D3+4 stub in D2 (2S)",100,-5.0,5.0);
    hist25 = new TH1F("h25","Tracklet in outer D3+4 stub in D5 (2S)",100,-0.5,0.5);
    hist26 = new TH1F("h26","Tracklet in outer D3+4 stub in D5 (2S)",100,-5.0,5.0);

    hist41 = new TH1F("h41","Tracklet in outer D3+4 stub in D1 (PS)",100,-0.5,0.5);
    hist42 = new TH1F("h42","Tracklet in outer D3+4 stub in D1 (PS)",100,-5.0,5.0);
    hist43 = new TH1F("h43","Tracklet in outer D3+4 stub in D2 (PS)",100,-0.5,0.5);
    hist44 = new TH1F("h44","Tracklet in outer D3+4 stub in D2 (PS)",100,-5.0,5.0);
    hist45 = new TH1F("h45","Tracklet in outer D3+4 stub in D5 (PS)",100,-0.5,0.5);
    hist46 = new TH1F("h46","Tracklet in outer D3+4 stub in D5 (PS)",100,-5.0,5.0);
  }
  else if (L4L5) {
    hist11 = new TH1F("h11","Tracklet in D4+5 stub in D1 (PS)",100,-0.5,0.5);
    hist12 = new TH1F("h12","Tracklet in D4+5 stub in D1 (PS)",100,-5.0,5.0);
    hist13 = new TH1F("h13","Tracklet in D4+5 stub in D2 (PS)",100,-0.5,0.5);
    hist14 = new TH1F("h14","Tracklet in D4+5 stub in D2 (PS)",100,-5.0,5.0);
    hist15 = new TH1F("h15","Tracklet in D4+5 stub in D3 (PS)",100,-0.5,0.5);
    hist16 = new TH1F("h16","Tracklet in D4+5 stub in D3 (PS)",100,-5.0,5.0);

    hist31 = new TH1F("h31","Tracklet in D4+5 stub in D1 (2S)",100,-0.5,0.5);
    hist32 = new TH1F("h32","Tracklet in D4+5 stub in D1 (2S)",100,-5.0,5.0);
    hist33 = new TH1F("h33","Tracklet in D4+5 stub in D2 (2S)",100,-0.5,0.5);
    hist34 = new TH1F("h34","Tracklet in D4+5 stub in D2 (2S)",100,-5.0,5.0);
    hist35 = new TH1F("h35","Tracklet in D4+5 stub in D3 (2S)",100,-0.5,0.5);
    hist36 = new TH1F("h36","Tracklet in D4+5 stub in D3 (2S)",100,-5.0,5.0);

    hist21 = new TH1F("h21","Tracklet in outer D4+5 stub in D1 (2S)",100,-0.5,0.5);
    hist22 = new TH1F("h22","Tracklet in outer D4+5 stub in D1 (2S)",100,-5.0,5.0);
    hist23 = new TH1F("h23","Tracklet in outer D4+5 stub in D2 (2S)",100,-0.5,0.5);
    hist24 = new TH1F("h24","Tracklet in outer D4+5 stub in D2 (2S)",100,-5.0,5.0);
    hist25 = new TH1F("h25","Tracklet in outer D4+5 stub in D3 (2S)",100,-0.5,0.5);
    hist26 = new TH1F("h26","Tracklet in outer D4+5 stub in D3 (2S)",100,-5.0,5.0);

    hist41 = new TH1F("h41","Tracklet in outer D4+5 stub in D1 (PS)",100,-0.5,0.5);
    hist42 = new TH1F("h42","Tracklet in outer D4+5 stub in D1 (PS)",100,-5.0,5.0);
    hist43 = new TH1F("h43","Tracklet in outer D4+5 stub in D2 (PS)",100,-0.5,0.5);
    hist44 = new TH1F("h44","Tracklet in outer D4+5 stub in D2 (PS)",100,-5.0,5.0);
    hist45 = new TH1F("h45","Tracklet in outer D4+5 stub in D3 (PS)",100,-0.5,0.5);
    hist46 = new TH1F("h46","Tracklet in outer D4+5 stub in D3 (PS)",100,-5.0,5.0);
  }


  hist11->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist12->GetXaxis()->SetTitle("#Delta r (cm)");
  hist13->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist14->GetXaxis()->SetTitle("#Delta r (cm)");
  hist15->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist16->GetXaxis()->SetTitle("#Delta r (cm)");

  hist21->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist22->GetXaxis()->SetTitle("#Delta r (cm)");
  hist23->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist24->GetXaxis()->SetTitle("#Delta r (cm)");
  hist25->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist26->GetXaxis()->SetTitle("#Delta r (cm)");

  hist31->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist32->GetXaxis()->SetTitle("#Delta r (cm)");
  hist33->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist34->GetXaxis()->SetTitle("#Delta r (cm)");
  hist35->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist36->GetXaxis()->SetTitle("#Delta r (cm)");


  hist41->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist42->GetXaxis()->SetTitle("#Delta r (cm)");
  hist43->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist44->GetXaxis()->SetTitle("#Delta r (cm)");
  hist45->GetXaxis()->SetTitle("#Delta r#phi (cm)");
  hist46->GetXaxis()->SetTitle("#Delta r (cm)");


  // ----------------------------------------------------------------------------------  
  // input file
  // ----------------------------------------------------------------------------------  

  ifstream in("diskmatch.txt");
  
  int count=0;
  

  // ----------------------------------------------------------------------------------  
  // loop 
  // ----------------------------------------------------------------------------------  

  while (in.good()){
    
    double rtracklet,ztracklet,rstub,zstub,deltarphi,deltar;
    in >>rtracklet>>ztracklet>>rstub>>zstub>>deltarphi>>deltar;

    count++;
    if (count%1000==0) cout << "count="<<count<<" "<<endl;
    

    // ----------------------------------------------------------------------------------  
    // z(D1) = 135 cm
    // z(D2) = 160 cm
    // z(D3) = 189 cm
    // z(D4) = 224 cm
    // z(D5) = 265 cm
    // ----------------------------------------------------------------------------------  

    // ----------------------------------------------------------------------------------  
    // L1-L2
    // ----------------------------------------------------------------------------------  
    if (L1L2) {
      if (fabs(fabs(ztracklet)-135)<10 && rtracklet<60.) {
	if (rstub<60.) { //seed in PS modules, project to PS
	  if (fabs(fabs(zstub)-189)<10) {
	    if (fabs(deltar)<3.0) hist11->Fill(deltarphi);
	    if (fabs(deltarphi)<0.12) hist12->Fill(deltar);
	  }
	  else if (fabs(fabs(zstub)-224)<10) {
	    if (fabs(deltar)<3.0) hist13->Fill(deltarphi);
	    if (fabs(deltarphi)<0.15) hist14->Fill(deltar);
	  }
	  else if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<3.0) hist15->Fill(deltarphi);
	    if (fabs(deltarphi)<0.2) hist16->Fill(deltar);
	  }
	}
	else { //seed in PS modules, project to 2S
	  if (fabs(fabs(zstub)-189)<10) {
	    if (fabs(deltar)<3.0) hist31->Fill(deltarphi);
	    if (fabs(deltarphi)<0.12) hist32->Fill(deltar);
	  }
	  else if (fabs(fabs(zstub)-224)<10) {
	    if (fabs(deltar)<3.0) hist33->Fill(deltarphi);
	    if (fabs(deltarphi)<0.15) hist34->Fill(deltar);
	  }
	  else if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<3.0) hist35->Fill(deltarphi);
	    if (fabs(deltarphi)<0.2) hist36->Fill(deltar);
	  }
	}
      }
      else if (fabs(fabs(ztracklet)-135)<10 && rtracklet>60.) { 
	if (rstub>60.) { //seed in 2S modules, project to 2S
	  if (fabs(fabs(zstub)-189)<10) {
	    if (fabs(deltar)<5.0) hist21->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist22->Fill(deltar);
	  }
	  else if (fabs(fabs(zstub)-224)<10) {
	    if (fabs(deltar)<5.0) hist23->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist24->Fill(deltar);
	  }
	  else if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<5.0) hist25->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist26->Fill(deltar);
	  }
	}
	else { //seed in 2S modules, project to PS
	  if (fabs(fabs(zstub)-189)<10) {
	    if (fabs(deltar)<5.0) hist41->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist42->Fill(deltar);
	  }
	  else if (fabs(fabs(zstub)-224)<10) {
	    if (fabs(deltar)<5.0) hist43->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist44->Fill(deltar);
	  }
	  else if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<5.0) hist45->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist46->Fill(deltar);
	  }
	}
      }
    }
    // ----------------------------------------------------------------------------------  
    // L2-L3
    // ----------------------------------------------------------------------------------  
    else if (L2L3) {
      if (fabs(fabs(ztracklet)-160)<10 && rtracklet<60.) {
	if (rstub<60.) { //seed in PS modules, project to PS
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<3.0) hist11->Fill(deltarphi);
	    if (fabs(deltarphi)<0.12) hist12->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-224)<10) {
	    if (fabs(deltar)<3.0) hist13->Fill(deltarphi);
	    if (fabs(deltarphi)<0.15) hist14->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<3.0) hist15->Fill(deltarphi);
	    if (fabs(deltarphi)<0.2) hist16->Fill(deltar);
	  }
	}
	else { //seed in PS modules, project to 2S
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<3.0) hist31->Fill(deltarphi);
	    if (fabs(deltarphi)<0.12) hist32->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-224)<10) {
	    if (fabs(deltar)<3.0) hist33->Fill(deltarphi);
	    if (fabs(deltarphi)<0.15) hist34->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<3.0) hist35->Fill(deltarphi);
	    if (fabs(deltarphi)<0.2) hist36->Fill(deltar);
	  }
	}
      }
      else if (fabs(fabs(ztracklet)-160)<10 && rtracklet>60.) { 
	if (rstub>60.) { //seed in 2S modules, project to 2S
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<5.0) hist21->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist22->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-224)<10) {
	    if (fabs(deltar)<5.0) hist23->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist24->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<5.0) hist25->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist26->Fill(deltar);
	  }
	}
	else { //seed in 2S modules, project to PS
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<5.0) hist41->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist42->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-224)<10) {
	    if (fabs(deltar)<5.0) hist43->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist44->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<5.0) hist45->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist46->Fill(deltar);
	  }
	}
      }
    }
    // ----------------------------------------------------------------------------------  
    // L3-L4
    // ----------------------------------------------------------------------------------  
    else if (L3L4) {
      if (fabs(fabs(ztracklet)-189)<10 && rtracklet<60.) {
	if (rstub<60.) { //seed in PS modules, project to PS
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<3.0) hist11->Fill(deltarphi);
	    if (fabs(deltarphi)<0.12) hist12->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-160)<10) {
	    if (fabs(deltar)<3.0) hist13->Fill(deltarphi);
	    if (fabs(deltarphi)<0.15) hist14->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<3.0) hist15->Fill(deltarphi);
	    if (fabs(deltarphi)<0.2) hist16->Fill(deltar);
	  }
	}
	else { //seed in PS modules, project to 2S
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<3.0) hist31->Fill(deltarphi);
	    if (fabs(deltarphi)<0.12) hist32->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-160)<10) {
	    if (fabs(deltar)<3.0) hist33->Fill(deltarphi);
	    if (fabs(deltarphi)<0.15) hist34->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<3.0) hist35->Fill(deltarphi);
	    if (fabs(deltarphi)<0.2) hist36->Fill(deltar);
	  }
	}
      }
      else if (fabs(fabs(ztracklet)-189)<10 && rtracklet>60.) { 
	if (rstub>60.) { //seed in 2S modules, project to 2S
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<5.0) hist21->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist22->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-160)<10) {
	    if (fabs(deltar)<5.0) hist23->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist24->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<5.0) hist25->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist26->Fill(deltar);
	  }
	}
	else { //seed in 2S modules, project to PS
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<5.0) hist41->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist42->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-160)<10) {
	    if (fabs(deltar)<5.0) hist43->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist44->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-265)<10) {
	    if (fabs(deltar)<5.0) hist45->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist46->Fill(deltar);
	  }
	}
      }
    }
    // ----------------------------------------------------------------------------------  
    // L4-L5
    // ----------------------------------------------------------------------------------  
    else if (L4L5) {
      if (fabs(fabs(ztracklet)-224)<10 && rtracklet<60.) {
	if (rstub<60.) { //seed in PS modules, project to PS
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<3.0) hist11->Fill(deltarphi);
	    if (fabs(deltarphi)<0.12) hist12->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-160)<10) {
	    if (fabs(deltar)<3.0) hist13->Fill(deltarphi);
	    if (fabs(deltarphi)<0.15) hist14->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-189)<10) {
	    if (fabs(deltar)<3.0) hist15->Fill(deltarphi);
	    if (fabs(deltarphi)<0.2) hist16->Fill(deltar);
	  }
	}
	else { //seed in PS modules, project to 2S
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<3.0) hist31->Fill(deltarphi);
	    if (fabs(deltarphi)<0.12) hist32->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-160)<10) {
	    if (fabs(deltar)<3.0) hist33->Fill(deltarphi);
	    if (fabs(deltarphi)<0.15) hist34->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-189)<10) {
	    if (fabs(deltar)<3.0) hist35->Fill(deltarphi);
	    if (fabs(deltarphi)<0.2) hist36->Fill(deltar);
	  }
	}
      }
      else if (fabs(fabs(ztracklet)-224)<10 && rtracklet>60.) { 
	if (rstub>60.) { //seed in 2S modules, project to 2S
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<5.0) hist21->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist22->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-160)<10) {
	    if (fabs(deltar)<5.0) hist23->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist24->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-189)<10) {
	    if (fabs(deltar)<5.0) hist25->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist26->Fill(deltar);
	  }
	}
	else { //seed in 2S modules, project to PS
	  if (fabs(fabs(zstub)-135)<10) {
	    if (fabs(deltar)<5.0) hist41->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist42->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-160)<10) {
	    if (fabs(deltar)<5.0) hist43->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist44->Fill(deltar);
	  }
	  if (fabs(fabs(zstub)-189)<10) {
	    if (fabs(deltar)<5.0) hist45->Fill(deltarphi);
	    if (fabs(deltarphi)<0.5) hist46->Fill(deltar);
	  }
	}
      }
    }

    
  }
  
  
  // draw
  c1->cd(1);
  hist11->Draw();
  c1->cd(2);
  hist12->Draw();
  c1->cd(3);
  hist13->Draw();
  c1->cd(4);
  hist14->Draw();
  c1->cd(5);
  hist15->Draw();
  c1->cd(6);
  hist16->Draw();

  c2->cd(1);
  hist21->Draw();
  c2->cd(2);
  hist22->Draw();
  c2->cd(3);
  hist23->Draw();
  c2->cd(4);
  hist24->Draw();
  c2->cd(5);
  hist25->Draw();
  c2->cd(6);
  hist26->Draw();

  c3->cd(1);
  hist31->Draw();
  c3->cd(2);
  hist32->Draw();
  c3->cd(3);
  hist33->Draw();
  c3->cd(4);
  hist34->Draw();
  c3->cd(5);
  hist35->Draw();
  c3->cd(6);
  hist36->Draw();

  c4->cd(1);
  hist41->Draw();
  c4->cd(2);
  hist42->Draw();
  c4->cd(3);
  hist43->Draw();
  c4->cd(4);
  hist44->Draw();
  c4->cd(5);
  hist45->Draw();
  c4->cd(6);
  hist46->Draw();
  
  if (L1L2) {
    c1->Print("matchdiskPStoPS_L1L2.pdf");
    c1->Print("matchdiskPStoPS_L1L2.png");
    c1->Print("matchdiskPStoPS_L1L2.eps");
    c3->Print("matchdiskPSto2S_L1L2.pdf");
    c3->Print("matchdiskPSto2S_L1L2.png");
    c3->Print("matchdiskPSto2S_L1L2.eps");
    c2->Print("matchdisk2Sto2S_L1L2.pdf");
    c2->Print("matchdisk2Sto2S_L1L2.png");
    c2->Print("matchdisk2Sto2S_L1L2.eps");
    c4->Print("matchdisk2StoPS_L1L2.pdf");
    c4->Print("matchdisk2StoPS_L1L2.png");
    c4->Print("matchdisk2StoPS_L1L2.eps");
  }
  else if (L2L3) {
    c1->Print("matchdiskPStoPS_L2L3.pdf");
    c1->Print("matchdiskPStoPS_L2L3.png");
    c1->Print("matchdiskPStoPS_L2L3.eps");
    c3->Print("matchdiskPSto2S_L2L3.pdf");
    c3->Print("matchdiskPSto2S_L2L3.png");
    c3->Print("matchdiskPSto2S_L2L3.eps");
    c2->Print("matchdisk2Sto2S_L2L3.pdf");
    c2->Print("matchdisk2Sto2S_L2L3.png");
    c2->Print("matchdisk2Sto2S_L2L3.eps");
    c4->Print("matchdisk2StoPS_L2L3.pdf");
    c4->Print("matchdisk2StoPS_L2L3.png");
    c4->Print("matchdisk2StoPS_L2L3.eps");
  }
  else if (L3L4) {
    c1->Print("matchdiskPStoPS_L3L4.pdf");
    c1->Print("matchdiskPStoPS_L3L4.png");
    c1->Print("matchdiskPStoPS_L3L4.eps");
    c3->Print("matchdiskPSto2S_L3L4.pdf");
    c3->Print("matchdiskPSto2S_L3L4.png");
    c3->Print("matchdiskPSto2S_L3L4.eps");
    c2->Print("matchdisk2Sto2S_L3L4.pdf");
    c2->Print("matchdisk2Sto2S_L3L4.png");
    c2->Print("matchdisk2Sto2S_L3L4.eps");
    c4->Print("matchdisk2StoPS_L3L4.pdf");
    c4->Print("matchdisk2StoPS_L3L4.png");
    c4->Print("matchdisk2StoPS_L3L4.eps");
  }
  else if (L4L5) {
    c1->Print("matchdiskPStoPS_L4L5.pdf");
    c1->Print("matchdiskPStoPS_L4L5.png");
    c1->Print("matchdiskPStoPS_L4L5.eps");
    c3->Print("matchdiskPSto2S_L4L5.pdf");
    c3->Print("matchdiskPSto2S_L4L5.png");
    c3->Print("matchdiskPSto2S_L4L5.eps");
    c2->Print("matchdisk2Sto2S_L4L5.pdf");
    c2->Print("matchdisk2Sto2S_L4L5.png");
    c2->Print("matchdisk2Sto2S_L4L5.eps");
    c4->Print("matchdisk2StoPS_L4L5.pdf");
    c4->Print("matchdisk2StoPS_L4L5.png");
    c4->Print("matchdisk2StoPS_L4L5.eps");
  }

}

