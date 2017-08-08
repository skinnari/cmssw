//This class implementes the match calculator
#ifndef FPGAMATCHCALCULATOR_H
#define FPGAMATCHCALCULATOR_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAMatchCalculator:public FPGAProcessBase{

public:

  FPGAMatchCalculator(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
    
    fullmatchesToPlus_=0;
    fullmatchesToMinus_=0;
    double dphi=two_pi/NSector;
    phimin_=iSector*dphi;
    phimax_=phimin_+dphi;
    if (phimin_>0.5*two_pi) phimin_-=two_pi;
    if (phimax_>0.5*two_pi) phimax_-=two_pi;
    phioffset_=phimin_-(dphi)/6.0;
    string subname=name.substr(8,2);
    fullmatchesToPlus_=0;
    fullmatchesToMinus_=0;
    layer_=0;
    if (subname=="L1") layer_=1;
    if (subname=="L2") layer_=2;
    if (subname=="L3") layer_=3;
    if (subname=="L4") layer_=4;
    if (subname=="L5") layer_=5;
    if (subname=="L6") layer_=6;
    disk_=0;
    if (subname=="F1") disk_=1;
    if (subname=="F2") disk_=2;
    if (subname=="F3") disk_=3;
    if (subname=="F4") disk_=4;
    if (subname=="F5") disk_=5;
    if (subname=="D1") disk_=1;
    if (subname=="D2") disk_=2;
    if (subname=="D3") disk_=3;
    if (subname=="D4") disk_=4;
    if (subname=="D5") disk_=5;
    if (subname=="B1") disk_=-1;
    if (subname=="B2") disk_=-2;
    if (subname=="B3") disk_=-3;
    if (subname=="B4") disk_=-4;
    if (subname=="B5") disk_=-5;
    if (layer_==0 && disk_==0) {
      cout << "name subname "<<name<<" "<<subname<<endl;
      assert(0);
    }
    
  }

  void addOutput(FPGAMemoryBase* memory,string output){
    if (writetrace) {
      cout << "In "<<name_<<" adding output to "<<memory->getName()
	   << " to output "<<output<<endl;
    }
    if (output=="matchout1"||
	output=="matchout2"||
	output=="matchout3"||
	output=="matchout4"||
	output=="matchout5"||
	output=="matchout6"||
	output=="matchout7"||
	output=="matchout8"
	){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      fullmatches_.push_back(tmp);
      return;
    }
    if (output=="matchoutplus"){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      assert(fullmatchesToPlus_==0);
      fullmatchesToPlus_=tmp;
      return;
    }
    if (output=="matchoutminus"){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      assert(fullmatchesToMinus_==0);
      fullmatchesToMinus_=tmp;
      return;
    }
    cout << "Count not fined output = "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="allstubin"){
      FPGAAllStubs* tmp=dynamic_cast<FPGAAllStubs*>(memory);
      assert(tmp!=0);
      allstubs_=tmp;
      return;
    }
    if (input=="allprojin"){
      FPGAAllProjections* tmp=dynamic_cast<FPGAAllProjections*>(memory);
      assert(tmp!=0);
      allprojs_=tmp;
      return;
    }
    if (input=="match1in"||
	input=="match2in"||
	input=="match3in"||
	input=="match4in"||
	input=="match5in"||
	input=="match6in"||
	input=="match7in"||
	input=="match8in"){
      FPGACandidateMatch* tmp=dynamic_cast<FPGACandidateMatch*>(memory);
      assert(tmp!=0);
      matches_.push_back(tmp);
      return;
    }
    assert(0);
  }

  std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > mergeMatches(vector<FPGACandidateMatch*>& candmatch) {

    std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > >  tmp;

	//cout << "FPGAMatchCalcultor::mergedMatches inputs"<<endl;
	std::vector<unsigned int> indexArray;
	indexArray.reserve(candmatch.size());
	for (unsigned int i=0;i<candmatch.size();i++) {
	  indexArray.push_back(0);
	}

	int bestIndex=-1;
	do {
      int bestSector=100;
	  int bestTCID=(1<<16);
	  bestIndex=-1;
	  for (unsigned int i=0;i<candmatch.size();i++) {
	    if (indexArray[i]>=candmatch[i]->nMatches()) {
		  // skip as we were at the end
		  continue;
		}
		int TCID=candmatch[i]->getFPGATracklet(indexArray[i])->TCID();
		int dSector=candmatch[i]->getFPGATracklet(indexArray[i])->homeSector()-iSector_;
		if (dSector>2) dSector-=NSector;
		if (dSector<-2) dSector+=NSector;
		assert(abs(dSector)<2);
		if (dSector==-1) dSector=2;
		if (dSector<bestSector) {
		  bestSector=dSector;
		  bestTCID=TCID;
		  bestIndex=i;
		}
		if (dSector==bestSector) {
		  if (TCID<bestTCID) {
		    bestTCID=TCID;
			bestIndex=i;
		  }
		}
	  }
	  if (bestIndex!=-1) {
	    tmp.push_back(candmatch[bestIndex]->getMatch(indexArray[bestIndex]));
		indexArray[bestIndex]++;
	  }
	} while (bestIndex!=-1);

	//cout << "In FPGAMatchCalculator::mergedMatches : "<<tmp.size()<<endl;
	for (unsigned int i=0;i<tmp.size();i++) {
	  //cout <<iSector_<<" "<<tmp[i].first->homeSector()<<" "<<tmp[i].first->TCID()<<endl;
	  if (i>0) {
	    //This allows for equal TCIDs. This means that we can e.g. have a track seeded
	    //in L1L2 that projects to both L3 and D4. The algorithm will pick up the first hit and
		//drop the second
		//cout << "Tracklet sector : "<<getName()<<" "<<iSector_<<" "
		//	 <<tmp[i-1].first->homeSector()<<" "<<tmp[i].first->homeSector()<<endl;

	    if (iSector_==tmp[i-1].first->homeSector()&&
			iSector_==tmp[i].first->homeSector()) {
		  assert(tmp[i-1].first->TCID()<=tmp[i].first->TCID());
		}
	  }
	}
	//cout << endl;

	return tmp;
  }
	
  void execute() {

    
    //again here we will cheat a little and use the information in matches

    assert(fullmatches_.size()!=0);
    //assert(fullmatchesToPlus_!=0);
    //assert(fullmatchesToMinus_!=0);

    unsigned int countall=0;
    unsigned int countsel=0;

    FPGATracklet* oldTracklet=0;

    std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > mergedMatches=mergeMatches(matches_);
    
    //for(unsigned int j=0;j<matches_.size();j++){
      //cout << "FPGAMatchCalculator::execute matches_.size() "<<matches_[j]->nMatches()<<endl;
      //for(unsigned int i=0;i<matches_[j]->nMatches();i++){

    for(unsigned int j=0;j<mergedMatches.size();j++){
	
      //if (debug1) {
      //  cout << getName() <<" has match from "<<matches_[j]->getName()<<endl;
      //}
	
	countall++;
	
	L1TStub* stub=mergedMatches[j].second.second;
	FPGAStub* fpgastub=mergedMatches[j].second.first;
	FPGATracklet* tracklet=mergedMatches[j].first;

	if (oldTracklet!=0) {
	  //cout << "Tracklet sector two : "<<iSector_<<" "
	  //     <<tracklet->homeSector()<<" "<<oldTracklet->homeSector()<<endl;
	  //allow equal here since we can have more than one cadidate match per tracklet projection
		if (iSector_==tracklet->homeSector()&&
			iSector_==oldTracklet->homeSector()) {
		  assert(oldTracklet->TCID()<=tracklet->TCID());
		}
	}
	oldTracklet=tracklet;
	
	if (layer_!=0) {
	  
	  double pttracklet=0.3*3.8/(tracklet->rinv()*100);
	  bool keep=fabs(1.0/pttracklet-1.0/stub->pt())<ptstubconsistencymatching;
            
	  if (!keep) continue;
      
	  int seedlayer=tracklet->layer();
	
	  double phi=stub->phi();
	  if (phi<0) phi+=two_pi;
	  phi-=phioffset_;
	  
	  double dr=stub->r()-tracklet->rproj(layer_);
	  //cout << "stub r and rproj :"<<stub->r()<<" "<<tracklet->rproj(layer_)<<" "<<layer_<<" "<<getName()<<" "<<seedlayer<<" "<<tracklet->disk()<<" tracklet "<<tracklet<<endl;
	  assert(fabs(dr)<drmax);

	  int ir=fpgastub->r().value();
	
	  int iphi=tracklet->fpgaphiproj(layer_).value();
	  int icorr=0;
	  
	  int lr=1;    

	  if (layer_<4) {
	    lr=2;
	    icorr=(ir*tracklet->fpgaphiprojder(layer_).value())>>(idrinvbits+2+phi0bitshift-rinvbitshift-phiderbitshift);
	    assert(idrinvbits+1+phi0bitshift-rinvbitshift-phiderbitshift>=0);
	  } else {
	    int shift=idrinvbits+phi0bitshift-rinvbitshift-phiderbitshift;
	    if (shift>=0){
	      icorr=(ir*tracklet->fpgaphiprojder(layer_).value())>>shift;
	    } else {
	      icorr=(ir*tracklet->fpgaphiprojder(layer_).value())<<(-shift);
	    }
	    //cout <<"idrinvbits phi0bitshift rinvbitshift phiderbitshift: "
	    //	 <<idrinvbits<<" "<<phi0bitshift<<" "
	    //	 <<rinvbitshift<<" "<<phiderbitshift<<endl;
	    //assert(idrinvbits-1+phi0bitshift-rinvbitshift-phiderbitshift>=0);
	  }             
	  
	  iphi+=icorr;
	    
	  int iz=tracklet->fpgazproj(layer_).value();

	  int izcor=0;

	  if (layer_<4) {
	    izcor=(ir*tracklet->fpgazprojder(layer_).value())>>(idrinvbits-1-zderbitshift-tbitshift);
	  } else {
	    izcor=(ir*tracklet->fpgazprojder(layer_).value())>>(idrinvbits-zderbitshift-tbitshift+(nbitszprojL123-nbitszprojL456));
	  }

	  //izcor>>=1; ?? not needed from binarytrackletOverlap
	  
	  iz+=izcor;

	  int ideltaz=fpgastub->z().value()-iz;


	  int ideltaphi=0;
	  if (layer_<4) {  //FIXME are we loosing precission here? Yes!!!!
	    ideltaphi=(fpgastub->phi().value()<<3)-(iphi<<(phi0bitshift+2)); 
	  } else {
	    assert(phi0bitshift-1>=0);
	    ideltaphi=fpgastub->phi().value()-(iphi<<(phi0bitshift-1));  
	  }


	  if (dumpmatch) {
	    int l=1;
            int fact=1;
	    if (layer_>3) {
	      l=8; 
	      fact=(1<<(nbitszprojL123-nbitszprojL456));
	    }

	    cout << "DUMPMATCH : layer = "<<layer_<<endl;
	    cout << "DUMPMATCH : phiproj  = "<<tracklet->phiproj(layer_)<<" "
		 << tracklet->fpgaphiproj(layer_).value()*kphiproj123/l<<endl;
	    cout << "DUMPMATCH : phi corr = "<< dr*tracklet->phiprojder(layer_)
		 << " "<<icorr*kphiproj123/l<<endl;
	    cout << "DUMPMATCH : dr       = "<<dr<<" "<<ir*kr*lr<<endl;
	    cout << "DUMPMATCH : zproj    = "<<tracklet->zproj(layer_)<<" "
		 <<  tracklet->fpgazproj(layer_).value()*kz*fact<<" "
		 <<  tracklet->zproj(layer_)-tracklet->fpgazproj(layer_).value()*kz*fact
		 << endl;
	    cout << "DUMPMATCH : z der   = "<<layer_<<" "<<tracklet->zprojder(layer_)
		 << " "<<tracklet->fpgazprojder(layer_).value()*kz/(kr*(1<<(idrinvbits-zderbitshift-tbitshift)))
		 <<" ("<<tracklet->fpgazprojder(layer_).value()<<")"<<endl;  
	    cout << "DUMPMATCH : z icorr corr  = "<<layer_<< " "<<izcor*kz*fact<<" "<<dr*tracklet->zprojder(layer_)<<endl;
	    double phistub=stub->phi();
	    double dphi=two_pi/NSector;
	    phistub+=dphi/6.0;
	    for(;phistub>dphi;phistub-=dphi);
	    for(;phistub<0;phistub+=dphi);
	    cout << "DUMPMATCH : phi stub project   = "<<phistub<<" "<<tracklet->phiproj(layer_)<<endl;
	    cout << "DUMPMATCH : deltaphi   = "<<ideltaphi*kphi1*rmean[layer_-1]*10000.0<<" "
		 << " "<<(phistub-(tracklet->phiproj(layer_)+dr*tracklet->phiprojder(layer_)))*rmean[layer_-1]*10000.0<<endl;
	    cout << "DUMPMATCH : ideltaz deltaz   = "<<layer_<<" "<<ideltaz*kz<<" "<<stub->z()-(tracklet->zproj(layer_)+
										 dr*tracklet->zprojder(layer_))<<endl;
		 
	  }
	  
	  assert(fabs(dr)<drmax);

	  
	  double dphi=phi-(tracklet->phiproj(layer_)+
			   dr*tracklet->phiprojder(layer_));

	  if (dphi>0.5*two_pi) dphi-=two_pi;
	  if (dphi<-0.5*two_pi) dphi+=two_pi;

	  
	  double dz=stub->z()-(tracklet->zproj(layer_)+
			       dr*tracklet->zprojder(layer_));

	  //if (layer_<4) {
	  //  cout << "Delta z :"<<layer_<<" "<<dz<<" "<<ideltaz*kz<<" "<<dr*tracklet->phiprojder(layer_)
	  //	 <<" "<<izcor*kz<<endl;
	  //}else{
	  //  cout << "Delta z :"<<layer_<<" "<<dz<<" "<<ideltaz*kz*16<<endl;
	  //}
	  
	  double dphiapprox=phi-(tracklet->phiprojapprox(layer_)+
				 dr*tracklet->phiprojderapprox(layer_));

	  if (dphiapprox>0.5*two_pi) dphiapprox-=two_pi;
	  if (dphiapprox<-0.5*two_pi) dphiapprox+=two_pi;


	  double dzapprox=stub->z()-(tracklet->zprojapprox(layer_)+
				     dr*tracklet->zprojderapprox(layer_));
	  
	  bool imatch=false;


	  if (layer_<4){
	    if (seedlayer==1) {
	      imatch=(fabs(ideltaphi)<0.1/(kphi1*rmean[layer_-1]))&&(fabs(ideltaz)<0.5/kz);
	    }
	    else {
	      imatch=(fabs(ideltaphi)<0.095/(kphi1*rmean[layer_-1]))&&(fabs(ideltaz)<14.0/kz);
	    }

	  }else{
	    int fact=(1<<(nbitszprojL123-nbitszprojL456));
	    if (seedlayer==1) {
	      imatch=(fabs(ideltaphi)<0.30/(kphi1*rmean[layer_-1]))&&(fabs(fact*ideltaz)<4.0/kz);
	      if (debug1) {  
		cout << "ideltaphi ideltaz :"<<ideltaphi<<" "<<ideltaz<<endl;
		cout << "deltaphi :"<<ideltaphi*kphi1*rmean[layer_-1]<<endl;
		cout << "deltaz :"<<fact*ideltaz*kz<<endl;
		cout << "imatch : "<<imatch<<endl;
	      }
	    }
	    else {
	      imatch=(fabs(ideltaphi)<0.2/(kphi1*rmean[layer_-1]))&&(fabs(fact*ideltaz)<14.0/kz);
	    }
	  }
	  
	  
	  if (imatch) {
	    //cout << "1: ideltaphi = "<<ideltaphi<<" "<<ideltaphi*kphi1*stub->r()<<endl;
	    std::pair<FPGAStub*,L1TStub*> tmp(fpgastub,stub);

	    countsel++;

	    tracklet->addMatch(layer_,ideltaphi,ideltaz,
			       dphi,dz,dphiapprox,dzapprox,
			       fpgastub->stubindex().value()+(1<<6)*(fpgastub->fedregion()-1),
			       stub->r(),tmp);
	    

	    if (debug1) {
	      cout << "Accepted full match in layer " <<getName()
		   << " "<<tracklet
		   << " "<<iSector_<<endl;	   
	    }
	      
	    if (tracklet->plusNeighbor(layer_)){
	      //cout << "FPGAMatchCaclulator "<<getName()<<" adding match to: "
	      //   <<fullmatchesToMinus_->getName()<<endl;
	      //cout << "Accepted full match to minus in layer " <<getName()<<" "<<tracklet
	      //	   <<" "<<fullmatchesToMinus_->getName()<<endl;
	      assert(fullmatchesToMinus_!=0);
		  
		  int nmatch=fullmatchesToMinus_->nMatches();
		  if (nmatch>1) {
		    //cout << "FPGAMatchCalculator nmatch = "<<nmatch<<endl;
		    assert(fullmatchesToMinus_->getFPGATracklet(nmatch-2)->TCID()<
				   fullmatchesToMinus_->getFPGATracklet(nmatch-1)->TCID());
		  }
		  
	      fullmatchesToMinus_->addMatch(tracklet,tmp);
	    } else if (tracklet->minusNeighbor(layer_)) {
	      //cout << "FPGAMatchCaclulator "<<getName()<<" adding match to: "
	      //   <<fullmatchesToPlus_->getName()<<endl;
		  // FIXME? fullmatchesToPlus_ ?
	      assert(fullmatchesToMinus_!=0);
	      fullmatchesToPlus_->addMatch(tracklet,tmp);
	      //cout << "Accepted full match to plus in layer " <<getName()<<" "<<tracklet
	      //	   <<" "<<fullmatchesToPlus_->getName()<<endl;
	    } else {
	      for (unsigned int l=0;l<fullmatches_.size();l++){
		//cout << "Adding match to: "<<fullmatches_[l]->getName()<<endl;
		if ((tracklet->layer()==1&&fullmatches_[l]->getName().substr(3,2)=="L1")||
		    (tracklet->layer()==3&&fullmatches_[l]->getName().substr(3,2)=="L3")||
		    (tracklet->layer()==5&&fullmatches_[l]->getName().substr(3,2)=="L5")){
		  //cout << "Accepted full match  in layer " <<getName()<<" "<<tracklet->layer()
		  //     <<" "<<fullmatches_[l]->getName()<<endl;
		  assert(tracklet->homeSector()==iSector_);
		  fullmatches_[l]->addMatch(tracklet,tmp);
		}
	      }
	    }


	  }
	  
	} else {


	  //Disk matches are calculated here
	  //hack for now to skip stubs in the wrong (F/B) part of the detector
	  if (stub->z()*tracklet->t()<0.0) continue;

	  int disk=disk_;
	  if (tracklet->t()<0) disk=-disk_;
	  
	  //cout << "FPGAMatchCalculator disk="<<disk<<" "<<getName()<<endl;
	  assert(disk!=0);
	  
	  double pttracklet=0.3*3.8/(tracklet->rinv()*100);
	  bool keep=fabs(1.0/pttracklet-1.0/stub->pt())<ptstubconsistencydiskmatching;

	  keep=true;

	  if (!keep) continue;

	  //cout << "FPGAMatchCalculator past pt cut disk="<<disk<<" "<<getName()<<endl;

	  double phi=stub->phi();
	  if (phi<0) phi+=two_pi;
	  phi-=phioffset_;

	  //cout << "stub->phi() phioffset_ "<<stub->phi()<<" "<<phioffset_<<endl;
	  
	  double dz=stub->z()-tracklet->zprojdisk(disk);

	  //cout << "disk tracklet->t "<<disk<<" "<<tracklet->t()<<endl;
	  //cout << "z of stub and tracklet proj : "<<stub->z()<<" "<<tracklet->zprojdisk(disk)<<endl;
	  
	  assert(fabs(dz)<dzmax);

	  int iz=fpgastub->z().value();

	  int iphi=tracklet->fpgaphiprojdisk(disk).value();

	  int shifttmp=t2bits+tbitshift+phi0bitshift+2-rinvbitshiftdisk-phiderdiskbitshift;

	  assert(shifttmp>=0);
	  int iphicorr=(iz*tracklet->fpgaphiprojderdisk(disk).value())>>shifttmp;

	  iphi+=iphicorr;
	  
	  double phicorr=dz*tracklet->phiprojderdisk(disk);

	  assert(fabs(tracklet->phiprojderdisk(disk))<0.1);
	  assert(fabs(phicorr)<0.1);

	  double phiproj=tracklet->phiprojdisk(disk)+phicorr;
	
	  int ir=tracklet->fpgarprojdisk(disk).value()+rmindisk/krprojshiftdisk;

	  
	  int shifttmp2=rprojdiskbitshift+t3shift-rderdiskbitshift;

	  assert(shifttmp2>=0);
	  int ircorr=(iz*tracklet->fpgarprojderdisk(disk).value())>>shifttmp2;

	  int dorcorr=1;

	  ir+=ircorr*dorcorr;

	  double rcorr=dz*tracklet->rprojderdisk(disk);

	  double rproj=tracklet->rprojdisk(disk)+rcorr*dorcorr;

	  int ideltaphi=fpgastub->phi().value()*kphi/kphiproj123-iphi; 

	  double deltar=stub->r()-rproj;

	  int irstub = fpgastub->r().value();
	  if(!stub->isPSmodule()){
	    double rstub = rDSS[irstub];
	    assert (rstub < rmaxdisk);
            irstub = (1<<nrbitsdisk)*(rstub)/(rmaxdisk-rmindisk);
	  } else {
	    irstub+=rmindisk/kr;
	  }
	  
	  int ideltar=(irstub*krdisk)/krprojshiftdisk-ir;

	  //ideltar>>=1;
	  
	  //cout << "stub kr krprojshiftdisk "<<kr<<" "<<krprojshiftdisk<<endl;
	  //cout << "stub->r() : "<<stub->r()<<" "<<irstub*krdisk<<" "<<ir*krprojshiftdisk<<endl;
	  


	  double dr=stub->r()-(tracklet->rprojdisk(disk)+
			      dz*tracklet->rprojderdisk(disk)*dorcorr);

	  double dphi=phi-(tracklet->phiprojdisk(disk)+
			   dz*tracklet->phiprojderdisk(disk));

	  if (dphi>0.5*two_pi) dphi-=two_pi;
	  if (dphi<-0.5*two_pi) dphi+=two_pi;

	  if (dphi>0.5*two_pi/NSector) dphi-=two_pi/NSector;
	  if (dphi<-0.5*two_pi/NSector) dphi+=two_pi/NSector;
	  
	  //cout << "dphi "<<dphi<<" "<<phi<<" "<<tracklet->phiprojdisk(disk)
	  //     <<" "<<dz*tracklet->phiprojderdisk(disk)<<endl;
	  
	  double dphiapprox=phi-(tracklet->phiprojapproxdisk(disk)+
				 dz*tracklet->phiprojderapproxdisk(disk));
	  
	  if (dphiapprox>0.5*two_pi) dphiapprox-=two_pi;
	  if (dphiapprox<-0.5*two_pi) dphiapprox+=two_pi;

	  if (dphiapprox>0.5*two_pi/NSector) dphiapprox-=two_pi/NSector;
	  if (dphiapprox<-0.5*two_pi/NSector) dphiapprox+=two_pi/NSector;

	  
	  double drapprox=stub->r()-(tracklet->rprojapproxdisk(disk)+
				    dz*tracklet->rprojderapproxdisk(disk)*dorcorr);

	  //cout << "dr : "<<stub->r()<<" "<<drapprox<<" "<<ideltar*krprojshiftdisk<<endl;
	  
	  double alpha=0.0;
	
	  if (!stub->isPSmodule()) {
	    alpha=stub->alpha(); 	
	    dphi+=dr*alpha;
	    dphiapprox+=drapprox*alpha;
	    ideltaphi+=ideltar*fpgastub->alpha().value()*krprojshiftdisk*kalpha/kphiproj123;  
	  }	



	  double drphicut=0.20;
	  double drcut=0.75; 
	  if (!stub->isPSmodule()) {
	    drcut=3.0; //1.9
	    drphicut=0.5; 
	  }

	  bool match=(fabs(dphi)<drphicut/stub->r())&&(fabs(deltar)<drcut);

	  bool imatch=(fabs(ideltaphi)<drphicut/(kphiproj123*stub->r()))&&(fabs(ideltar)<drcut/krprojshiftdisk);

	if (debug1) {
	  cout << "imatch disk: "<<imatch<<" "<<fabs(ideltaphi)<<" "<<drphicut/(kphiproj123*stub->r())<<" "
	       <<fabs(ideltar)<<" "<<drcut/krprojshiftdisk<<" r = "<<stub->r()<<endl;
	}
	  
	if (writeDiskMatch1) {

	  static ofstream out1("diskmatch1.txt");

	  out1 << disk<<" "
	       << phiproj<<" "
	       << rproj<<" "
	       << dphi<<" "
	       << deltar<<"    "
	       << iphi*kphiprojdisk<<" "
	       << ir*krprojshiftdisk<<"  "
	       << ideltaphi*kphiprojdisk<<" "
	       << ideltar*krprojshiftdisk<<" "
	       << endl;
	  
	}
	  
	if (dumpmatch) {

	  cout << "DUMPMATCHDISK "<<tracklet->layer()<<" "<<tracklet->disk()<<" "<<stub->r()
	       <<" dr (float) = "<<deltar<<" "
	       <<" dr (int) = "<<ideltar*krprojshiftdisk
	       <<" dr res = "<<krprojshiftdisk
	       <<" rcorr = "<<rcorr
	       <<" ircorr = "<<ircorr*krprojshiftdisk
	       <<endl;

		 cout << "DUMPMATCHDISK rproj(float) = "<<tracklet->rprojdisk(disk)
		      <<" rproj(int) = "<<(tracklet->fpgarprojdisk(disk).value()+rmindisk/krprojshiftdisk)*krprojshiftdisk
		      <<endl;
		 
	  cout << "DUMPMATCHDISK1 : "<<disk<<" "
	       <<tracklet->phiprojdisk(disk)<<" "
	       <<dz<<" "
	    //<<tracklet->phiprojderdisk(disk)<<" "
	       <<phicorr<<" "
	       <<phiproj<<" | "
	       <<tracklet->rprojdisk(disk)<<" "
	       <<tracklet->rprojderdisk(disk)<<" "
	       <<rcorr<<" "
	       <<rproj<<" dphi= "
	       <<dphi<<" dr= "
	       <<deltar<<" "
	       <<match
	       <<endl;
	  cout << "DUMPMATCHDISK2 : "<<disk<<" "
	       <<tracklet->fpgaphiprojdisk(disk).value()*kphiproj123<<" "
	       <<iz*kzdisk<<" "
	    //<<tracklet->fpgaphiprojderdisk(disk).value()*kphiprojderdiskshift<<" "
	       <<iphicorr*kphiproj123<<" ("<<iphicorr<<") "
	       <<iphi*kphiproj123<<" | "
	       <<tracklet->fpgarprojdisk(disk).value()*krprojshiftdisk<<" "
	       <<tracklet->fpgarprojderdisk(disk).value()*krprojderdiskshift<<" "
	       <<ircorr*krprojshiftdisk<<" "
	       <<ir*krprojshiftdisk<<" dphi= "
	       <<ideltaphi*kphiproj123<<" dr= "
	       <<ideltar*krprojshiftdisk<<" "
	       <<imatch
	       <<endl;
	}

	  
	  if (imatch) {

	    //cout << "2: ideltaphi = "<<ideltaphi<<endl;

	    std::pair<FPGAStub*,L1TStub*> tmp(fpgastub,stub);

	    countsel++;

	    if (debug1) {
	      cout << "FPGAMatchCalculator found match in disk "<<getName()<<endl;
	    }

	    assert(fabs(dphi)<0.2);
	    assert(fabs(dphiapprox)<0.2);
	    
	    tracklet->addMatchDisk(disk,ideltaphi,ideltar,
				   dphi,dr,dphiapprox,drapprox,
				   stub->alphatruncated(),
				   fpgastub->stubindex().value()+((fpgastub->fedregion()-1)<<6),
				   stub->z(),tmp);

	    
	    if (tracklet->plusNeighborDisk(disk)){
	      fullmatchesToMinus_->addMatch(tracklet,tmp);
	      if (debug1) {
		cout << "Accepted full match to minus in disk " <<getName()<<" "<<tracklet
		     <<" "<<fullmatchesToMinus_->getName()<<endl;
	      }
		  int nmatch=fullmatchesToMinus_->nMatches();
		  if (nmatch>1) {
		    assert(fullmatchesToMinus_->getFPGATracklet(nmatch-2)->TCID()<
				   fullmatchesToMinus_->getFPGATracklet(nmatch-1)->TCID());
		  }
	    } else if (tracklet->minusNeighborDisk(disk)) {
	      fullmatchesToPlus_->addMatch(tracklet,tmp);
	      if (debug1) {
		cout << "Accepted full match to plus in disk " <<getName()<<" "<<tracklet
		     <<" "<<fullmatchesToPlus_->getName()<<endl;
	      }
	    } else {
	      //cout << "Found new match:"<<endl;
	      for (unsigned int l=0;l<fullmatches_.size();l++){
		//cout << "Try to add full match in disk " <<getName()<<" layer="<<tracklet->layer()<<" disk="<<tracklet->disk()
		//     <<" "<<fullmatches_[l]->getName()<<endl;
		if (((abs(tracklet->disk())==1&&tracklet->layer()==1)&&fullmatches_[l]->getName().substr(3,4)=="D1L1")||
		    (tracklet->layer()==2&&fullmatches_[l]->getName().substr(3,4)=="D1L2")||    //dangerous to check only layer!!!
		    ((abs(tracklet->disk())==1&&tracklet->layer()==0)&&fullmatches_[l]->getName().substr(3,4)=="D1D2")||
		    ((tracklet->disk()==0&&tracklet->layer()==1)&&fullmatches_[l]->getName().substr(3,4)=="L1L2")||
                    ((tracklet->disk()==0&&tracklet->layer()==3)&&fullmatches_[l]->getName().substr(3,4)=="L3L4")||
		    ((abs(tracklet->disk())==3&&tracklet->layer()==0)&&fullmatches_[l]->getName().substr(3,4)=="D3D4")){
		  fullmatches_[l]->addMatch(tracklet,tmp);
		  if (debug1) {
		    cout << "In "<<getName()<<" added match to "<<fullmatches_[l]->getName()<<endl;
		  }
		}
	      }
	    }
	  }
	}
 	//if (countall>=MAXMC) break;
    //  }
      if (countall>=MAXMC) break;
    }


    if (writeMatchCalculator) {
      static ofstream out("matchcalculator.txt");
      out << getName()<<" "<<countall<<" "<<countsel<<endl;
    }


  }
    
private:

  int layer_;
  int disk_;
  double phimin_;
  double phimax_;
  double phioffset_;

  FPGAAllStubs* allstubs_;
  FPGAAllProjections* allprojs_;

  vector<FPGACandidateMatch*> matches_;

  vector<FPGAFullMatch*> fullmatches_;
  FPGAFullMatch* fullmatchesToPlus_;
  FPGAFullMatch* fullmatchesToMinus_;

};

#endif
