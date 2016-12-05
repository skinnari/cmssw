//Holds the candidate matches
#ifndef FPGACANDIDATEMATCH_H
#define FPGACANDIDATEMATCH_H

#include "FPGATracklet.hh"
#include "FPGAMemoryBase.hh"
#include "FPGAStub.hh"
#include "L1TStub.hh"

using namespace std;

class FPGACandidateMatch:public FPGAMemoryBase{

public:

  FPGACandidateMatch(string name, unsigned int iSector, 
		     double phimin, double phimax):
    FPGAMemoryBase(name,iSector){
    phimin_=phimin;
    phimax_=phimax;
  }

  void addMatch(FPGATracklet* tracklet,std::pair<FPGAStub*,L1TStub*> stub) {
    std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > tmp(tracklet,stub);
    matches_.push_back(tmp);
  }

  unsigned int nMatches() const {return matches_.size();}

  FPGATracklet* getFPGATracklet(unsigned int i) const {return matches_[i].first;}
  std::pair<FPGAStub*,L1TStub*> getStub(unsigned int i) const {return matches_[i].second;}
  std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > getMatch(unsigned int i) const { return matches_[i];}
  
  void clean() {
    matches_.clear();
  }

  void writeCandidate(bool first) {

    string subname=getName().substr(8,2);
    int layer=0;
    int disk=0;

    if (subname=="L1") layer=1;
    if (subname=="L2") layer=2;
    if (subname=="L3") layer=3;
    if (subname=="L4") layer=4;
    if (subname=="L5") layer=5;
    if (subname=="L6") layer=6;

    
    if (subname=="F1") disk=1;
    if (subname=="F2") disk=2;
    if (subname=="F3") disk=3;
    if (subname=="F4") disk=4;
    if (subname=="F5") disk=5;

    
    assert(layer!=0||disk!=0);


    std::string fname="./MemPrints/Matches/CandidateMatches_";
    fname+=getName();
    fname+="_";
    ostringstream oss;
    oss << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=oss.str();
    fname+=".dat";
    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else{
      out_.open(fname.c_str(),std::ofstream::app);
    
    }
    out_ << "BX = "<<(bitset<3>)bx_ << " Event : " << event_ << endl;
    for (unsigned int j=0;j<matches_.size();j++){
      if (j<16) out_ <<"0";
      out_ << hex << j << dec ;
      if (layer!=0){
	out_ <<" "<<matches_[j].first->trackletprojstrlayer(layer) <<" "<<matches_[j].second.first->strbare()<< endl;
      }
      if (disk!=0) {
	out_ <<" "<<matches_[j].first->trackletprojstrdisk(disk) <<" "<<matches_[j].second.first->strbare()<< endl;
      }
    }
    out_.close();
    bx_++;
    event_++;
    if (bx_>7) bx_=0;

  }
private:

  double phimin_;
  double phimax_;
  std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > matches_;

};

#endif
