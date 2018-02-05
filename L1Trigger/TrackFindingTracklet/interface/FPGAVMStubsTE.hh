//This class holds the reduced VM stubs
#ifndef FPGAVMSTUBSTE_H
#define FPGAVMSTUBSTE_H

#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGATETableOuter.hh"
#include "FPGATETableInner.hh"
#include "FPGATETableOuterDisk.hh"
#include "FPGATETableInnerDisk.hh"
#include "FPGATETableInnerOverlap.hh"
#include "FPGAMemoryBase.hh"

using namespace std;

class FPGAVMStubsTE:public FPGAMemoryBase{

public:

  FPGAVMStubsTE(string name, unsigned int iSector, 
	      double phimin, double phimax):
    FPGAMemoryBase(name,iSector){
    phimin_=phimin;
    phimax_=phimax;
    string subname=name.substr(6,2);
    layer_ = 0;
    disk_  = 0;
    phibin_=-1;
    other_=0;
    if (subname=="L1") layer_=1;
    if (subname=="L2") layer_=2;
    if (subname=="L3") layer_=3;
    if (subname=="L4") layer_=4;
    if (subname=="L5") layer_=5;
    if (subname=="L6") layer_=6;
    if (subname=="D1") disk_=1;
    if (subname=="D2") disk_=2;
    if (subname=="D3") disk_=3;
    if (subname=="D4") disk_=4;
    if (subname=="D5") disk_=5;
    if (layer_==0&&disk_==0) {
      cout << name<<" subname = "<<subname<<" "<<layer_<<" "<<disk_<<endl;
    }
    assert((layer_!=0)||(disk_!=0));
    overlap_=false;
    subname=name.substr(11,1);
    if (subname=="X") overlap_=true;
    if (subname=="Y") overlap_=true;
    if (subname=="W") overlap_=true;
    if (subname=="Q") overlap_=true;

    subname=name.substr(12,2);
    phibin_=subname[0]-'0';
    if (subname[1]!='n') {
      phibin_=10*phibin_+(subname[1]-'0');
    }
    
  }
  
  bool addStub(std::pair<FPGAStub*,L1TStub*> stub) {

    static bool first=true;
    static FPGATETableInner innerTableL1_;
    static FPGATETableInner innerTableL3_;
    static FPGATETableInner innerTableL5_;
    static FPGATETableOuter outerTableL2_;
    static FPGATETableOuter outerTableL4_;
    static FPGATETableOuter outerTableL6_;

    static FPGATETableInnerDisk innerTableD1_;
    static FPGATETableInnerDisk innerTableD3_;
    static FPGATETableOuterDisk outerTableD2_;
    static FPGATETableOuterDisk outerTableD4_;

    static FPGATETableInnerOverlap innerTableOverlapL1_;
    static FPGATETableInnerOverlap innerTableOverlapL2_;
    static FPGATETableOuterDisk outerTableOverlapD1_;

    
    if (first) {
      //layer-layer seeding
      innerTableL1_.init(1,2,7,4);
      innerTableL3_.init(3,4,7,4);
      innerTableL5_.init(5,6,7,4);
      outerTableL2_.init(2,7,4);
      outerTableL4_.init(4,7,4);
      outerTableL6_.init(6,7,4);

      //disk-disk seeding
      innerTableD1_.init(1,2,7,3);
      innerTableD3_.init(3,4,7,3);
      outerTableD2_.init(2,7,3);
      outerTableD4_.init(4,7,3);

      //overlap seeding
      innerTableOverlapL1_.init(1,1,7,3);
      innerTableOverlapL2_.init(2,1,7,3);
      outerTableOverlapD1_.init(1,7,3);
      
      first=false;
    }
   

    if (overlap_) {
      if (stub.first->isBarrel()){
	if (debug1) cout << "Have barrel stub for overlap in layer = "<<layer_<<endl;
	FPGAWord z=stub.first->z();
	FPGAWord r=stub.first->r();
	int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-7);
	int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-3);
	int binlookup=-1;
	switch (layer_) {
	  case 1 : binlookup=innerTableOverlapL1_.lookup(zbin,rbin);
	    break;
	  case 2 : binlookup=innerTableOverlapL2_.lookup(zbin,rbin);
	    break;
	  case 0 : return false; //FIXME - why are we here???
	    break;
	default: assert(0);
	}
	if (binlookup<0) {
	  if (debug1) cout << getName() << " binlookup<0 returning"<<endl;
	  return false;
	}
	if (debug1) cout << getName()<<" Stub with lookup = "<<binlookup<<" in layer = "<<layer_<<endl;
	stub.first->setVMBits(binlookup);
      } else {
	if (debug1) cout << "Have disk stub for overlap in disk = "<<disk_<<endl;
	FPGAWord r=stub.first->r();
	FPGAWord z=stub.first->z();
	int rbin=(r.value())>>(r.nbits()-7);
	int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-3);
	bool negdisk=stub.first->disk().value()<0.0;
	if (negdisk) zbin=7-zbin; //Should this be separate table?
	int binlookup=-1;
	if (disk_==1) {
	  binlookup=outerTableOverlapD1_.lookup(rbin,zbin);
	  assert(binlookup>=0);
	  stub.first->setVMBits(binlookup);
	  int bin=binlookup/8;
	  assert(bin<4);
	  if (negdisk) bin+=4;
	  stubsbinned_[bin].push_back(stub);
	  if (debug1) cout << getName()<<" Stub with lookup = "<<binlookup
			   <<" in disk = "<<disk_<<"  in bin = "<<bin<<endl;
	}
      }
    } else {
      if (stub.first->isBarrel()){
	FPGAWord z=stub.first->z();
	FPGAWord r=stub.first->r();
	int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-7);
	int rbin=(r.value()+(1<<(r.nbits()-1)))>>(r.nbits()-4);
	int binlookup=-1;
	switch (layer_) {
	  case 2 : binlookup=outerTableL2_.lookup(zbin,rbin);
	    break;
	  case 4 : binlookup=outerTableL4_.lookup(zbin,rbin);
	    break;
	  case 6 : binlookup=outerTableL6_.lookup(zbin,rbin);
	    break;
	  case 1 : binlookup=innerTableL1_.lookup(zbin,rbin);
	    break;
	  case 3 : binlookup=innerTableL3_.lookup(zbin,rbin);
	    break;
	  case 5 : binlookup=innerTableL5_.lookup(zbin,rbin);
	    break;
       	  default : assert(0);
	}

	if (layer_%2==0) {
	  assert(binlookup>=0);
	  int bin=binlookup/8;
	  stubsbinned_[bin].push_back(stub);
	} else {
	  if (binlookup<0) return false;
	}
	stub.first->setVMBits(binlookup);

	
      } else {
	FPGAWord r=stub.first->r();
	FPGAWord z=stub.first->z();
	int rbin=(r.value())>>(r.nbits()-7);
	int zbin=(z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-3);
	bool negdisk=stub.first->disk().value()<0.0;
	if (negdisk) zbin=7-zbin; //Should this be separate table?
	int binlookup=-1;

	switch (disk_) {
	  case 2 : binlookup=outerTableD2_.lookup(rbin,zbin);
	    break;
	  case 4 : binlookup=outerTableD4_.lookup(rbin,zbin);
	    break;
	  case 1 : binlookup=innerTableD1_.lookup(rbin,zbin);
	    break;
	  case 3 : binlookup=innerTableD3_.lookup(rbin,zbin);
	    break;
	  default : assert(0);  
	}

	if (disk_%2==0) {
	  assert(binlookup>=0);
	  int bin=binlookup/8;
	  assert(bin<4);
	  if (negdisk) bin+=4;
	  stubsbinned_[bin].push_back(stub);
	} else {
	  if (binlookup<0) return false;
	}
	stub.first->setVMBits(binlookup);
	
      }
    }
    stubs_.push_back(stub);
    return true;
  }

  unsigned int nStubs() const {return stubs_.size();}
  unsigned int nStubsBinned(unsigned int bin) const {return stubsbinned_[bin].size();}

  FPGAStub* getFPGAStub(unsigned int i) const {return stubs_[i].first;}
  L1TStub* getL1TStub(unsigned int i) const {return stubs_[i].second;}
  std::pair<FPGAStub*,L1TStub*> getStub(unsigned int i) const {return stubs_[i];}

  FPGAStub* getFPGAStubBinned(unsigned int bin, unsigned int i) const {return stubsbinned_[bin][i].first;}
  L1TStub* getL1TStubBinned(unsigned int bin, unsigned int i) const {return stubsbinned_[bin][i].second;}
  std::pair<FPGAStub*,L1TStub*> getStubBinned(unsigned int bin,unsigned int i) const {return stubsbinned_[bin][i];}


  
  void clean() {
    stubs_.clear();
    for (unsigned int i=0;i<NLONGVMBINS;i++){
      stubsbinned_[i].clear();
    }
  }

  void writeStubs(bool first) {

    std::string fname="MemPrints/VMStubsTE/VMStubs_";
    fname+=getName();
    //get rid of duplicates
    int len = fname.size();
    if(fname[len-2]=='n'&& fname[len-1]>'1'&&fname[len-1]<='9') return;
    //
    fname+="_";
    ostringstream oss;
    oss << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=oss.str();
    fname+=".dat";
    if (first) {
      bx_ = 0;
      event_ = 1;
      out_.open(fname.c_str());
    }
    else
      out_.open(fname.c_str(),std::ofstream::app);

    out_ << "BX = "<<(bitset<3>)bx_ << " Event : " << event_ << endl;

    if (layer_%2==1) {
      for (unsigned int j=0;j<stubs_.size();j++){
	string stub=stubs_[j].first->stubindex().str();
	stub+="|";
    stub+=stubs_[j].first->stubaddressaste().str();
    stub+="|";
	FPGAWord tmp;
	tmp.set(stubs_[j].first->getVMBits().value(),10,true,__LINE__,__FILE__);
	stub+=tmp.str();  
	if (j<16) out_ <<"0";
	out_ << hex << j << dec ;
	out_ <<" "<<stub<< endl;
      }
    } else if (layer_!=0) {
      for (unsigned int i=0;i<NLONGVMBINS;i++) {
	for (unsigned int j=0;j<stubsbinned_[i].size();j++){
	  string stub=stubsbinned_[i][j].first->stubindex().str();
      stub+="|";
      stub+=stubsbinned_[i][j].first->stubaddressaste().str();
	  out_ << hex << i << " " << j << dec << " "<<stub<<endl;
	}
      }
    }

    //Need to be done smarter to handle overlap region
    if (!overlap_) {
      if (disk_%2==1) {
	for (unsigned int j=0;j<stubs_.size();j++){
	  string stub=stubs_[j].first->stubindex().str();
	  stub+="|";
      stub+=stubs_[j].first->stubaddressaste().str();
      stub+="|";
	  FPGAWord tmp;
	  tmp.set(stubs_[j].first->getVMBits().value(),10,true,__LINE__,__FILE__);
	  stub+=tmp.str();  
	  if (j<16) out_ <<"0";
	  out_ << hex << j << dec ;
	  out_ <<" "<<stub<< endl;
	}
      } else if (disk_!=0) {
	for (unsigned int i=0;i<NLONGVMBINS;i++) {
	  for (unsigned int j=0;j<stubsbinned_[i].size();j++){
	    string stub=stubsbinned_[i][j].first->stubindex().str();
        stub+="|";
        stub+=stubsbinned_[i][j].first->stubaddressaste().str();
	    out_ << hex << i << " " << j << dec << " "<<stub<<endl;
	  }
	}
      }
    }
    

    out_.close();

    bx_++;
    event_++;
    if (bx_>7) bx_=0;

  }

  int phibin() const {
    return phibin_;
  }

  void getPhiRange(double &phimin, double &phimax) {

    if (layer_==1 || layer_==3 || layer_==5) {
      int nphibin=24;
      double dphi=dphisector/nphibin;
      phimax=phibin()*dphi;
      phimin=phimax-dphi;
      return;
    }

    if (layer_==2 || layer_==4 || layer_==6) {
      int nphibin=12;
      double dphi=dphisector/nphibin;
      phimax=phibin()*dphi-dphisector/6.0;
      phimin=phimax-dphi;
      return;
    }

    assert(0); //not implemented yet
    
  }
  
  void setother(FPGAVMStubsTE* other){
    other_=other;
  }
  
  FPGAVMStubsTE* other() const {
    return other_;
  }

  

private:

  int layer_;
  int disk_;
  int phibin_;
  FPGAVMStubsTE* other_;
  bool overlap_;
  double phimin_;
  double phimax_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubs_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubsbinned_[NLONGVMBINS];

};

#endif
