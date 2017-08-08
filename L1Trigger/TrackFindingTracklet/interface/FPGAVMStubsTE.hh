//This class holds the reduced VM stubs
#ifndef FPGAVMSTUBSTE_H
#define FPGAVMSTUBSTE_H

#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGATEBin.hh"
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
    
  }
  
  void addStub(std::pair<FPGAStub*,L1TStub*> stub) {
    //cout << "addStub "<<getName()<<endl;
    stubs_.push_back(stub);
    if (stub.first->isBarrel()){
      FPGAWord z=stub.first->z();
      int bin = ((z.value()+(1<<(z.nbits()-1)))>>(z.nbits()-NLONGVMBITS));
      assert(bin>=0);
      assert((unsigned int)bin<NLONGVMBINS);
      stubsbinned_[bin].push_back(stub);
    } else {
      int ir=stub.first->r().value();
      assert(NLONGVMRBITS==3); //This code is not generic...
      int ir1=(rinnerdisk+(routerPSdisk-rinnerdisk)*0.25)/kr;
      int ir2=(rinnerdisk+(routerPSdisk-rinnerdisk)*0.50)/kr;
      int ir3=(rinnerdisk+(routerPSdisk-rinnerdisk)*0.75)/kr;
      int offset=0;
      if (stub.first->disk().value()<0) offset=4;
      int bin=3+offset;
      if (ir<ir3) bin=2+offset;
      if (ir<ir2) bin=1+offset;
      if (ir<ir1) bin=0+offset;
      //cout << "Adding stub in bin " << bin << endl;
      stubsbinned_[bin].push_back(stub);
    }
  }

  unsigned int nStubs() const {return stubs_.size();}
  unsigned int nStubsBinned(unsigned int bin) const {return stubsbinned_[bin].size();}

  FPGAStub* getFPGAStub(unsigned int i) const {return stubs_[i].first;}
  L1TStub* getL1TStub(unsigned int i) const {return stubs_[i].second;}
  std::pair<FPGAStub*,L1TStub*> getStub(unsigned int i) const {return stubs_[i];}

  FPGAStub* getFPGAStubBinned(unsigned int bin, unsigned int i) const {return stubsbinned_[bin][i].first;}
  L1TStub* getL1TStubBinned(unsigned int bin, unsigned int i) const {return stubsbinned_[bin][i].second;}
  std::pair<FPGAStub*,L1TStub*> getStubBinned(unsigned int bin,unsigned int i) const {return stubsbinned_[bin][i];}

  int getzbin(unsigned int i,int overlap=false) {
    FPGAStub* stub=stubs_[i].first;
    //cout << "r of stub "<<stubs_[i].second->r()<<endl;
    //cout << "z of stub "<<stubs_[i].second->z()<<endl;
    int bin = ((stub->z().value()+(1<<(stub->z().nbits()-1)))>>(stub->z().nbits()-NLONGVMODDLAYERBITS));
    int layer=stub->layer().value()+1;
    assert (layer%2==1); //can only call this on layer 1, 3, and 5. Remember that FPGA stub has layer
                                        //counting from zero

    if (layer==1) {
      static FPGATEBin *layer1=0;
      if (layer1==0) {
	layer1= new FPGATEBin;
	layer1->initLayerLayer(1,false);
	if (writeTEBinTable) layer1->writeTable("TEBinTableLayer1ToLayer2.txt");
      }
      return layer1->lookup(bin);
    }

    if (layer==3&&(!overlap)) {
      static FPGATEBin *layer3=0;
      if (layer3==0) {
	layer3= new FPGATEBin;
	layer3->initLayerLayer(3,false);
	if (writeTEBinTable) layer3->writeTable("TEBinTableLayer3ToLayer4.txt");

      }
      return layer3->lookup(bin);
    }

    if (layer==3&&overlap) {
      static FPGATEBin *layer3overlap=0;
      if (layer3overlap==0) {
	layer3overlap= new FPGATEBin;
	layer3overlap->initLayerLayer(3,true);
	if (writeTEBinTable) layer3overlap->writeTable("TEBinTableLayer3ToLayer2.txt");
      }
      return layer3overlap->lookup(bin);
    }

    if (layer==5) {
      static FPGATEBin *layer5=0;
      if (layer5==0) {
	layer5= new FPGATEBin;
	layer5->initLayerLayer(5,false);
	if (writeTEBinTable) layer5->writeTable("TEBinTableLayer5ToLayer6.txt");

      }
      return layer5->lookup(bin);
    }

    assert(0);
    return 0;
    
  }

  int getrbin(unsigned int i) {
    FPGAStub* stub=stubs_[i].first;
    int disk=abs(stub->disk().value());
    assert (disk%2==1); //can only call this on disk 1, 3.
    int bin = (stub->r().value())>>(stub->r().nbits()-NLONGVMODDDISKBITS+1);
    assert(bin>=0);
    int ibin=bin;
    if (stub->disk().value()<0) {
      ibin+=(1<<(NLONGVMODDDISKBITS-1));
    }
    assert((unsigned int)bin<(1<<NLONGVMODDDISKBITS));

    if (disk==1) {
      static FPGATEBin *disk1=0;
      if (disk1==0) {
	disk1= new FPGATEBin;
	disk1->initDiskDisk(1);
	if (writeTEBinTable) disk1->writeTable("TEBinTableDisk1ToDisk2.txt");
      }
      return disk1->lookup(ibin);
    }

    if (disk==3) {
      static FPGATEBin *disk3=0;
      if (disk3==0) {
	disk3= new FPGATEBin;
	disk3->initDiskDisk(3);
	if (writeTEBinTable) disk3->writeTable("TEBinTableDisk3ToDisk4.txt");

      }
      return disk3->lookup(ibin);
    }

    assert(0);
    return 0;

    
    
  }


  int getrbinfromlayer(unsigned int i) {
    FPGAStub* stub=stubs_[i].first;
    assert(stub->isBarrel());
    int layer=stub->layer().value()+1;
    int bin = ((stub->z().value()+(1<<(stub->z().nbits()-1)))>>(stub->z().nbits()-NLONGVMODDLAYERBITS));

    assert(layer==1);
    
    if (layer==1) {
      static FPGATEBin *layer1overlap=0;
      if (layer1overlap==0) {
	layer1overlap= new FPGATEBin;
	layer1overlap->initLayertoDisk(1);
	if (writeTEBinTable) layer1overlap->writeTable("TEBinTableLayer1ToDisk1.txt");

      }
      return layer1overlap->lookup(bin);
    }

    assert(0);
    return 0;
    
  }

  int getzbinfromdisk(unsigned int i, unsigned int layer=1) {
    FPGAStub* stub=stubs_[i].first;
    assert(stub->isDisk());
    int disk=abs(stub->disk().value());
    assert (disk==1); //can only call this on disk 1

    //cout << "r value and nbits : "<<stub->r().value()<<" "<<stub->r().nbits()<<endl;

    int bin = (stub->r().value())>>(stub->r().nbits()-NLONGVMODDDISKBITS+1);
    assert(bin>=0);

    int ibin=bin;
    if (stub->disk().value()<0) {
      ibin+=(1<<(NLONGVMODDDISKBITS-1));
    }
    assert((unsigned int)bin<(1<<NLONGVMODDDISKBITS));

    if (layer==1) {
      static FPGATEBin *layer1=0;
      if (layer1==0) {
	layer1= new FPGATEBin;
	layer1->initDisktoLayer(1);
	if (writeTEBinTable) layer1->writeTable("TEBinTableDisk1ToLayer1.txt");
      }
      return layer1->lookup(ibin);
    }

    if (layer==2) {
      static FPGATEBin *layer2=0;
      if (layer2==0) {
	layer2= new FPGATEBin;
	layer2->initDisktoLayer(2);
	if (writeTEBinTable) layer2->writeTable("TEBinTableDisk1ToLayer2.txt");
      }
      return layer2->lookup(ibin);
    }

    assert(0);
    return 0;
    
  }

  
  void clean() {
    stubs_.clear();
    for (unsigned int i=0;i<NLONGVMBINS;i++){
      stubsbinned_[i].clear();
    }
  }

  void writeStubs(bool first) {

    std::string fname="VMStubs_";
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
	FPGAWord tmp;
	tmp.set(getzbin(j),4,true,__LINE__,__FILE__);
	stub+=tmp.str();  
	if (j<16) out_ <<"0";
	out_ << hex << j << dec ;
	out_ <<" "<<stub<< endl;
      }
    } else if (layer_!=0) {
      for (unsigned int i=0;i<NLONGVMBINS;i++) {
	for (unsigned int j=0;j<stubsbinned_[i].size();j++){
	  string stub=stubsbinned_[i][j].first->stubindex().str();
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
	  FPGAWord tmp;
	  tmp.set(getrbin(j),4,true,__LINE__,__FILE__);
	  stub+=tmp.str();  
	  if (j<16) out_ <<"0";
	  out_ << hex << j << dec ;
	  out_ <<" "<<stub<< endl;
	}
      } else if (disk_!=0) {
	for (unsigned int i=0;i<NLONGVMBINS;i++) {
	  for (unsigned int j=0;j<stubsbinned_[i].size();j++){
	    string stub=stubsbinned_[i][j].first->stubindex().str();
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



private:

  int layer_;
  int disk_;
  bool overlap_;
  double phimin_;
  double phimax_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubs_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubsbinned_[NLONGVMBINS];

};

#endif
