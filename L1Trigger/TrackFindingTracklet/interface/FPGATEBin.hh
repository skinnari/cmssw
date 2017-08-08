#ifndef FPGATEBIN_H
#define FPGATEBIN_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>


using namespace std;

class FPGATEBin{

public:

  FPGATEBin() {
   
  }

  ~FPGATEBin() {

  }


  void initLayerLayer(int layer, bool overlap=false){

    for(unsigned int bin=0;bin<(1<<NLONGVMODDLAYERBITS);bin++) { 
      double zmin=bin*2*zlength/(1<<NLONGVMODDLAYERBITS)-zlength;
      double zmax=(bin+1)*2*zlength/(1<<NLONGVMODDLAYERBITS)-zlength;
      //cout << "zmin zmax : "<<zmin<<" "<<zmax<<endl;
      double z1=z0cut+rmean[layer]*(zmin-z0cut)/rmean[layer-1];
      double z2=-z0cut+rmean[layer]*(zmax+z0cut)/rmean[layer-1];
      if (overlap) {
	z1=-z0cut+rmean[layer-2]*(zmin+z0cut)/rmean[layer-1];
	z2=z0cut+rmean[layer-2]*(zmax-z0cut)/rmean[layer-1];
	assert(z1<z2);
      }
      //cout << "z0cut "<<z0cut<<endl;
      //cout << "z1 z2 : "<<z1<<" "<<z2<<endl;
      if (z1<-zlength) z1=-zlength+1.0;
      if (z1>zlength) z1=zlength-1.0;
      if (z2<-zlength) z2=-zlength+1.0;
      if (z2>zlength) z2=zlength-1.0;
      int iz1=z1/kz;
      int iz2=z2/kz;
      int bin1 = ((iz1+(1<<(nbitszL123-1)))>>(nbitszL123-NLONGVMBITS));
      assert(bin1>=0);
      assert((unsigned int)bin1<NLONGVMBINS);
      int bin2 = ((iz2+(1<<(nbitszL123-1)))>>(nbitszL123-NLONGVMBITS));
      //cout << "bin1 bin2 : "<<bin1<<" "<<bin2<<endl;
      assert(bin2>=0);
      assert((unsigned int)bin2<NLONGVMBINS);
      assert(bin1<=bin2);
      if (bin1==bin2){
	table_.push_back(bin1<<1);
      } else {
	assert(bin1+1==bin2);
	table_.push_back((bin1<<1)+1);
      }
    }

  }
	    
  void initDiskDisk(int disk){

    for(unsigned int bin=0;bin<(1<<NLONGVMODDDISKBITS);bin++) { 

      //cout << "bin "<<bin<<endl;
      
      unsigned int ibin=bin&((1<<(NLONGVMODDDISKBITS-1))-1);
      double rmin=ibin*rmaxdisk/(1<<(NLONGVMODDDISKBITS-1));
      double rmax=(ibin+1)*rmaxdisk/(1<<(NLONGVMODDDISKBITS-1));
      //cout << "rmin rmax : "<<rmin<<" "<<rmax<<endl;
      double r1=rmin*(zmean[disk]+z0cut)/(zmean[disk-1]+z0cut);
      double r2=rmax*(zmean[disk]-z0cut)/(zmean[disk-1]-z0cut);
      //cout << "bin ibin r1 r2 : "<<bin<<" "<<ibin<<" "<<r1<<" "<<r2<<endl;
      assert(r1>=0.0);
      assert(r2>r1);
      if (r1>routerPSdisk) r1=routerPSdisk-1.0;
      if (r2>routerPSdisk) r2=routerPSdisk-1.0;
      int bin1 = (r1-rinnerdisk)*(1<<(NLONGVMRBITS-1))/(routerPSdisk-rinnerdisk);
      //cout << "bin1 : "<<bin1<<endl;
      if (bin1<0) bin1=0;
      assert((unsigned int)bin1<NLONGVMRBINS/2);
      int bin2 = (r2-rinnerdisk)*(1<<(NLONGVMRBITS-1))/(routerPSdisk-rinnerdisk);
      if (bin2<0) bin2=0;
      assert((unsigned int)bin2<NLONGVMRBINS/2);
      assert(bin1<=bin2);
      if (bin!=ibin) {
	bin1+=NLONGVMRBINS/2;
	bin2+=NLONGVMRBINS/2;
      }
      assert(bin1<=bin2);
      if (bin1==bin2){
	table_.push_back(bin1<<1);
      } else {
	assert(bin1+1==bin2);
	table_.push_back((bin1<<1)+1);
      }
    }

  }
	    

  void initLayertoDisk(int layer){

    for(unsigned int bin=0;bin<(1<<NLONGVMODDLAYERBITS);bin++) { 
      double zmin=bin*2*zlength/(1<<NLONGVMODDLAYERBITS)-zlength;
      double zmax=(bin+1)*2*zlength/(1<<NLONGVMODDLAYERBITS)-zlength;
      //cout << "zmin zmax : "<<" "<<zmin<<" "<<zmax<<endl;
      //cout << "rmean layer : "<<rmean[layer-1]<<" "<<layer<<endl;
      double zm=zmean[0];
      if (zmin<0.0) zm=-zm;
      double r1=rmean[layer-1]*(zm-z0cut)/(zmin-z0cut);
      double r2=rmean[layer-1]*(zm+z0cut)/(zmax+z0cut);
      //cout << "r1 r2 : "<<r1<<" "<<r2<<endl;
      assert(r1>rinnerdisk);
      assert(r2>rinnerdisk);
      if (r1>routerPSdisk) r1=routerPSdisk-1.0;
      if (r2>routerPSdisk) r2=routerPSdisk-1.0;
      int bin1 = (r1-rinnerdisk)*(1<<(NLONGVMRBITS-1))/(routerPSdisk-rinnerdisk);
      assert(bin1>=0);
      assert((unsigned int)bin1<NLONGVMRBINS/2);
      int bin2 = (r2-rinnerdisk)*(1<<(NLONGVMRBITS-1))/(routerPSdisk-rinnerdisk);
      //cout << "bin1 bin2 : "<<bin1<<" "<<bin2<<endl;
      assert(bin2>=0);
      assert((unsigned int)bin2<NLONGVMRBINS/2);
      assert(bin1<=bin2);
      if (zmin<0.0) {
	bin1+=NLONGVMRBINS/2;
	bin2+=NLONGVMRBINS/2;
      }
      if (bin1==bin2){
	table_.push_back(bin1<<1);
      } else {
	table_.push_back((bin1<<1)+1);
      }
    }
  }


    void initDisktoLayer(int layer){

    for(unsigned int bin=0;bin<(1<<NLONGVMODDDISKBITS);bin++) { 

      //cout << "bin "<<bin<<endl;
      
      unsigned int ibin=bin&((1<<(NLONGVMODDDISKBITS-1))-1);

      assert((unsigned int)ibin<(1<<(NLONGVMODDDISKBITS-1)));
      double rmin=ibin*rmaxdisk/(1<<(NLONGVMODDDISKBITS-1));
      double rmax=(ibin+1)*rmaxdisk/(1<<(NLONGVMODDDISKBITS-1));
      //cout << "rmin rmax : "<<rmin<<" "<<rmax<<endl;
      double z1=z0cut+rmax*(zmean[0]-z0cut)/rmean[layer-1];
      double z2=-z0cut+rmin*(zmean[0]+z0cut)/rmean[layer-1];
      //cout << "z0cut "<<z0cut<<endl;
      if (rmin<rmean[layer-1]) {
	z1=zlength-1.0;
	z2=zlength-1.0;
      }
      if (ibin!=bin) {
	z1=-z1;
	z2=-z2;
      }
      //cout << "z1 z2 : "<<z1<<" "<<z2<<endl;
      if (z1<-zlength) z1=-zlength+1.0;
      if (z1>zlength) z1=zlength-1.0;
      if (z2<-zlength) z2=-zlength+1.0;
      if (z2>zlength) z2=zlength-1.0;
      int iz1=z1/kz;
      int iz2=z2/kz;
      int bin1 = ((iz1+(1<<(nbitszL123-1)))>>(nbitszL123-NLONGVMBITS));
      assert(bin1>=0);
      assert((unsigned int)bin1<NLONGVMBINS);
      int bin2 = ((iz2+(1<<(nbitszL123-1)))>>(nbitszL123-NLONGVMBITS));
      //cout << "bin1 bin2 : "<<bin1<<" "<<bin2<<endl;
      assert(bin2>=0);
      assert((unsigned int)bin2<NLONGVMBINS);
    
      assert(bin1<=bin2);
      if (bin1==bin2){
	table_.push_back(bin1<<1);
      } else {
	assert(bin1+1==bin2);
	table_.push_back((bin1<<1)+1);
      }
    }

  }
	    

  
  
  void writeTable(std::string fname) {

    ofstream out(fname.c_str());

    for (unsigned int i=0;i<table_.size();i++){
      out << i<<" "<<table_[i]<<endl;
    }
    out.close();
      
  }

  int lookup(unsigned int bin) const {
    assert(bin<table_.size());
    return table_[bin];
  }

private:

  vector<int> table_;
  
};



#endif



