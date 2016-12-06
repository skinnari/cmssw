//This class implementes the projection tranceiver
#ifndef FPGAMATCHTRANSCEIVER_H
#define FPGAMATCHTRANSCEIVER_H

#include "FPGAProcessBase.hh"

using namespace std;

class FPGAMatchTransceiver:public FPGAProcessBase{

public:

  FPGAMatchTransceiver(string name, unsigned int iSector):
    FPGAProcessBase(name,iSector){
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
	output=="matchout8"||
	output=="matchout9"||
	output=="matchout10"||
	output=="matchout11"||
	output=="matchout12"||
	output=="matchout13"||
	output=="matchout14"||
	output=="matchout15"||
	output=="matchout16"||
	output=="matchout17"||
	output=="matchout18"||
	output=="matchout19"||
	output=="matchout20"||
	output=="matchout21"
	){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      outputmatches_.push_back(tmp);
      return;
    }
    cout << "In FPGAMatchTransceiver: Did not find output  = "<<output<<endl;
    assert(0);
  }

  void addInput(FPGAMemoryBase* memory,string input){
    if (writetrace) {
      cout << "In "<<name_<<" adding input from "<<memory->getName()
	   << " to input "<<input<<endl;
    }
    if (input=="matchin1"||
	input=="matchin2"||
	input=="matchin3"||
	input=="matchin4"||
	input=="matchin5"||
	input=="matchin6"||
	input=="matchin7"||
	input=="matchin8"||
	input=="matchin9"||
	input=="matchin10"||
	input=="matchin11"||
	input=="matchin12"||
	input=="matchin13"||
	input=="matchin14"||
	input=="matchin15"||
	input=="matchin16"||
	input=="matchin17"||
	input=="matchin18"||
	input=="matchin19"||
	input=="matchin20"||
	input=="matchin21"||
	input=="matchin22"||
	input=="matchin23"||
	input=="matchin24"||
	input=="matchin25"||
	input=="matchin26"||
	input=="matchin27"||
	input=="matchin28"||
	input=="matchin29"||
	input=="matchin30"||
	input=="matchin31"||
	input=="matchin32"||
	input=="matchin33"||
	input=="matchin34"||
	input=="matchin35"||
	input=="matchin36"||
	input=="matchin37"||
	input=="matchin38"||
	input=="matchin39"||
	input=="matchin40"||
	input=="matchin41"||
	input=="matchin42"||
	input=="matchin43"||
	input=="matchin44"||
	input=="matchin45"||
	input=="matchin46"||
	input=="matchin47"||
	input=="matchin48"||
	input=="matchin49"||
	input=="matchin50"
	){
      FPGAFullMatch* tmp=dynamic_cast<FPGAFullMatch*>(memory);
      assert(tmp!=0);
      inputmatches_.push_back(tmp);

      string basename=tmp->getName().substr(0,10);
      int add=0;
      for (unsigned int i=0;i<inputmatchesgroup_.size();i++){
	if (inputmatchesgroup_[i][0]->getName().substr(0,10)==basename) {
	  add++;
	  inputmatchesgroup_[i].push_back(tmp);
	}
      }
      assert(add<=1);
      if (add==0) {
	vector<FPGAFullMatch*> tmp1;
	tmp1.push_back(tmp);
	inputmatchesgroup_.push_back(tmp1);
      }
      return;
    }
    
    cout << "In FPGAMatchTransceiver: Did not find input  = "<<input<<endl;
    assert(0);
  }


  std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > orderedMatches(vector<FPGAFullMatch*>& fullmatch) {
    
    std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > tmp;


    std::vector<unsigned int> indexArray;
    for (unsigned int i=0;i<fullmatch.size();i++) {
      indexArray.push_back(0);
      //cout << "FPGAMatchTransceiver::orderedMatches "<<iSector_
      //   <<" "<<fullmatch[i]->getName()<<" "<<fullmatch[i]->nMatches();
      //for (unsigned int j=0;j<fullmatch[i]->nMatches();j++){
	//cout <<" "<<fullmatch[i]->getFPGATracklet(j)->TCID()
	//     <<" "<<fullmatch[i]->getFPGATracklet(j)->TCIDName()
	//     <<"("<<fullmatch[i]->getFPGATracklet(j)->homeSector()<<")";
	//}
      //cout<<endl;
    }

    

    int bestIndex=-1;
    do {
      int bestTCID=(1<<16);
      bestIndex=-1;
      for (unsigned int i=0;i<fullmatch.size();i++) {
	if (indexArray[i]>=fullmatch[i]->nMatches()) {
	  //skip as we were at the end
	  continue;
	}
	int TCID=fullmatch[i]->getFPGATracklet(indexArray[i])->TCID();
	if (TCID<bestTCID) {
	  bestTCID=TCID;
	bestIndex=i;
	}
      }
      if (bestIndex!=-1) {
	tmp.push_back(fullmatch[bestIndex]->getMatch(indexArray[bestIndex]));
	indexArray[bestIndex]++;
      }
    } while (bestIndex!=-1);

    //cout << "In FPGAMatchTransceiver::orderedMatches : "<<tmp.size()<<endl;
    for (unsigned int i=0;i<tmp.size();i++) {
      //cout <<" "<<tmp[i].first->TCID()<<endl;
      if (i>0) {
	//This allows for equal TCIDs. This means that we can e.g. have a track seeded
	//in L1L2 that projects to both L3 and D4. The algorithm will pick up the first hit and
	//drop the second
	assert(tmp[i-1].first->TCID()<=tmp[i].first->TCID());
      }
    }
    //cout << endl;

    return tmp;
    
  }
  

  
  
  //this->inputmatches_ to other->outputmatches_ 
  void execute(FPGAMatchTransceiver* other){
    int count=0;
    for(unsigned int i=0;i<inputmatchesgroup_.size();i++){
      string basename=inputmatchesgroup_[i][0]->getName().substr(0,10);
      std::vector<std::pair<FPGATracklet*,std::pair<FPGAStub*,L1TStub*> > > inputmatchesordered=orderedMatches(inputmatchesgroup_[i]);
      for(unsigned int l=0;l<inputmatchesordered.size();l++){
	int layer=inputmatchesordered[l].first->layer();
	int disk=inputmatchesordered[l].first->disk();
	string seed="";
	if (layer==1&&disk==0) seed="L1L2";
	if (layer==3&&disk==0) seed="L3L4";
	if (layer==5&&disk==0) seed="L5L6";
	if (layer==1&&disk==1) seed="F1L1";
	if (layer==1&&disk==-1) seed="B1L1";
	if (layer==0&&disk==1) seed="F1F2";
	if (layer==0&&disk==3) seed="F3F4";
	if (layer==0&&disk==-1) seed="B1B2";
	if (layer==0&&disk==-3) seed="B3B4";
	assert(seed!="");
	seed+=basename.substr(7,10);
	//cout << getName()<<" basename seed "<<basename<<" "<<seed<<endl;
	int nAdd=0;
	for(unsigned int j=0;j<other->outputmatches_.size();j++){
	  std::size_t found = other->outputmatches_[j]->getName().find(seed);
	  if (found!=std::string::npos){
	    if (debug1) {
	      cout << getName()<<" layer = "<<layer<<" disk = "<<disk<<" "<<seed<< " to = "
		   <<other->outputmatches_[j]->getName()<<" "<<j
		   << " "<<iSector_
		   << " "<<inputmatchesordered[l].first<<endl;
	    }      
	    nAdd++;
	    other->outputmatches_[j]->addMatch(inputmatchesordered[l]);
	  }
	}
	//cout << "nAdd = "<<nAdd<<endl;
	assert(nAdd==1);
      }
      count+=inputmatchesordered.size();
     }
    if (writeMatchTransceiver) {
      static ofstream out("matchtransceiver.txt");
      out << getName() << " " 
	  << count << endl;
    }
  }
  
  
private:

  vector<FPGAFullMatch*> inputmatches_;
  vector<vector<FPGAFullMatch*> > inputmatchesgroup_;

  vector<FPGAFullMatch*> outputmatches_;

};

#endif
