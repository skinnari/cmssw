//
// Integer representation of floating point arithmetic suitable for FPGA designs
// 
// Author: Yuri Gershtein 
// Date:   March 2018
//

#include "L1Trigger/TrackFindingTracklet/interface/imath.h"

std::string var_base::itos(int i) 
{
    std::ostringstream os;
    os << i;
    return os.str();
}

std::string var_base::get_kstring()
{

  char s[1024];
  std::string t="";
  std::map<std::string,int>::iterator it;
  for(it = Kmap_.begin(); it != Kmap_.end(); ++it){
    sprintf(s,"^(%i)",it->second);
    std::string t0(s);
    t = t + it->first + t0;
  }

  return t;
}

void var_base::analyze()
{
  if(!readytoanalyze_) return;
  
  double u = maxval_;
  if(u < -minval_)
      u = -minval_;

  int iu = log2(get_range()/u);
  if(iu>1){
    printf("analyzing %s: range %g is much larger then %g. suggest cutting by a factor of 2^%i\n",name_.c_str(),get_range(),u,iu);
  }
#ifdef IMATH_ROOT
  if(h_){
    double eff = h_->Integral()/h_->GetEntries();
    if(eff<0.99) {
      printf("analyzing %s: range is too small, contains %f\n",name_.c_str(),eff);
      h_->Print();
    }
    h_file_->cd();
    TCanvas *c = new TCanvas();
    c->cd();
    h_->Draw("colz");
    h_->Write();
  }
  else{
    if(use_root) printf("analyzing %s: no histogram!\n",name_.c_str());
  }
#endif

  if(p1_) p1_->analyze();
  if(p2_) p2_->analyze();

  readytoanalyze_ = false;
}

std::string var_base::dump()
{
  char s[1024];
  std::string u = get_kstring();
  sprintf(s,"Name = %s \t Op = %s \t nbits = %i \n       ival = %li \t fval = %g \t K = %g Range = %f\n       units = %s\n",  
	  name_.c_str(), op_.c_str(), nbits_, ival_, fval_, K_, get_range(), u.c_str());
  std::string t(s);
  return t;
}

void var_base::dump_cout()
{
  char s[2048];
  std::string u = get_kstring();
  sprintf(s,"Name = %s \t Op = %s \t nbits = %i \n       ival = %li \t fval = %g \t K = %g Range = %f\n       units = %s\n       step = %i, latency = %i\n", name_.c_str(), op_.c_str(), nbits_, ival_, fval_, K_, get_range(), u.c_str(), step_, latency_);
  std::string t(s);
  std::cout<<t;
  if(p1_) p1_->dump_cout();
  if(p2_) p2_->dump_cout();
}

bool var_base::calculate(int debug_level)
{
  bool ok1 = true;
  bool ok2 = true;
  bool ok3 = true;
  
  if(p1_) ok1 = p1_->calculate(debug_level);
  if(p2_) ok2 = p2_->calculate(debug_level);
  if(p3_) ok3 = p3_->calculate(debug_level);

  long int ival_prev = ival_;
  local_calculate();

  bool all_ok = ok1 && ok2 && ok3 && debug_level;

  if(fval_ > maxval_) maxval_ = fval_;
  if(fval_ < minval_) minval_ = fval_;
#ifdef IMATH_ROOT
  if(use_root){
    if(h_==0){
      h_file_->cd();
      std::string hname = "h_"+name_;
      h_ = (TH2F*) h_file_->Get(hname.c_str());
      if(h_ == 0){
	h_precision_ = 0.5*h_nbins_*K_;
	std::string st = name_+";fval;fval-ival*K";
	h_ = new TH2F(hname.c_str(),name_.c_str(),
		      h_nbins_,-get_range(), get_range(),
		      h_nbins_,-h_precision_, h_precision_);
	if(debug_level==3) std::cout<<" booking histogram "<<hname<<"\n";
      }
    }
    if(ival_ != ival_prev || op_=="def" || op_=="const") h_->Fill(fval_, K_*ival_-fval_);
  }
#endif

  bool todump = false;
  int nmax = sizeof(long int)*8;
  int ns = nmax - nbits_;
  long int itest = ival_;
  itest = itest<<ns;
  itest = itest>>ns;
  if(itest!=ival_){
    if(debug_level == 3 || (ival_!=ival_prev && all_ok)){
      std::cout<<"imath: truncated value mismatch!! "<<ival_<<" != "<<itest<<"\n";
      todump = true;
    }
    all_ok = false;
  }
  
  val_ = ival_ * K_;
  float ftest = val_ ;
  float tolerance = 0.1 * fabs(fval_);
  if(tolerance < 2 * K_) tolerance = 2 * K_;
  if(fabs(ftest-fval_)> tolerance){
    if( debug_level == 3 || (ival_!=ival_prev &&(all_ok && (op_!="inv" ||debug_level>=2 )))){
      std::cout<<"imath: **GROSS** value mismatch!! "<<fval_<<" != "<<ftest<<"\n";
      if(op_=="inv") std::cout<<p1_->dump()<<"\n-----------------------------------\n";
      todump = true;
    }
    all_ok = false;
  }
  
  if(todump)
    std::cout<<dump();

  return all_ok;
}

void var_inv::writeLUT(std:: ofstream& fs)
{
  for(int i=0; i<Nelements_; ++i)
      fs<<std::hex<<(LUT[i]&((1<<nbits_)-1))<<std::dec<<"\n";
}

void var_adjustK::adjust(double Knew, double epsilon, bool do_assert, int nbits)
{
  //WARNING!!!
  //THIS METHID CAN BE USED ONLY FOR THE FINAL ANSWER
  //THE CHANGE IN CONSTANT CAN NOT BE PROPAGATED UP THE CALCULATION TREE
  
    K_     = p1_->get_K();
    Kmap_  = p1_->get_Kmap();
    double r = Knew / K_;

    lr_ = (r>1)? log2(r)+epsilon : log2(r);
    K_ = K_ * pow(2,lr_);
    if(do_assert) assert(fabs(Knew/K_ - 1)<epsilon);
    
    if(nbits>0)
      nbits_ = nbits;
    else
      nbits_ = p1_->get_nbits()-lr_;

    Kmap_["2"] = Kmap_["2"] + lr_;
    
}

//
//  local calculations
//

void var_adjustK::local_calculate()
{
  fval_ = p1_->get_fval();
  ival_ = p1_->get_ival();
  if(lr_>0)
    ival_ = ival_ >> lr_;
  else if(lr_<0)
    ival_ = ival_ <<(-lr_);
}
void var_adjustKR::local_calculate()
{
  fval_ = p1_->get_fval();
  ival_ = p1_->get_ival();
  if(lr_>0)
    ival_ = ((ival_ >> (lr_-1))+1)>>1; //rounding
  else if(lr_<0)
    ival_ = ival_ <<(-lr_);
}
void var_add::local_calculate()
{
  fval_ = p1_->get_fval() + p2_->get_fval();
  long int i1 = p1_->get_ival();
  long int i2 = p2_->get_ival();
  if(shift1>0) i1 = i1 << shift1;
  if(shift2>0) i2 = i2 << shift2;
  ival_ = i1 + i2;
  if(ps_>0) ival_ = ival_ >> ps_;
}
void var_subtract::local_calculate()
{
  fval_ = p1_->get_fval() - p2_->get_fval();
  long int i1 = p1_->get_ival();
  long int i2 = p2_->get_ival();
  if(shift1>0) i1 = i1 << shift1;
  if(shift2>0) i2 = i2 << shift2;
  ival_ = i1 - i2;
  if(ps_>0) ival_ = ival_ >> ps_;
}
void var_nounits::local_calculate()
{
  fval_ = p1_->get_fval();
  ival_ = (p1_->get_ival() * cI_)>>ps_;
}
void var_timesC::local_calculate()
{
  fval_ = p1_->get_fval() * cF_;
  ival_ = (p1_->get_ival() * cI_)>>ps_;
}
void var_neg::local_calculate()
{
  fval_ = -p1_->get_fval();
  ival_ = -p1_->get_ival();
}
void var_shift::local_calculate()
{
  fval_ = p1_->get_fval() * pow(2,-shift_);
  ival_ = p1_->get_ival();
  if(shift_>0) ival_ = ival_>>shift_;
  if(shift_<0) ival_ = ival_<<(-shift_);
}
void var_shiftround::local_calculate()
{
  fval_ = p1_->get_fval() * pow(2,-shift_);
  ival_ = p1_->get_ival();
  if(shift_>0) ival_ = ((ival_>>(shift_-1))+1)>>1;
  if(shift_<0) ival_ = ival_<<(-shift_);
}
void var_mult::local_calculate()
{
  fval_ = p1_->get_fval() * p2_->get_fval();
  ival_ = (p1_->get_ival() * p2_->get_ival())>>ps_;
}

void var_DSP_postadd::local_calculate()
{
  fval_ = p1_->get_fval() * p2_->get_fval() + p3_->get_fval();
  ival_ = p3_->get_ival();
  if(shift3_>0) ival_ = ival_<<shift3_;
  if(shift3_<0) ival_ = ival_>>(-shift3_);
  ival_ += p1_->get_ival() * p2_->get_ival();
  ival_ = ival_>>ps_;
}

void var_inv::initLUT(double offset)
{
  offset_ = offset;
  double offsetI = round_int(offset_ / p1_->get_K());
  for(int i=0; i<Nelements_; ++i){
    int i1 = addr_to_ival(i);
    LUT[i] = gen_inv(offsetI+i1);
  }
}
void var_inv::local_calculate()
{
  fval_ = 1./(offset_ + p1_->get_fval());
  ival_ = LUT[ival_to_addr(p1_->get_ival())];
}
//
// print functions
//

void var_base::makeready()
{
  pipe_counter_ = 0;
  pipe_delays_.clear();
  readytoprint_   = true;
  readytoanalyze_ = true;
  usedasinput_    = false;
  if(p1_) p1_->makeready();
  if(p2_) p2_->makeready();
  if(p3_) p3_->makeready();
}

bool var_base::has_delay(int i)
{
  //dumb sequential search
  for(unsigned int j=0; j<pipe_delays_.size(); ++j)
    if(pipe_delays_[j] == i) return true;
  return false;
}

std::string var_base::pipe_delay(var_base *v, int nbits, int delay)
{
  //have we been delayed by this much already?
  if(v->has_delay(delay)) return "";
  v->add_delay(delay);
  std::string name = v->get_name();
  std::string name_delayed = name+"_delay"+itos(delay);
  std::string out = "wire signed ["+itos(nbits-1)+":0] "+name_delayed+";\n";
  out = out + pipe_delay_wire(v, name_delayed, nbits, delay);
  return out;
}
std::string var_base::pipe_delays(const int step)
{
  std::string out = "";
  if (p1_) out += p1_->pipe_delays (step);
  if (p2_) out += p2_->pipe_delays (step);
  if (p3_) out += p3_->pipe_delays (step);

  int l = step - latency_ - step_;
  return (out + pipe_delay(this,get_nbits(),l));
}
std::string var_base::pipe_delay_wire(var_base *v, std::string name_delayed, int nbits, int delay)
{
  std::string name = v->get_name();
  std::string name_pipe    = name+"_pipe"+itos(v->pipe_counter());
  v->pipe_increment();
  std::string out = "pipe_delay #(.STAGES("+itos(delay)+"), .WIDTH("+itos(nbits)+")) "
    + name_pipe + "(.clk(clk), .val_in("+name+"), .val_out("+name_delayed+"));\n";
  return out;
}

void var_base::print_step(int step, std::ofstream& fs){
  if(!readytoprint_) return;
  if(step > step_) return;
  int l1 = 0;
  int l2 = 0;
  int l3 = 0;
  if(p1_) {
    p1_->print_step(step, fs);
    l1 = step_ - p1_->get_latency() - p1_->get_step();
  }
  if(p2_) {
    p2_->print_step(step, fs);
    l2 = step_ - p2_->get_latency() - p2_->get_step();
  }
  if(p3_) {
    p3_->print_step(step, fs);
    l3 = step_ - p3_->get_latency() - p3_->get_step();
  }
  if(step==step_){
    if(l1<0 || l2<0 ||l3<0 || (l1>0&&l2>0&&l3>0) ){
      printf("%s::print_step(%i): something wrong with latencies! %i %i %i\n",name_.c_str(),step,l1, l2, l3);
      dump_cout();
      assert(0);      
    }
    if(l1>0) {
      if(p1_->get_op()!="const")
	fs<<pipe_delay(p1_,p1_->get_nbits(),l1);
      else
	l1 = 0;
    }
    if(l2>0) {
      if(p2_->get_op()!="const")
	fs<<pipe_delay(p2_,p2_->get_nbits(),l2);
      else
	l2 = 0;
    }
    if(l3>0) {
      if(p3_->get_op()!="const")
	fs<<pipe_delay(p3_,p3_->get_nbits(),l3);
      else
	l3 = 0;
    }

    if (op_ == "flag"){
      for (const auto &cut : cuts_)
        fs<<cut->get_cut_var()->pipe_delays(step_);
    }
  
    print(fs, l1, l2, l3);
    readytoprint_ = false;
  }
}

void var_base::print_all(std::ofstream& fs)
{
  for(int i=0; i<=step_; ++i){
    fs<<"//\n// STEP "<<i<<"\n\n";
    print_step(i,fs);
  }
}

void var_base::print_truncation(std::string &t, const std::string &o1, const int ps) const
{
  if (ps > 0) {
    t += "wire signed ["+itos(nbits_-1)+":0]"+name_+";\n";
    t += "reg signed  ["+itos(nbits_+ps-1)+":0]"+name_+"_tmp;\n";
    t += "always @(posedge clk) "+name_+"_tmp <= "+o1+";\n";
    t += "assign "+name_+" = "+name_+"_tmp["+itos(nbits_+ps-1)+":"+itos(ps)+"];\n";
  }
  else {
    t += "reg signed  ["+itos(nbits_-1)+":0]"+name_+";\n";
    t += "always @(posedge clk) "+name_+" <= "+o1+";\n";
  }
}

void var_base::get_inputs(std::vector<var_base*> *vd)
{
  if(op_ == "def" && !usedasinput_){
    usedasinput_ = true;
    vd->push_back(this);
  }
  else{
    if(p1_) p1_->get_inputs(vd);
    if(p2_) p2_->get_inputs(vd);
    if(p3_) p3_->get_inputs(vd);
  }
}

#ifdef IMATH_ROOT
TTree* var_base::AddToTree(var_base* v, char *s)
{
  if(h_file_==0){
    h_file_ = new TFile("imath.root","RECREATE");
    printf("recreating file imath.root\n");
  }
  h_file_->cd();
  TTree *tt = (TTree*) h_file_->Get("tt");
  if(tt==0){
    tt = new TTree("tt","");
    printf("creating TTree tt\n");
  }
  std::string si = v->get_name()+"_i";
  std::string sf = v->get_name()+"_f";
  std::string sv = v->get_name();
  if(s!=0){
    std::string prefix(s);
    si = prefix + si;
    sf = prefix + sf;
    sv = prefix + sv;
  }
  if(!tt->GetBranchStatus(si.c_str())){
    tt->Branch(si.c_str(),(Long64_t*) &(v->ival_));
    tt->Branch(sf.c_str(),&(v->fval_));
    tt->Branch(sv.c_str(),&(v->val_));
  }

  if(v->p1_) AddToTree(v->p1_, s);
  if(v->p2_) AddToTree(v->p2_, s);
  if(v->p3_) AddToTree(v->p3_, s);
  
  return tt;
}
TTree* var_base::AddToTree(double* v, char *s)
{
  if(h_file_==0){
    h_file_ = new TFile("imath.root","RECREATE");
    printf("recreating file imath.root\n");
  }
  h_file_->cd();
  TTree *tt = (TTree*) h_file_->Get("tt");
  if(tt==0){
    tt = new TTree("tt","");
    printf("creating TTree tt\n");
  }
  tt->Branch(s,v);
  return tt;
}
TTree* var_base::AddToTree(int* v, char *s)
{
  if(h_file_==0){
    h_file_ = new TFile("imath.root","RECREATE");
    printf("recreating file imath.root\n");
  }
  h_file_->cd();
  TTree *tt = (TTree*) h_file_->Get("tt");
  if(tt==0){
    tt = new TTree("tt","");
    printf("creating TTree tt\n");
  }
  tt->Branch(s,v);
  return tt;
}
void var_base::FillTree()
{
  if(h_file_==0) return;
  h_file_->cd();
  TTree *tt = (TTree*) h_file_->Get("tt");
  if(tt==0) return;
  tt->Fill();
}
void var_base::WriteTree()
{
  if(h_file_==0) return;
  h_file_->cd();
  TTree *tt = (TTree*) h_file_->Get("tt");
  if(tt==0) return;
  tt->Write();
}
  
#endif


void var_base::Verilog_print(std::vector<var_base*> v, std::ofstream& fs)
{

  //step at which all the outputs should be valid
  int maxstep = 0;
  
  //header of the module

  //inputs
  std::vector<var_base*> vd;
  vd.clear();
  int imax = v.size();
  for(int i=0; i<imax; ++i){
    (v[i])->get_inputs(&vd);
    int step = v[i]->get_step() + v[i]->get_latency();
    if( step > maxstep) maxstep = step;
  }

  //print header
  fs<<"module \n";
  fs<<"(\n";
  fs<<"   input clk,\n";
  fs<<"   input reset,\n\n";

  imax = vd.size();
  for(int i=0; i<imax; ++i)
    fs<<"   input ["<<(vd[i])->get_nbits()-1<<":0] "<<(vd[i])->get_name()<<"_wire,\n";
  fs<<"\n";

  imax = v.size()-1;
  for(int i=0; i<imax; ++i)
    if (v[i]->get_nbits() > 1)
      fs<<"   output ["<<(v[i])->get_nbits()-1<<":0] "<<(v[i])->get_name()<<"_wire,\n";
    else
      fs<<"   output "<<(v[i])->get_name()<<"_wire,\n";
  if(imax>=0){
    if (v[imax]->get_nbits() > 1)
      fs<<"   output ["<<(v[imax])->get_nbits()-1<<":0] "<<(v[imax])->get_name()<<"_wire\n";
    else
      fs<<"   output "<<(v[imax])->get_name()<<"_wire\n";
  }
  fs<<");\n\n";

  //body of the module
  imax = v.size();
  for(int i=0; i<imax; ++i){
    fs<<"\n//\n";
    fs<<"// calculating "<<(v[i])->get_name()<<"\n";
    fs<<"//\n";
    (v[i])->print_all(fs);
  }
  fs<<"\n";

  //trailer
  fs<<"\n";  
  fs<<"\n//\n";
  fs<<"// wiring the outputs \n";
  fs<<"// latency = "<<maxstep<<"\n";
  fs<<"//\n";
  for(int i=0; i<imax; ++i){
    std::string n = v[i]->get_name()+"_wire";
    int delay = maxstep - v[i]->get_step()-v[i]->get_latency();
    if(delay==0)
      fs<<"assign "<< n <<" = "<<(v[i])->get_name()<<";\n";
    else
      fs<<pipe_delay_wire(v[i],n, v[i]->get_nbits(), delay);
  }
  
  fs<<"endmodule\n";
}

void var_cut::local_passes(std::map<const var_base *,std::vector<bool> > &passes, const std::map<const var_base *,std::vector<bool> > * const previous_passes) const {
  const int lower_cut = lower_cut_/cut_var_->get_K();
  const int upper_cut = upper_cut_/cut_var_->get_K();
  if (!previous_passes || (previous_passes && !previous_passes->count(cut_var_))){
    if (!passes.count(cut_var_)) passes[cut_var_];
    passes.at(cut_var_).push_back(cut_var_->get_ival() > lower_cut && cut_var_->get_ival() < upper_cut);
  }
}

bool var_base::local_passes() const {
  bool passes = false;
  for (const auto &cut : cuts_){
    const var_cut * const cast_cut = (var_cut *) cut;
    const int lower_cut = cast_cut->get_lower_cut()/K_;
    const int upper_cut = cast_cut->get_upper_cut()/K_;
    passes = passes || (ival_ > lower_cut && ival_ < upper_cut);
  }
  return passes;
}

void var_base::passes(std::map<const var_base *,std::vector<bool> > &passes, const std::map<const var_base *,std::vector<bool> > * const previous_passes) const {
  if (p1_) p1_->passes(passes,previous_passes);
  if (p2_) p2_->passes(passes,previous_passes);
  if (p3_) p3_->passes(passes,previous_passes);

  for (const auto &cut : cuts_){
    const var_cut * const cast_cut = (var_cut *) cut;
    const int lower_cut = cast_cut->get_lower_cut()/K_;
    const int upper_cut = cast_cut->get_upper_cut()/K_;
    if (!previous_passes || (previous_passes && !previous_passes->count(this))){
      if (!passes.count(this)) passes[this];
      passes.at(this).push_back(ival_ > lower_cut && ival_ < upper_cut);
    }
  }
}

void var_cut::print(std::map<const var_base *,std::set<std::string> > &cut_strings, const int step, const std::map<const var_base *,std::set<std::string> > * const previous_cut_strings) const {
  int l = step - cut_var_->get_latency() - cut_var_->get_step();
  std::string name = cut_var_->get_name();
  if (l>0) name += "_delay"+itos(l);

  const int lower_cut = lower_cut_/cut_var_->get_K();
  const int upper_cut = upper_cut_/cut_var_->get_K();
  if (!previous_cut_strings || (previous_cut_strings && !previous_cut_strings->count(cut_var_))){
    if (!cut_strings.count(cut_var_)) cut_strings[cut_var_];
    cut_strings.at(cut_var_).insert(name+" > "+itos(lower_cut)+" && "+name+" < "+itos(upper_cut));
  }
}

void var_base::print_cuts(std::map<const var_base *,std::set<std::string> > &cut_strings, const int step, const std::map<const var_base *,std::set<std::string> > * const previous_cut_strings) const {
  if (p1_) p1_->print_cuts(cut_strings,step,previous_cut_strings);
  if (p2_) p2_->print_cuts(cut_strings,step,previous_cut_strings);
  if (p3_) p3_->print_cuts(cut_strings,step,previous_cut_strings);

  int l = step - latency_ - step_;
  std::string name = name_;
  if (l>0) name += "_delay"+itos(l);

  for (const auto &cut : cuts_){
    const var_cut * const cast_cut = (var_cut *) cut;
    const int lower_cut = cast_cut->get_lower_cut()/K_;
    const int upper_cut = cast_cut->get_upper_cut()/K_;
    if (!previous_cut_strings || (previous_cut_strings && !previous_cut_strings->count(this))){
      if (!cut_strings.count(this)) cut_strings[this];
      cut_strings.at(this).insert(name+" > "+itos(lower_cut)+" && "+name+" < "+itos(upper_cut));
    }
  }
}

void var_base::add_cut(var_cut *cut, const bool call_set_cut_var){
  cuts_.push_back(cut);
  if (call_set_cut_var) cut->set_cut_var(this,false);
}

void var_cut::set_cut_var(var_base *cut_var, const bool call_add_cut){
  cut_var_ = cut_var;
  if (call_add_cut) cut_var->add_cut(this,false);
  if (parent_flag_) parent_flag_->calculate_step();
}

void var_flag::add_cut(var_base *cut, const bool call_set_parent_flag){
  cuts_.push_back(cut);
  if (cut->get_op() == "cut" && call_set_parent_flag){
    var_cut * const cast_cut = (var_cut *) cut;
    cast_cut->set_parent_flag(this,false);
  }
  calculate_step();
}

void var_cut::set_parent_flag(var_flag *parent_flag, const bool call_add_cut){
  parent_flag_ = parent_flag;
  if (call_add_cut) parent_flag->add_cut(this,false);
}

var_base * var_base::get_cut_var() {
  if (op_ == "cut")
    return cut_var_;
  else
    return this;
}


void var_adjustK::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(l2==0);
  assert(l3==0);

  std::string shift = "";
  if(lr_>0)
    shift = " >>> " + itos(lr_);
  else if(lr_<0)
    shift = " << " + itos(-lr_);

  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);

  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"; 
  std::string t = "wire signed ["+itos(nbits_-1)+":0]"+name_+";\n";
  t += "assign "+name_+" = "+n1+shift;
  fs<<t<<"; \n";
  
}

void var_adjustKR::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(l2==0);
  assert(l3==0);

  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);

  std::string  o1 = n1;
  if(lr_==1) o1 = "("+o1+"+1)>>>1";
  if(lr_>1)  o1 = "( ("+o1 + ">>>" + itos(lr_-1)+")+1)>>>1";
  if(lr_<0)  o1 = o1 + "<<" + itos(-lr_);
  
  std::string t = "reg signed ["+itos(nbits_-1)+":0]"+name_+";\n";
  t += "always @(posedge clk) "+name_+" <= "+o1;
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t<<";\n"; 
  
}

void var_def::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(l1==0);
  assert(l2==0);
  assert(l3==0);

  std::string t = "reg signed  ["+itos(nbits_-1)+":0]"+name_+";\n";
  t = t + "always @(posedge clk) "+name_+" <= "+name_+"_wire;\n";
  fs<<"// units "<<get_kstring()<<"\t"<<K_<<"\n"<<t; 
}

void var_param::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(l1==0);
  assert(l2==0);
  assert(l3==0);
  std::string t = "parameter "+ name_ + " = ";
  if(ival_<0)
    t = t + "- " + itos(nbits_) + "\'sd"+ itos(-ival_);
  else
    t = t + itos(nbits_) + "\'sd"+ itos(ival_);
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t<<";\n"; 
}

void var_add::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(p2_);
  assert(l3==0);
  std::string o1 = p1_->get_name();
  if(l1>0) o1 += "_delay"+itos(l1);
  if(shift1>0) {
    o1 += "<<"+itos(shift1);
    o1 = "("+o1+")";
  }
  
  std::string o2 = p2_->get_name();
  if(l2>0) o2 += "_delay"+itos(l2);
  if(shift2>0) {
    o2 += "<<"+itos(shift2);
    o2 = "("+o2+")";
  }

  o1 = o1 + " + " + o2; 
  
  std::string t = "";
  print_truncation(t, o1, ps_);
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t;
}

void var_subtract::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(p2_);
  assert(l3==0);
  std::string o1 = p1_->get_name();
  if(l1>0) o1 += "_delay"+itos(l1);
  if(shift1>0) {
    o1 += "<<"+itos(shift1);
    o1 = "("+o1+")";
  }
  
  std::string o2 = p2_->get_name();
  if(l2>0) o2 += "_delay"+itos(l2);
  if(shift2>0) {
    o2 += "<<"+itos(shift2);
    o2 = "("+o2+")";
  }

  o1 = o1 + " - " + o2; 
  
  std::string t = "";
  print_truncation(t, o1, ps_);
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t;
}

void var_nounits::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(l2==0);
  assert(l3==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string o1 = "(" + n1 + " * " + itos(cI_) + ")";

  std::string t = "";
  print_truncation(t, o1, ps_);
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t;
}

void var_timesC::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(l2==0);
  assert(l3==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string o1 = "(" + n1 + " * " + itos(cI_) + ")";

  std::string t = "";
  print_truncation(t, o1, ps_);
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t;
}
void var_neg::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(l2==0);
  assert(l3==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);

  std::string t = "reg signed  ["+itos(nbits_-1)+":0]"+name_+";\n";
  t += "always @(posedge clk) "+name_+" <= - "+n1;
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t<<";\n"; 
}
void var_shiftround::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(l2==0);
  assert(l3==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string  o1 = n1;
  if(shift_==1) o1 = "("+o1+"+1)>>>1";
  if(shift_>1)  o1 = "( ("+o1 + ">>>" + itos(shift_-1)+")+1)>>>1";
  if(shift_<0)  o1 = o1 + "<<" + itos(-shift_);

  std::string t = "reg signed ["+itos(nbits_-1)+":0]"+name_+";\n";
  t += "always @(posedge clk) "+name_+" <= "+o1;
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t<<";\n"; 
}
void var_shift::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(l2==0);
  assert(l3==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string o1  = n1;
  if(shift_>0) o1 = o1 + ">>>" + itos(shift_);
  if(shift_<0) o1 = o1 + "<<" + itos(-shift_);

  std::string t = "wire signed ["+itos(nbits_-1)+":0]"+name_+";\n";
  t += "assign "+name_+" = "+o1;
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t<<";\n"; 
}
void var_mult::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(l3==0);
  assert(p1_);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  assert(p2_);
  std::string n2 = p2_->get_name();
  if(l2>0) n2 = n2 + "_delay"+itos(l2);
  std::string o1 =  n1 + " * " + n2;

  std::string t = "";
  print_truncation(t, o1, ps_);
  fs<<"// "<<nbits_ <<" bits \t "<<get_kstring()<<"\t"<<K_<<"\n"<<t;
}
void var_inv::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(l2==0);
  assert(l3==0);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  //first calculate address
  std::string t1 = "addr_" + name_;
  std::string t = "wire ["+itos(nbaddr_-1)+":0] "+t1+";\n";
  t =   t + "assign "+t1+" = ";
  if(shift_>0)
    t = t + "(" + n1 + ">>>"+itos(shift_)+") & "+itos(mask_);
  else
    t = t + n1 + " & "+itos(mask_);
  fs<<t<<"; // address for the LUT\n"; 

  t = "wire signed ["+ itos(nbits_-1) + ":0] "+ name_ + ";\n";
  fs<<t;

  std::string t2 = "LUT_" + name_;

  fs<<"Memory #( \n";
  fs<<"         .RAM_WIDTH("<<nbits_<<"),            // Specify RAM data width \n";
  fs<<"         .RAM_DEPTH("<<Nelements_<<"),                     // Specify RAM depth (number of entries) \n";
  fs<<"         .RAM_PERFORMANCE(\"HIGH_PERFORMANCE\"), // \"HIGH_PERFORMANCE\" = 2 clks latency \n";
  fs<<"         .INIT_FILE() \n";
  fs<<"       ) "<< t2 <<" ( \n";
  fs<<"         .addra("<<itos(nbaddr_)<<"\'b0),    // Write address bus, width determined from RAM_DEPTH  \n";
  fs<<"         .addrb("<<t1<<" ),                   // Read address bus, width determined from RAM_DEPTH  \n";
  fs<<"         .dina("<<itos(nbits_)<<"\'b0),      // RAM input data, width determined from RAM_WIDTH   \n";
  fs<<"         .clka(clk),      // Write clock \n";
  fs<<"         .clkb(clk),      // Read clock  \n";
  fs<<"         .wea(1\'b0),        // Write enable  \n";
  fs<<"         .enb(1\'b1),        // Read Enable, for additional power savings, disable when not in use  \n";
  fs<<"         .rstb(reset),      // Output reset (does not affect memory contents)                      \n";
  fs<<"         .regceb(1\'b1),  // Output register enable                                                \n";
  fs<<"         .doutb("<<name_<<")     // RAM output data,                                                \n";
  fs<<"     ); \n";


}

void var_DSP_postadd::print(std::ofstream& fs, int l1, int l2, int l3)
{
  assert(p1_);
  assert(p2_);
  assert(p3_);
  std::string n1 = p1_->get_name();
  if(l1>0) n1 = n1 + "_delay"+itos(l1);
  std::string n2 = p2_->get_name();
  if(l2>0) n2 = n2 + "_delay"+itos(l2);
  std::string n3 = p3_->get_name();
  if(l3>0) n3 = n3 + "_delay"+itos(l3);

  if(shift3_ >0) n3 = n3 + "<<"+itos(shift3_);
  if(shift3_ <0) n3 = n3 + ">>>"+itos(-shift3_);

  std::string n4 = "";
  if(ps_>0) n4 = ">>>"+itos(ps_);
  
  fs<<name_+" = DSP_postadd("+n1+", "+n2+", "+n3+")"+n4+";";
  
}

void var_flag::calculate_step(){
  int max_step = 0;
  for (const auto &cut : cuts_){
    if (!cut->get_cut_var()) continue;
    if (cut->get_cut_var()->get_latency()+cut->get_cut_var()->get_step() > max_step)
      max_step = cut->get_cut_var()->get_latency()+cut->get_cut_var()->get_step();
  }
  step_ = max_step;
}

bool var_flag::passes(){
  std::map<const var_base *,std::vector<bool> > passes0, passes1;
  for (const auto &cut : cuts_){
    if (cut->get_op() != "cut") continue;
    const var_cut * const cast_cut = (var_cut *) cut;
    cast_cut->local_passes(passes0);
  }
  for (const auto &cut : cuts_){
    if (cut->get_op() != "cut")
      cut->passes(passes1,&passes0);
    else{
      if (cut->get_cut_var()->get_p1()) cut->get_cut_var()->get_p1()->passes(passes1,&passes0);
      if (cut->get_cut_var()->get_p2()) cut->get_cut_var()->get_p2()->passes(passes1,&passes0);
      if (cut->get_cut_var()->get_p3()) cut->get_cut_var()->get_p3()->passes(passes1,&passes0);
    }
  }

  bool passes = true;
  for (const auto &cut_var : passes0){
    bool local_passes = false;
    for (const auto &pass : cut_var.second)
      local_passes = local_passes || pass;
    passes = passes && local_passes;
  }
  for (const auto &cut_var : passes1){
    bool local_passes = false;
    for (const auto &pass : cut_var.second)
      local_passes = local_passes || pass;
    passes = passes && local_passes;
  }

  return passes;
}

void var_flag::print(std::ofstream& fs, int l1, int l2, int l3){
  assert(l1==0);
  assert(l2==0);
  assert(l3==0);

  fs<<"wire "<<name_<<";"<<std::endl;
  fs<<"assign "<<name_<<" = (";
  std::map<const var_base *,std::set<std::string> > cut_strings0, cut_strings1;
  for (const auto &cut : cuts_){
    if (cut->get_op() != "cut") continue;
    const var_cut * const cast_cut = (var_cut *) cut;
    cast_cut->print(cut_strings0,step_);
  }
  for (const auto &cut : cuts_){
    if (cut->get_op() != "cut")
      cut->print_cuts(cut_strings1,step_,&cut_strings0);
    else{
      if (cut->get_cut_var()->get_p1()) cut->get_cut_var()->get_p1()->print_cuts(cut_strings1,step_,&cut_strings1);
      if (cut->get_cut_var()->get_p2()) cut->get_cut_var()->get_p2()->print_cuts(cut_strings1,step_,&cut_strings1);
      if (cut->get_cut_var()->get_p3()) cut->get_cut_var()->get_p3()->print_cuts(cut_strings1,step_,&cut_strings1);
    }
  }

  std::string separator = "";
  for (const auto &cut_var : cut_strings0){
    separator += "((";
    for (const auto &cut_string : cut_var.second){
      fs<<separator<<cut_string;
      separator = ") || (";
    }
    separator = ")) && ";
  }
  for (const auto &cut_var : cut_strings1){
    separator += "((";
    for (const auto &cut_string : cut_var.second){
      fs<<separator<<cut_string;
      separator = ") || (";
    }
    separator = ")) && ";
  }

  fs<<")));";
}


