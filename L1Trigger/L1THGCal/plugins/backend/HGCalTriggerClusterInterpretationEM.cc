
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterInterpreterBase.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"


class HGCalTriggerClusterInterpretationEM : public HGCalTriggerClusterInterpreterBase {

public:
  HGCalTriggerClusterInterpretationEM();
  ~HGCalTriggerClusterInterpretationEM() override{};
  void initialize(const edm::ParameterSet& conf) final;
  void eventSetup(const edm::EventSetup& es)  final;
  void interpret(l1t::HGCalMulticlusterBxCollection& multiclusters) const final;

private:
  std::vector<double> layer_containment_corrs_;
  std::vector<double> scale_corrections_coeff_;
  std::vector<double> dr_bylayer_;

  HGCalTriggerTools triggerTools_;
};

DEFINE_HGC_TPG_CLUSTER_INTERPRETER(HGCalTriggerClusterInterpretationEM, "HGCalTriggerClusterInterpretationEM");


HGCalTriggerClusterInterpretationEM::HGCalTriggerClusterInterpretationEM() {}

void HGCalTriggerClusterInterpretationEM::initialize(const edm::ParameterSet& conf) {
  layer_containment_corrs_ = conf.getParameter<std::vector<double>>("layer_containment_corr");
  scale_corrections_coeff_ = conf.getParameter<std::vector<double>>("scale_correction_coeff");
  dr_bylayer_ = conf.getParameter<std::vector<double>>("dr_bylayer");
}


void HGCalTriggerClusterInterpretationEM::eventSetup(const edm::EventSetup& es) {
  triggerTools_.eventSetup(es);
}


void HGCalTriggerClusterInterpretationEM::interpret(l1t::HGCalMulticlusterBxCollection& multiclusters) const {

  for(auto cluster3d: multiclusters) {
    const GlobalPoint& cluster3d_position = cluster3d.centreProj();
    // const
    double energy = 0.;

    for(const auto& cluster2d: cluster3d.constituents()) {
      const unsigned layer = triggerTools_.layerWithOffset(cluster2d.first);
      if(layer <= 28) {
        double dr = (cluster3d_position - cluster2d.second->centreProj()).mag();
        const unsigned layer_index = (layer-1)/2;
        if(dr <= dr_bylayer_.at(layer_index)) {
          energy += layer_containment_corrs_.at(layer_index)*cluster2d.second->energy();
        }
      }
    }
    energy -= scale_corrections_coeff_.at(0)*fabs(cluster3d.eta())+scale_corrections_coeff_.at(1);
    cluster3d.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM, max(energy, 0.));
    std::cout << cluster3d.detId()
              << " E (GeV): " << cluster3d.energy()
              << " pT (GeV): " << cluster3d.pt()
              << " E_EM (GeV): " << cluster3d.iEnergy(l1t::HGCalMulticluster::EnergyInterpretation::EM)
              << " pT_EM (GeV): " << cluster3d.iPt(l1t::HGCalMulticluster::EnergyInterpretation::EM) << std::endl;

  }
}
