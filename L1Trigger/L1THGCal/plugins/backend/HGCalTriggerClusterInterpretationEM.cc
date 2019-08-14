
// #include <limits>
// #include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterInterpreterBase.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
// #include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
// #include "CommonTools/MVAUtils/interface/TMVAEvaluator.h"


class HGCalTriggerClusterInterpretationEM : public HGCalTriggerClusterInterpreterBase {

public:
  HGCalTriggerClusterInterpretationEM();
  ~HGCalTriggerClusterInterpretationEM() override{};
  void initialize(const edm::ParameterSet& conf) final;
  void interpret(// const std::vector<std::pair<GlobalPoint, double>>& seedPositionsEnergy,
                         // const HGCalTriggerGeometryBase& triggerGeometry,
                         l1t::HGCalMulticlusterBxCollection& multiclusters) const final;

private:
  std::vector<double> layer_containment_corrs_;
  std::vector<double> scale_corrections_coeff_;
};

DEFINE_HGC_TPG_CLUSTER_INTERPRETER(HGCalTriggerClusterInterpretationEM, "HGCalTriggerClusterInterpretationEM");


HGCalTriggerClusterInterpretationEM::HGCalTriggerClusterInterpretationEM() {}

void HGCalTriggerClusterInterpretationEM::initialize(const edm::ParameterSet& conf) {
  layer_containment_corrs_ = conf.getParameter<std::vector<double>>("layer_containment_corr");
  scale_corrections_coeff_ = conf.getParameter<std::vector<double>>("scale_correction_coeff");
}

void HGCalTriggerClusterInterpretationEM::interpret(// const std::vector<std::pair<GlobalPoint, double>>& seedPositionsEnergy,
                                                    // const HGCalTriggerGeometryBase& triggerGeometry,
                                                    l1t::HGCalMulticlusterBxCollection& multiclusters) const {
                                                      for(auto cluster3d: multiclusters) {
                                                        double energy = cluster3d.energy();
                                                        cluster3d.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM, energy);
                                                        std::cout << cluster3d.detId()
                                                                  << " E (GeV): " << cluster3d.energy()
                                                                  << " E_EM (GeV): " << cluster3d.iEnergy(l1t::HGCalMulticluster::EnergyInterpretation::EM) << std::endl;

                                                      }
                                                    }
