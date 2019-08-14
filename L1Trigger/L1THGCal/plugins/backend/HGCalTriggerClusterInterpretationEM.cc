
// #include <limits>
// #include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterInterpreterBase.h"
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

};

DEFINE_HGC_TPG_CLUSTER_INTERPRETER(HGCalTriggerClusterInterpretationEM, "HGCalTriggerClusterInterpretationEM");


HGCalTriggerClusterInterpretationEM::HGCalTriggerClusterInterpretationEM() {}

void HGCalTriggerClusterInterpretationEM::initialize(const edm::ParameterSet& conf) {}
void HGCalTriggerClusterInterpretationEM::interpret(// const std::vector<std::pair<GlobalPoint, double>>& seedPositionsEnergy,
                                                    // const HGCalTriggerGeometryBase& triggerGeometry,
                                                    l1t::HGCalMulticlusterBxCollection& multiclusters) const {
                                                      for(auto cluster3d: multiclusters) {
                                                        std::cout << cluster3d.detId() << " e (GeV): " << cluster3d.energy() << std::endl;
                                                      }
                                                    }
