#ifndef __L1Trigger_L1THGCal_HGCalTriggerClusterInterpreterBase_h__
#define __L1Trigger_L1THGCal_HGCalTriggerClusterInterpreterBase_h__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"

class HGCalTriggerClusterInterpreterBase {
public:
  HGCalTriggerClusterInterpreterBase(){};
  virtual ~HGCalTriggerClusterInterpreterBase(){};
  virtual void initialize(const edm::ParameterSet& conf) = 0;
  virtual void interpret(// const std::vector<std::pair<GlobalPoint, double>>& seedPositionsEnergy,
                         // const HGCalTriggerGeometryBase& triggerGeometry,
                         l1t::HGCalMulticlusterBxCollection& multiclusters) const = 0;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"
typedef edmplugin::PluginFactory<HGCalTriggerClusterInterpreterBase*()> HGCalTriggerClusterInterpreterFactory;

#define DEFINE_HGC_TPG_CLUSTER_INTERPRETER(type, name) DEFINE_EDM_PLUGIN(HGCalTriggerClusterInterpreterFactory, type, name)

#endif
