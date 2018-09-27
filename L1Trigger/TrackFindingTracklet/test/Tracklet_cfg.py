# define basic process
import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1Tracklet")
 

# import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff') ## this needs to match the geometry you are running on
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')     ## this needs to match the geometry you are running on

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


# input
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
Source_Files = cms.untracked.vstring(
"/store/relval/CMSSW_9_3_7/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/5A8CFF7F-1E2D-E811-A7B0-0242AC130002.root"
)
process.source = cms.Source("PoolSource", fileNames = Source_Files,
                            inputCommands = cms.untracked.vstring( #this is needed to run on 9_3_X samples in 10_X_Y releases
                                'keep *_*_*_*',
                                'drop l1tEMTFHit2016*_*_*_*',
                                'drop l1tEMTFTrack2016*_*_*_*'
                                )
)


## L1 tracking: floating-point
#process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
#process.TTTracks = cms.Path(process.L1TrackletTracks)                         #run only the tracking (no MC truth associators)
#process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators) #run the tracking AND MC truth associators)

# L1 tracking: integer emulation
process.load("L1Trigger.TrackFindingTracklet.L1TrackletEmulationTracks_cff")
process.TTTracksEmulation = cms.Path(process.L1TrackletEmulationTracks)                         #run only the tracking (no MC truth associators)
process.TTTracksEmulationWithTruth = cms.Path(process.L1TrackletEmulationTracksWithAssociators) #run the tracking AND MC truth associators)


# output module
process.out = cms.OutputModule( "PoolOutputModule",
                                fileName = cms.untracked.string("Tracklets.root"),
                                fastCloning = cms.untracked.bool( False ),
                                outputCommands = cms.untracked.vstring('drop *',
                                                                       'keep *_TTTrack*_Level1TTTracks_*', 
#                                                                       'keep *_TTCluster*_*_*',
#                                                                       'keep *_TTStub*_*_*'
)
)
process.FEVToutput_step = cms.EndPath(process.out)

#process.schedule = cms.Schedule(process.TTTracksWithTruth,process.FEVToutput_step)
process.schedule = cms.Schedule(process.TTTracksEmulationWithTruth,process.FEVToutput_step)

