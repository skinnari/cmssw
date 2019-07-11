############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os
process = cms.Process("L1TrackElectrons")


############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D41_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
Source_Files = cms.untracked.vstring(
    "/store/relval/CMSSW_10_6_0_pre3/RelValElectronGunPt2To100/GEN-SIM-DIGI-RAW/105X_upgrade2023_realistic_v5_2023D41noPU-v2/10000/011E8CAC-87DC-A041-83EC-2CB1F8863E45.root"
)
    
process.source = cms.Source("PoolSource", 
                            fileNames = Source_Files,
                            inputCommands = cms.untracked.vstring(
                              'keep *_*_*_*',
                              'drop l1tEMTFHit2016*_*_*_*',
                              'drop l1tEMTFTrack2016*_*_*_*'
                              )
                            )

process.TFileService = cms.Service("TFileService", fileName = cms.string('ntuple_L1TkElec.root'), closeFileFast = cms.untracked.bool(True))



############################################################
# L1 tracking
############################################################

from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *

## floating-point simulation
#process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
#process.TTTracks = cms.Path(process.L1TrackletTracks)
#process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)

## emulation 
process.load("L1Trigger.TrackFindingTracklet.L1TrackletEmulationTracks_cff")
process.TTTracksEmulation = cms.Path(process.L1TrackletEmulationTracks)
process.TTTracksEmulationWithTruth = cms.Path(process.L1TrackletEmulationTracksWithAssociators)

## Extended (displaced) emulation
# process.load("L1Trigger.TrackFindingTracklet.L1ExtendedTrackletEmulationTracks_cff")
# process.TTTracksExtendedEmulation = cms.Path(process.L1ExtendedTrackletEmulationTracks)
# process.TTTracksExtendedEmulationWithTruth = cms.Path(process.L1ExtendedTrackletEmulationTracksWithAssociators)


############################################################
# Electron producers
############################################################

# barrel (crystal-level EG)
process.load('L1Trigger.L1CaloTrigger.L1EGammaCrystalsEmulatorProducer_cfi')
process.pL1EG_EB = cms.Path( process.L1EGammaClusterEmuProducer )

# endcap (HGC)
process.load('L1Trigger.L1CaloTrigger.l1EGammaEEProducer_cfi')
process.pL1EG_HGC = cms.Path( process.l1EGammaEEProducer )

             
############################################################
# Track + Electron Ntuple
############################################################

process.L1TkElNtuple = cms.EDAnalyzer('L1TrackElectronNtupler',
                                       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       #L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                       L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
                                       #L1TrackInputTag = cms.InputTag("TTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"),
                                       GenParticleInputTag = cms.InputTag("genParticles",""),
                                       L1EGammaCrystalInputTag = cms.InputTag("L1EGammaClusterEmuProducer","L1EGammaCollectionBXVEmulator"),
                                       L1EGammaHGCInputTag = cms.InputTag("l1EGammaEEProducer","L1EGammaCollectionBXVWithCuts")
    )
process.ana = cms.Path(process.L1TkElNtuple)

# output module (use this to just dump the output collections instead of making the ntuple)
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("Tracklets.root"),
                               fastCloning = cms.untracked.bool( False ),
#                               outputCommands = cms.untracked.vstring('keep *',
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_TTTrack*_Level1TTTracks_*', 
                                                                      'keep *_L1EGammaClusterEmu*_L1EG*_*', 
                                                                      'keep *_l1EGammaEE*_L1EG*_*', 
)
)
process.FEVToutput_step = cms.EndPath(process.out)


process.schedule = cms.Schedule(process.TTTracksEmulationWithTruth, process.pL1EG_EB, process.pL1EG_HGC, process.ana)
#process.schedule = cms.Schedule(process.TTTracksEmulationWithTruth, process.pL1EG_EB, process.pL1EG_HGC, process.FEVToutput_step)

