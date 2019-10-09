############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing ##parsing argument
import FWCore.Utilities.FileUtils as FileUtils
import os
from Configuration.StandardSequences.Eras import eras

options = VarParsing.VarParsing ('analysis')
# get and parse the command line arguments
options.parseArguments()

options.outputFile = 'D35_single_electron.root'
options.inputFiles= '/store/mc/PhaseIIMTDTDRAutumn18DR/SingleE_FlatPt-2to100/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/70000/F9B9F776-3DB1-5040-B16D-9B55CCCD3F82.root '
options.maxEvents = 100  # -1 means all events

process = cms.Process("L1TrackElectrons",eras.Phase2C4_trigger)


############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')  
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
Source_Files = cms.untracked.vstring(
    "/store/mc/PhaseIIMTDTDRAutumn18DR/SingleE_FlatPt-2to100/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/70000/F9B9F776-3DB1-5040-B16D-9B55CCCD3F82.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/SingleE_FlatPt-2to100/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/70000/F580DB8B-1012-024B-9E38-04B88DA9413C.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/SingleE_FlatPt-2to100/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/70000/EB25AB24-84EF-3146-B357-FECB5FA6DA20.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/SingleE_FlatPt-2to100/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/70000/E404FD33-007B-BB4E-965E-FCDC1E4CA11E.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/SingleE_FlatPt-2to100/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/70000/D1DC2664-4D92-BB46-94D6-F4AA215AE64F.root"
)
    
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(options.inputFiles),
                            inputCommands = cms.untracked.vstring(
                              'keep *_*_*_*',
                              'drop l1tEMTFHit2016*_*_*_*',
                              'drop l1tEMTFTrack2016*_*_*_*'
                              )
                            )

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile), closeFileFast = cms.untracked.bool(True))



############################################################
# L1 tracking
############################################################

from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

## floating-point simulation
#process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
#process.TTTracks = cms.Path(process.offlineBeamSpot*process.L1TrackletTracks)
#process.TTTracksWithTruth = cms.Path(process.offlineBeamSpot*process.L1TrackletTracksWithAssociators)

# emulation 
process.load("L1Trigger.TrackFindingTracklet.L1TrackletEmulationTracks_cff")
process.TTTracksEmulation = cms.Path(process.offlineBeamSpot*process.L1TrackletEmulationTracks)
process.TTTracksEmulationWithTruth = cms.Path(process.offlineBeamSpot*process.L1TrackletEmulationTracksWithAssociators)

## Extended (displaced) emulation
#process.load("L1Trigger.TrackFindingTracklet.L1ExtendedTrackletEmulationTracks_cff")
#process.TTTracksExtendedEmulation = cms.Path(process.offlineBeamSpot*process.L1ExtendedTrackletEmulationTracks)
#process.TTTracksExtendedEmulationWithTruth = cms.Path(process.offlineBeamSpot*process.L1ExtendedTrackletEmulationTracksWithAssociators)


############################################################
# Electron producers
############################################################

# calorimeter TPs

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)

process.load('SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cff')
process.EcalEBtp_step = cms.Path(process.simEcalEBTriggerPrimitiveDigis)

# barrel (crystal-level EG)
process.load('L1Trigger.L1CaloTrigger.L1EGammaCrystalsEmulatorProducer_cfi')
process.pL1EG_EB = cms.Path( process.L1EGammaClusterEmuProducer )

# endcap (HGC)
process.load('L1Trigger.L1CaloTrigger.l1EGammaEEProducer_cfi')
process.pL1EG_HGC = cms.Path( process.l1EGammaEEProducer )



############################################################
# add standard L1TkElectron producers
############################################################
process.load('L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi')

from L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi import L1TkElectronsCrystal
L1TkElectronsCrystal.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")
process.pL1TkElectronsCrystal = cms.Path( L1TkElectronsCrystal )

from L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi import L1TkIsoElectronsCrystal
L1TkIsoElectronsCrystal.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")
process.pL1TkIsoElectronsCrystal = cms.Path( L1TkIsoElectronsCrystal )

from L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi import L1TkElectronsHGC
L1TkElectronsHGC.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")
process.pL1TkElectronsHGC = cms.Path( L1TkElectronsHGC )

from L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi import L1TkIsoElectronsHGC
L1TkIsoElectronsHGC.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")
process.pL1TkIsoElectronsHGC = cms.Path( L1TkIsoElectronsHGC )


# Elliptical Matching
from L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi import L1TkElectronsCrystalEllipse
L1TkElectronsCrystalEllipse.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")
process.pL1TkElectronsCrystalEllipse = cms.Path( L1TkElectronsCrystalEllipse )

from L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi import L1TkIsoElectronsCrystalEllipse
L1TkIsoElectronsCrystalEllipse.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")
process.pL1TkIsoElectronsCrystalEllipse = cms.Path( L1TkIsoElectronsCrystalEllipse )

from L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi import L1TkElectronsHGCEllipse
L1TkElectronsHGCEllipse.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")
process.pL1TkElectronsHGCEllipse = cms.Path( L1TkElectronsHGCEllipse )

from L1Trigger.L1TTrackMatch.L1TkElectronTrackProducer_cfi import L1TkIsoElectronsHGCEllipse
L1TkIsoElectronsHGCEllipse.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks")
process.pL1TkIsoElectronsHGCEllipse = cms.Path( L1TkIsoElectronsHGCEllipse )


             
############################################################
# Track + Electron Ntuple
############################################################

process.L1TkElNtuple = cms.EDAnalyzer('L1TrackElectronNtupler',
                                      L1Tk_nPar = cms.int32(5), #4 or 5-parameter fit?
                                      #L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                      L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
                                      #L1TrackInputTag = cms.InputTag("TTTracksFromExtendedTrackletEmulation", "Level1TTTracks"),
                                      MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"),
                                      GenParticleInputTag = cms.InputTag("genParticles",""),
                                      L1EGammaCrystalInputTag = cms.InputTag("L1EGammaClusterEmuProducer","L1EGammaCollectionBXVEmulator"),
                                      L1EGammaHGCInputTag = cms.InputTag("l1EGammaEEProducer","L1EGammaCollectionBXVWithCuts"),
                                      tkEGBarrelInputTag = cms.InputTag("L1TkElectronsCrystal","EG"),
                                      tkEGIsoBarrelInputTag = cms.InputTag("L1TkIsoElectronsCrystal","EG"),
                                      tkEGHGCInputTag = cms.InputTag("L1TkElectronsHGC","EG"),
                                      tkEGIsoHGCInputTag = cms.InputTag("L1TkIsoElectronsHGC","EG"),
                                      tkEGBarrelEllipseInputTag = cms.InputTag("L1TkElectronsCrystalEllipse","EG"),
                                      tkEGIsoBarrelEllipseInputTag = cms.InputTag("L1TkIsoElectronsCrystalEllipse","EG"),
                                      tkEGHGCEllipseInputTag = cms.InputTag("L1TkElectronsHGCEllipse","EG"),
                                      tkEGIsoHGCEllipseInputTag = cms.InputTag("L1TkIsoElectronsHGCEllipse","EG")
)
process.ana = cms.Path(process.L1TkElNtuple)


# output module (use this to just dump the output collections instead of making the ntuple)
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("Tracklets.root"),
                               fastCloning = cms.untracked.bool( False ),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_TTTrack*_Level1TTTracks_*', 
                                                                      'keep *_L1EGammaClusterEmu*_L1EG*_*', 
                                                                      'keep *_l1EGammaEE*_L1EG*_*', 
                                                                      'keep *_L1Tk*_*_*', 
)
)
process.FEVToutput_step = cms.EndPath(process.out)

# this generates the trigger primitives from barrel ECAL + HGCal, runs 3 versions of L1 tracking (floating-point tracklet simulation, hybrid, extended hybrid), and the analyzers
process.schedule = cms.Schedule(process.hgcl1tpg_step, process.EcalEBtp_step, process.pL1EG_EB, process.pL1EG_HGC,
                                process.pL1TkElectronsCrystal, process.pL1TkIsoElectronsCrystal, process.pL1TkElectronsHGC,
                                process.pL1TkIsoElectronsHGC, process.pL1TkElectronsCrystalEllipse,
                                process.pL1TkIsoElectronsCrystalEllipse,
                                process.pL1TkElectronsHGCEllipse, process.pL1TkIsoElectronsHGCEllipse,
                                process.TTTracksEmulationWithTruth, process.ana)


