# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: repr --processName=REPR --python_filename=reprocess_test_10_5_0_pre1.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 2 --era Phase2 --eventcontent FEVTDEBUGHLT --filein root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root --conditions 103X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D28 --fileout file:step2_2ev_reprocess_slim.root
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing ##parsing argument
from Configuration.StandardSequences.Eras import eras

options = VarParsing.VarParsing ('analysis')
# get and parse the command line arguments
options.parseArguments()

options.outputFile = 'D35_single_electron.root'
options.inputFiles= '/store/mc/PhaseIIMTDTDRAutumn18DR/SingleE_FlatPt-2to100/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/70000/F9B9F776-3DB1-5040-B16D-9B55CCCD3F82.root '
options.maxEvents = 100  # -1 means all events

process = cms.Process('REPR',eras.Phase2C4_trigger)
#process = cms.Process('REPR',eras.Phase2C4_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/NeutrinoGun_E_10GeV/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/280000/CF7695FC-46FE-F649-9A7E-47ABE65D3861.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/NeutrinoGun_E_10GeV/FEVT/NoPU_103X_upgrade2023_realistic_v2_ext1-v1/60000/F620EB2E-E6E2-8343-B59E-109D2A52C764.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/SingleE_FlatPt-2to100/FEVT/NoPU_103X_upgrade2023_realistic_v2-v1/70000/F580DB8B-1012-024B-9E38-04B88DA9413C.root'),

    secondaryFileNames = cms.untracked.vstring(),

#hs#
    inputCommands = cms.untracked.vstring(
        'keep *_*_*_*',
        'drop l1tEMTFHit2016*_*_*_*',
        'drop l1tEMTFTrack2016*_*_*_*'
    )
#hs#

)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('repr nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step2_2ev_reprocess_slim.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')


############################################################
# hs Customization >>>
############################################################

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile), closeFileFast = cms.untracked.bool(True))

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
#process.load("L1Trigger.TrackFindingTracklet.L1ExtendedTrackletEmulationTracks_cff")
#process.TTTracksExtendedEmulation = cms.Path(process.L1ExtendedTrackletEmulationTracks)
#process.TTTracksExtendedEmulationWithTruth = cms.Path(process.L1ExtendedTrackletEmulationTracksWithAssociators)


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
                                                                      'keep *_l1EGammaEE*_L1EG*_*'
                                                                  )
)

#process.FEVToutput_step = cms.EndPath(process.out)
#process.schedule = cms.Schedule(process.TTTracksEmulationWithTruth, process.pL1EG_EB, process.pL1EG_HGC, process.ana)
############################################################
# hs Customization <<<<
############################################################


# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
#hs#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step)
process.schedule = cms.Schedule(process.L1simulation_step,process.TTTracksEmulationWithTruth, process.pL1EG_EB, process.pL1EG_HGC,process.ana,process.endjob_step)
#process.schedule = cms.Schedule(process.L1simulation_step,process.TTTracksExtendedEmulationWithTruth, process.pL1EG_EB, process.pL1EG_HGC,process.ana,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

from L1Trigger.Configuration.customiseUtils import L1TrackTriggerTracklet,configureCSCLCTAsRun2
process = L1TrackTriggerTracklet(process)
process = configureCSCLCTAsRun2(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
