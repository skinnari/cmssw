############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os
process = cms.Process("SKIM")

 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# D49 geometry (T15 tracker)
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')

process.load('Configuration.StandardSequences.EndOfProcess_cff')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

# using this file:
# /store/relval/CMSSW_11_1_0_pre1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/110X_mcRun4_realistic_v2_2026D49noPU-v1/20000/DFF1257E-6DB1-434F-BCC2-64EC6DEFCDAF.root
inputMC = ['file:TTbar_PU0_D49_GEN-SIM-DIGI-RAW.root']
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(*inputMC),
                            inputCommands = cms.untracked.vstring(
                              'keep *_*_*_*',
                              'drop l1tEMTFHit2016*_*_*_*',
                              'drop l1tEMTFTrack2016*_*_*_*'
                              )
                            )


process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring('drop *'),
    fileName = cms.untracked.string('skimmedForCI.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    )
)

process.output.outputCommands.append('keep  *_*_*Level1TTTracks*_*')
process.output.outputCommands.append('keep  *_*_*StubAccepted*_*')
process.output.outputCommands.append('keep  *_*_*ClusterAccepted*_*')
process.output.outputCommands.append('keep  *_*_*MergedTrackTruth*_*')
process.output.outputCommands.append('keep  *_genParticles_*_*')

process.pd = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.pd)
