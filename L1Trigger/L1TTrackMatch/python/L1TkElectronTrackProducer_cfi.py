import FWCore.ParameterSet.Config as cms

L1TkElectrons = cms.EDProducer("L1TkElectronTrackProducer",
	    label = cms.string("EG"),	# labels the collection of L1TkEmParticleProducer that is produced.
					# (not really needed actually)
        L1EGammaInputTag = cms.InputTag("simCaloStage2Digis",""),
        ETmin = cms.double( -1.0 ),             # Only the L1EG objects that have ET > ETmin in GeV
                                                # are considered. ETmin < 0 means that no cut is applied.
     	L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
        # Quality cuts on Track and Track L1EG matching criteria                                
        TrackChi2           = cms.double(1e10), # minimum Chi2 to select tracks
        TrackMinPt          = cms.double(10.0), # minimum Pt to select tracks                                     
	    useTwoStubsPT       = cms.bool( False ),
        TrackEGammaDeltaPhi = cms.vdouble(0.07, 0.0, 0.0), # functional Delta Phi cut parameters to match Track with L1EG objects
        TrackEGammaDeltaR   = cms.vdouble(0.08, 0.0, 0.0), # functional Delta R cut parameters to match Track with L1EG objects
        TrackEGammaDeltaEta = cms.double(1e10), # Delta Eta cutoff to match Track with L1EG objects
        TrackEGammaDeltaPhimax = cms.vdouble(0.0, 0.0), # functional Delta Phi maxium parameters to match Track with L1EG objects
        TrackEGammaDeltaEtamax   = cms.vdouble(0.0, 0.0), # functional Delta Eta maxium parameters to match Track with L1EG objects
                                                # are considered. (unused in default configuration)
        TrackEGammaPhiOffset = cms.double(0.0),
        TrackEGammaEtaOffset = cms.double(0.0),
	    RelativeIsolation = cms.bool( True ),	# default = True. The isolation variable is relative if True,
						# else absolute.
        EllipticalMatching = cms.bool( False ), #Elliptical Matching
        ClusterpT = cms.bool( False ), #Use Track pT or Cluster pT to calculate corrections
        IsoCut = cms.double( -0.10 ), 		# Cut on the (Trk-based) isolation: only the L1TkEmParticle for which
                                                # the isolation is below RelIsoCut are written into
                                                # the output collection. When RelIsoCut < 0, no cut is applied.
						# When RelativeIsolation = False, IsoCut is in GeV.
        # Determination of the isolation w.r.t. L1Tracks :
        PTMINTRA = cms.double( 2. ),	# in GeV
	    DRmin = cms.double( 0.03),
	    DRmax = cms.double( 0.2 ),
        maxChi2IsoTracks = cms.double(1e10), # max chi2 cut for a track to be considered for isolation computation
        minNStubsIsoTracks = cms.int32(0), # min cut on # of stubs for a track to be considered for isolation computation
	    DeltaZ = cms.double( 0.6 )    # in cm. Used for tracks to be used isolation calculation
)
L1TkIsoElectrons = L1TkElectrons.clone()
L1TkIsoElectrons.IsoCut = cms.double( 0.10 )
# for  LowPt Electron
L1TkElectronsLoose = L1TkElectrons.clone()
L1TkElectronsLoose.TrackEGammaDeltaPhi = cms.vdouble(0.07, 0.0, 0.0)
L1TkElectronsLoose.TrackEGammaDeltaR = cms.vdouble(0.12, 0.0, 0.0)
L1TkElectronsLoose.TrackMinPt = cms.double( 3.0 )


#### Additional collections that right now only the menu team is using - to be renamed/redefined by the EGamma group
# The important change is the EG seed -> PhaseII instead of PhaseI

#barrel
L1TkElectronsCrystal = L1TkElectrons.clone()
L1TkElectronsCrystal.L1EGammaInputTag = cms.InputTag("L1EGammaClusterEmuProducer","L1EGammaCollectionBXVEmulator")
L1TkElectronsCrystal.IsoCut = cms.double(-0.1)

L1TkElectronsCrystalEllipse = L1TkElectronsCrystal.clone()
L1TkElectronsCrystalEllipse.EllipticalMatching = cms.bool(True)
L1TkElectronsCrystalEllipse.ClusterpT = cms.bool(True)
L1TkElectronsCrystalEllipse.TrackEGammaDeltaPhimax = cms.vdouble(0.045, -0.0002)
L1TkElectronsCrystalEllipse.TrackEGammaDeltaEtamax = cms.vdouble(0.018, 1.6667)
L1TkElectronsCrystalEllipse.TrackEGammaPhiOffset = cms.double(0)
L1TkElectronsCrystalEllipse.TrackEGammaEtaOffset = cms.double(0)
L1TkElectronsCrystalEllipse.TrackMinPt = cms.double( 5.0 )# Use 10 GeV cut or 5 GeV cut

L1TkIsoElectronsCrystal = L1TkElectronsCrystal.clone()
L1TkIsoElectronsCrystal.IsoCut = cms.double(0.12)

L1TkIsoElectronsCrystalEllipse = L1TkElectronsCrystalEllipse.clone()
L1TkIsoElectronsCrystalEllipse.IsoCut = cms.double(0.12)

L1TkElectronsLooseCrystal = L1TkElectronsCrystal.clone()
L1TkElectronsLooseCrystal.TrackEGammaDeltaPhi = cms.vdouble(0.07, 0.0, 0.0)
L1TkElectronsLooseCrystal.TrackEGammaDeltaR = cms.vdouble(0.12, 0.0, 0.0)
L1TkElectronsLooseCrystal.TrackMinPt = cms.double( 3.0 )


#endcap
L1TkElectronsHGC=L1TkElectrons.clone()
L1TkElectronsHGC.L1EGammaInputTag = cms.InputTag("l1EGammaEEProducer","L1EGammaCollectionBXVWithCuts")
L1TkElectronsHGC.IsoCut = cms.double(-0.1)

L1TkElectronsHGCEllipse=L1TkElectronsHGC.clone()
L1TkElectronsHGCEllipse.EllipticalMatching = cms.bool(True)
L1TkElectronsHGCEllipse.ClusterpT = cms.bool(True)
L1TkElectronsHGCEllipse.TrackEGammaDeltaPhimax = cms.vdouble(0.025, 0)
L1TkElectronsHGCEllipse.TrackEGammaDeltaEtamax = cms.vdouble(0.01, 1)
L1TkElectronsHGCEllipse.TrackEGammaPhiOffset = cms.double(0)
L1TkElectronsHGCEllipse.TrackEGammaEtaOffset = cms.double(0)
L1TkElectronsHGCEllipse.TrackMinPt = cms.double( 5.0 )

L1TkIsoElectronsHGC=L1TkElectronsHGC.clone()
L1TkIsoElectronsHGC.DRmax = cms.double(0.4)
L1TkIsoElectronsHGC.DeltaZ = cms.double(1.0)
L1TkIsoElectronsHGC.maxChi2IsoTracks = cms.double(100)
L1TkIsoElectronsHGC.minNStubsIsoTracks = cms.int32(4)
L1TkIsoElectronsHGC.IsoCut = cms.double(0.30)

L1TkIsoElectronsHGCEllipse=L1TkElectronsHGCEllipse.clone()
L1TkIsoElectronsHGCEllipse.DRmax = cms.double(0.4)
L1TkIsoElectronsHGCEllipse.DeltaZ = cms.double(1.0)
L1TkIsoElectronsHGCEllipse.maxChi2IsoTracks = cms.double(100)
L1TkIsoElectronsHGCEllipse.minNStubsIsoTracks = cms.int32(4)
L1TkIsoElectronsHGCEllipse.IsoCut = cms.double(0.30)

L1TkElectronsLooseHGC = L1TkElectronsHGC.clone()
L1TkElectronsLooseHGC.TrackEGammaDeltaPhi = cms.vdouble(0.07, 0.0, 0.0)
L1TkElectronsLooseHGC.TrackEGammaDeltaR = cms.vdouble(0.12, 0.0, 0.0)
L1TkElectronsLooseHGC.TrackMinPt = cms.double( 3.0 )


