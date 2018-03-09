import ROOT as r
import os.path
r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L L1TrackNtuplePlot.C++");

treeName = '_TMTT_KF4ParamsComb'
truncationOptions = [
'WithTruncation',
'WithoutTruncation'
]

PUs = [
'0',
'140',
'200'
]

pGunSamples = {
	'MuonPt10' : 13,
	'MuonPt100' : 13,
	'ElectronPt10' : 11,
	'ElectronPt35' : 11,
}

run = r.gROOT.ProcessLine


print 'Running ttbar samples'
for PU in PUs:
	for truncation in truncationOptions:
		inputFile = 'TTbar_{PU}_{truncation}'.format( PU=PU, truncation = truncation)
		if not os.path.isfile( inputFile+'.root' ): 
			continue
		print '---> ',inputFile
 		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,0,0,3)'.format(inputFile=inputFile, treeName=treeName) )
 		print '---------> Muons'
 		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,13,0,3)'.format(inputFile=inputFile, treeName=treeName) )
 		print '---------> Electrons'
 		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,11,0,3)'.format(inputFile=inputFile, treeName=treeName) )


		print 'In jet'
		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",1,0,0,3)'.format(inputFile=inputFile, treeName=treeName) )
		print 'In high pt jet'
		run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",2,0,0,3)'.format(inputFile=inputFile, treeName=treeName) )

		print "Done"

print 'PGun samples'

for sample, pdg in pGunSamples.iteritems():
	print '---> ',sample
	for PU in PUs:
		for truncation in truncationOptions:
			inputFile = '{sample}_{PU}_{truncation}'.format( sample=sample, PU=PU, truncation = truncation)
			if not os.path.isfile( inputFile+'.root' ): 
				continue
			print '--------->',inputFile
 			run ( 'L1TrackNtuplePlot("{inputFile}","{treeName}",0,{pdg},0,3)'.format(inputFile=inputFile,pdg=pdg, treeName=treeName) )

