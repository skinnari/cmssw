# Recipe for producing template slides

## First time setup
```
git clone git@github.com:EmyrClement/beamer-plot-slides.git
cd beamer-plot-slides/
git checkout templateSlides
cd ../
```

## Run L1TrackNtuplePlot.C for all samples and options
```
python runAllForTemplateSlides.py
```

## Produce overlay plots
```
python compareEfficiency.py 
python compareResolution.py
```

## Produce slides
```
cd beamer-plot-slides/
./make_slides.py ../template.json
cd ../
ls template_slides.pdf
```
## Samples
These are the samples for the template slides

```
/RelValTTbar_14TeV/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/MINIAODSIM
/RelValTTbar_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140-v1/MINIAODSIM
/RelValTTbar_14TeV/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/MINIAODSIM

/RelValSingleElectronPt10Extended/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/MINIAODSIM
/RelValSingleElectronPt10Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140-v1/MINIAODSIM
/RelValSingleElectronPt10Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/MINIAODSIM

/RelValSingleElectronPt35Extended/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/MINIAODSIM
/RelValSingleElectronPt35Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140-v1/MINIAODSIM
/RelValSingleElectronPt35Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/MINIAODSIM

/RelValSingleMuPt10Extended/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/MINIAODSIM
/RelValSingleMuPt10Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140-v1/MINIAODSIM
/RelValSingleMuPt10Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/MINIAODSIM

/RelValSingleMuPt100Extended/CMSSW_9_3_2-93X_upgrade2023_realistic_v2_2023D17noPU-v1/MINIAODSIM
/RelValSingleMuPt100Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140-v1/MINIAODSIM
/RelValSingleMuPt100Extended/CMSSW_9_3_2-PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200-v1/MINIAODSIM
```
