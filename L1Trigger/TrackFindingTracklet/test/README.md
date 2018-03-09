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
