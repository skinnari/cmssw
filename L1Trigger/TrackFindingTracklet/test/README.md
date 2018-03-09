# Recipe for producing template slides


## Run L1TrackNtuplePlot.C for all samples and options
```
python runAllForTemplateSlides.py
```

### Produce overlay plots
```
python compareEfficiency.py 
python compareResolution.py
```

### Produce slides
```
cd beamer-plot-slides/
./make_slides.py ../template.json
cd ../
ls template_slides.pdf
```
