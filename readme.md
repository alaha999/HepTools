### Example Shape datacard with dummy systematics

1) creating dummy shape histograms:`example_makeHistograms.ipynb`
2) datacard maker module (class based implementation): `datacardMaker.py`
3) example datacard producer notebook: `datacardProducer.ipynb`


>>> input is : combined shape histograms
>>> output is: datacards

A few examples,

```txt
imax * number of channels
jmax * number of backgrounds
kmax * number of nuisance parameters (sources of systematical uncertainties)
shapes * * SHAPES_VLLM400_ch3L2018_LTMET.root $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC


bin          ch3L2018_LTMET 
observation  52.0     

bin      ch3L2018_LTMET  ch3L2018_LTMET  ch3L2018_LTMET  ch3L2018_LTMET  
process  VLLM400   DY   TT   WZ   
process  0     1     2     3     
rate     11.0 18.0 9.0 24.0 
----------------------------------------------------------------------------------------------
lumi lnN 1.025 1.025 1.025 1.025 
WZnorm lnN - - - 1.2 
CMS_btag shape 1 - 1 - 
CMS_JES shape - 1 1 - 
CMS_POGID shape 1 1 1 1 


* autoMCStats 0 0 1

```

### Example datacard (stat Only)

```
imax * number of channels
jmax * number of backgrounds
kmax * number of nuisance parameters (sources of systematical uncertainties)
shapes * * SHAPES_VLLM400_ch3L2018_LTMET.root $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC


bin          ch3L2018_LTMET 
observation  52.0     

bin      ch3L2018_LTMET  ch3L2018_LTMET  ch3L2018_LTMET  ch3L2018_LTMET  
process  VLLM400   DY   TT   WZ   
process  0     1     2     3     
rate     11.0 18.0 9.0 24.0 
----------------------------------------------------------------------------------------------

* autoMCStats 0 0 1

```


