### Example Shape datacard with dummy systematics

1) creating dummy shape histograms:`example_makeHistograms.ipynb`
2) datacard maker module (class based implementation): `datacardMaker.py`
3) example datacard producer notebook: `datacardProducer.ipynb`


- Input : `combined shape histograms`
- Output: `datacards`

Content of the input shape histograms should be:

```bash
$ root shape_histograms/SHAPES_VLLM400_ch3L2018_LTMET.root
$ _file0->ls()

Name: ch3L2018_LTMET_VLLM400 Title: VLLM400_LTMET
Name: ch3L2018_LTMET_VLLM400_CMS_btagUp Title: VLLM400_LTMET
Name: ch3L2018_LTMET_VLLM400_CMS_btagDown Title: VLLM400_LTMET
Name: ch3L2018_LTMET_VLLM400_CMS_JESUp Title: VLLM400_LTMET
Name: ch3L2018_LTMET_VLLM400_CMS_JESDown Title: VLLM400_LTMET
Name: ch3L2018_LTMET_VLLM400_CMS_POGIDUp Title: VLLM400_LTMET
Name: ch3L2018_LTMET_VLLM400_CMS_POGIDDown Title: VLLM400_LTMET
Name: ch3L2018_LTMET_DY Title: DY_LTMET
Name: ch3L2018_LTMET_DY_CMS_btagUp Title: DY_LTMET
Name: ch3L2018_LTMET_DY_CMS_btagDown Title: DY_LTMET
Name: ch3L2018_LTMET_DY_CMS_JESUp Title: DY_LTMET
Name: ch3L2018_LTMET_DY_CMS_JESDown Title: DY_LTMET
Name: ch3L2018_LTMET_DY_CMS_POGIDUp Title: DY_LTMET
Name: ch3L2018_LTMET_DY_CMS_POGIDDown Title: DY_LTMET
Name: ch3L2018_LTMET_TT Title: TT_LTMET
Name: ch3L2018_LTMET_TT_CMS_btagUp Title: TT_LTMET
Name: ch3L2018_LTMET_TT_CMS_btagDown Title: TT_LTMET
Name: ch3L2018_LTMET_TT_CMS_JESUp Title: TT_LTMET
Name: ch3L2018_LTMET_TT_CMS_JESDown Title: TT_LTMET
Name: ch3L2018_LTMET_TT_CMS_POGIDUp Title: TT_LTMET
Name: ch3L2018_LTMET_TT_CMS_POGIDDown Title: TT_LTMET
Name: ch3L2018_LTMET_WZ Title: WZ_LTMET
Name: ch3L2018_LTMET_WZ_CMS_btagUp Title: WZ_LTMET
Name: ch3L2018_LTMET_WZ_CMS_btagDown Title: WZ_LTMET
Name: ch3L2018_LTMET_WZ_CMS_JESUp Title: WZ_LTMET
Name: ch3L2018_LTMET_WZ_CMS_JESDown Title: WZ_LTMET
Name: ch3L2018_LTMET_WZ_CMS_POGIDUp Title: WZ_LTMET
Name: ch3L2018_LTMET_WZ_CMS_POGIDDown Title: WZ_LTMET
Name: ch3L2018_LTMET_data_obs Title: data_obs_LTMET
```

# Rules and nomenclatures

### Naming nomenclature:
```
For each process: Shapes_<process>_ch<channel><year>_<variable>.root
Combined shape histograms: SHAPES_<signal>_ch<channel><year>_<variable>.root

NB: for combined shape histograms <signal> is important, as it distinguishes from other signal mass points.

Example: SHAPES_VLLtauM400_ch3L2018_LTMET.root
It means the shape histograms are for:
 - VLLtauM400 signal, channel=3L, year= 2018, discriminating variable = LTMET

This way of naming would be easy to combine channel wise, yearwise later at datacard level
```
### Systematics Naming nomenclature
```
- nominal: ch<channel><year>_<variable> 
- Syst Up: ch<channel><year>_<variable>_<SystName>Up
- Syst Down: ch<channel><year>_<variable>_<SystName>Down

example:
- nominal: ch3L2018_LTMET_DY
- CMS_btagUp: ch3L2018_LTMET_DY_CMS_btagUp
- CMS_btagDown: ch3L2018_LTMET_DY_CMS_btagDown

*data histograms should be always named as: data_obs*
*for other process: use some sensible names, not too long, not too confusing*

See the example_makeHistograms.ipynb notebook for an example setup 

```

### Example datacard (with systematics)

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