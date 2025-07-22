#class implementation of datacard maker
#Arnab Laha
#
#Shape datacard

import ROOT
import os,sys
class DatacardMaker(object):
    
    def __setattr__(self,key,value):
        if not hasattr(self,key):
            raise TypeError("%r is not a valid pset parameter" % key)
        object.__setattr__(self,key,value)
    
    def __init__(self):
        pass

    channel  = ''
    year     = ''
    variable = ''
    signal   = ''
    processlist=[]    
    inputfile= ''
    statonly=True
    hprefix=None
    outputdir=os.getcwd()
    syst_table=[]
    debug = 0
    
    def setChannel(self,channel):
        self.channel = channel

    def setYear(self,year):
        self.year = year

    def setVariable(self,variable):
        self.variable = variable

    def setSignal(self,signal):
        self.signal=signal
        self.processlist.append(signal)
        
    def setBackground(self,bkglist):
        self.processlist.extend(bkglist)
        
    def setData(self,data):
        self.processlist.append(data)
        
    def setInputFile(self,inputfile):
        self.inputfile = inputfile

    def setStatOnly(self,statonly):
        self.statonly=statonly
        
    def setHistogramPrefix(self,hprefix):
        self.hprefix = hprefix
        
    def addSystematics(self,name,model,variation,process):
        self.syst_table.append((name,model,variation,process))
        

    def show(self):
        print('datacard settings')
        print('channel  :',self.channel)
        print('year     :',self.year)
        print('variable :',self.variable)
        print('inputfile:',self.inputfile)
        print('processes:',self.processlist)
        print('signal   :',self.signal)
        print('StatOnly :',self.statonly)
        print('outputdir:',self.outputdir)
        if not self.statonly:            
            print('\n\033[32msyst table\033[0m\n')
            for name,model,variation,sysprocess in self.syst_table:
                if '*' in sysprocess:sysprocess = self.processlist.copy()
                if 'data_obs' in sysprocess:sysprocess.remove('data_obs')
                print(name,model,variation,sysprocess)
        print()
        
    def make(self):
        print('datacard maker::begin')
        ###############
        ## open file ##
        ###############
        rfile = ROOT.TFile.Open(self.inputfile,'READ')
        if(self.debug>1):
            for keys in rfile.GetListOfKeys():
                print(keys)
                
        ##set histogram prefix
        if self.hprefix==None:
            self.hprefix=f'ch{self.channel}{self.year}_{self.variable}'
            print()
            print("No histogram prefix (hprefix) chosen")
            print("Nomenclature::")
            print("histogram prefix: ch<channel><year>_<variable>")
            print("histogram name  : <hprefix>_<process>_<syst>")
            print("Use setHistogramPrefix() to set externally")            
            print("histogram prefix set to: ",self.hprefix)
            print()
        else:
            print("histogram prefix: ",self.hprefix)
            print()

        ##fetch nominal histograms    
        histodict={}
        for process in self.processlist:
            histodict[process]=rfile.Get(self.hprefix+'_'+process)

        print("\nprocess/hname nominal/yield-integral")
        for k,v in histodict.items():
            print(f"{k}/{v.GetName()}/{v.Integral()}")

        ##directories
        main_dir     = os.path.join(self.outputdir,f"SHAPES_{self.signal}_{self.variable}")
        datacard_dir = os.path.join(main_dir,f"SHAPES_{self.signal}_ch{self.channel}{self.year}_{self.variable}")
        if not os.path.isdir(main_dir):os.makedirs(main_dir)
        if not os.path.isdir(datacard_dir):os.makedirs(datacard_dir)
        print("\ncreated directories to save datacards")
        print(f"main dir    : {main_dir} :: save yearwise datacard directories")
        print(f"datacard dir: {datacard_dir} :: save datacards + combined shape histograms")
        print()
        #################
        
        print("\ncreating datacard text file\n")
        if(self.statonly):
            dcfilename=f"{datacard_dir}/datacardStatOnly_{self.signal}_ch{self.channel}{self.year}_{self.variable}.txt"
        else:
            dcfilename=f"{datacard_dir}/datacard_{self.signal}_ch{self.channel}{self.year}_{self.variable}.txt"            
            
        ##datacard content
        bundles=self.processlist.copy()
        bundles.remove('data_obs')        
        bundles.remove(self.signal)
        bundles= [self.signal]+bundles
        print("\nprocesses (without data_obs): ",bundles)
        
        if not self.inputfile.startswith("SHAPES"):
            rfilename = self.inputfile.split('/')[-1]
            print("SHAPE histogram filename: ",rfilename)
            print()
        else:
            rfilename = self.inputfile

        #action starts   
        dcfile = open(dcfilename,'w')
        sys.stdout = dcfile
        print('imax * number of channels')
        print('jmax * number of backgrounds')
        print('kmax * number of nuisance parameters (sources of systematical uncertainties)')

        print(f'shapes * * {rfilename} $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC')
        print()
        print()
        print('bin          ',end ='')
        print(self.hprefix, end=' ')
        print()
        print('observation  ',end='')
        print(str(histodict['data_obs'].Integral())+'     ',end='')
        print()
        print()
        print('bin      ',end='')
        for process in range(len(bundles)):print(self.hprefix+'  ',end='')
        print()
        print('process  ',end='')
        for process in bundles: print(process+'   ',end='')
        print()
        print('process  ',end='')
        for iprocess in range(len(bundles)): print(str(iprocess)+'     ',end='')
        print()
        print('rate     ',end='')
        for process in bundles: print(str(histodict[process].Integral())+' ',end='')
        print()
        print("----------------------------------------------------------------------------------------------")
        if not (self.statonly):
            for systname,systmodel,systvalue,systproclist in self.syst_table:            
                print(f'{systname} {systmodel} ',end='')
                for iprocess in range(len(bundles)):
                    if '*' in systproclist or bundles[iprocess] in systproclist:
                        print(str(systvalue)+' ',end='')
                    else:
                        print('-'+' ',end='')
                print()
            print()
        ##stat-error
        print()
        print('* autoMCStats 0 0 1', end='')
        print()
        print()
        dcfile.close()
        sys.stdout=sys.__stdout__
        
        os.system(f'cp {self.inputfile} {datacard_dir}')
        print(f'datacard file    >>>',dcfilename)
        print(f'shape histograms >>> {datacard_dir}/{rfilename}')
        print()
        print("datacard maker::end")
        
        
