import ROOT
import cmsstyle as CMS

#Useful functions
def AddHist(histlist):
    h_add = histlist[0].Clone('h_add')
    for h in histlist[1:]:
        h_add.Add(h)
    return h_add

#------------------------------------#
#Read the file from the disk
inputDir = "input/"    
file_DY    = ROOT.TFile.Open(inputDir + "ttbar3LCR2018_DY_MC_May13_v5.root","READ")
file_ttbar = ROOT.TFile.Open(inputDir + "ttbar3LCR2018_ttbar_MC_May13_v5.root","READ")
file_ZZ    = ROOT.TFile.Open(inputDir + "ttbar3LCR2018_ZZ_MC_May13_v5.root","READ")
file_WZ    = ROOT.TFile.Open(inputDir + "ttbar3LCR2018_WZ_MC_May13_v5.root","READ")
file_TTZ   = ROOT.TFile.Open(inputDir + "ttbar3LCR2018_TTZ_MC_May13_v5.root","READ")
file_TTW   = ROOT.TFile.Open(inputDir + "ttbar3LCR2018_TTW_MC_May13_v5.root","READ")
file_Data  = ROOT.TFile.Open(inputDir + "ttbar3LCR2018_Data_May13_v5.root","READ")

print(" \nfile opened in ROOT successfully..")


## Get the histograms

plotname = "mass12_all"
h_dy   = file_DY.Get(plotname); 
h_tt   = file_ttbar.Get(plotname);
h_zz   = file_ZZ.Get(plotname);
h_wz   = file_WZ.Get(plotname);  
h_ttz  = file_TTZ.Get(plotname);  
h_ttw  = file_TTW.Get(plotname); 
h_data = file_Data.Get(plotname);

##Scale
dtlumi = 59.8*1000;
ttlumi = 28701360/88.29;
ttzlumi = 13280000/0.2432;
ttwlumi = 4911941/0.2149;
dylumi = 99717900/5765.0;  
wzlumi =10527550/5.052;
zzlumi =98613000/1.325;

h_dy.Scale(dtlumi/dylumi)
h_tt.Scale(dtlumi/ttlumi)
h_wz.Scale(dtlumi/wzlumi)
h_zz.Scale(dtlumi/zzlumi)
h_ttw.Scale(dtlumi/ttwlumi)
h_ttz.Scale(dtlumi/ttzlumi)

print(f"Histograms read from file: {plotname} ")

####
rebin = 50

MChist={}
MChist['DY']=h_dy.Clone().Rebin(rebin)
MChist['TT']=h_tt.Clone().Rebin(rebin)
MChist['ZZ']=h_zz.Clone().Rebin(rebin)
MChist['WZ']=h_wz.Clone().Rebin(rebin)
MChist['TTW']=h_ttw.Clone().Rebin(rebin)
MChist['TTZ']=h_ttz.Clone().Rebin(rebin)

data = h_data.Clone('data')
data.Rebin(rebin)

bkgh=[]
for bkg,h in MChist.items():
    bkgh.append(h)

totbkg = AddHist(bkgh)

MChist = dict(sorted(MChist.items(), key=lambda item: item[1].Integral()))

import math
#Styling
CMS.SetLumi(7.98)
CMS.SetEnergy("13.6")
CMS.SetExtraText("Work in progress")
square=CMS.kSquare
iPos=0

# Create canvas
canv_name = "stack_ratio"
dicanv = CMS.cmsDiCanvas(canv_name, 0, 1000, 0.1, 1E5, 0, 2.0, plotname, "Events", "", square=square, extraSpace=0.1, iPos=iPos,)

# Upper pad
dicanv.cd(1)
ROOT.gPad.SetLogy()
stack = ROOT.THStack("stack", "Stacked")
h_err = totbkg.Clone("h_err")
leg = CMS.cmsLeg(0.60, 0.89 - 0.05 * 5, 0.99, 0.89, textSize=0.04,columns=2)
CMS.cmsDrawStack(stack, leg, MChist)
#leg.AddEntry(h1.Clone(), "h1", "l")
leg.AddEntry(data.Clone(), "Data", "pe")
CMS.cmsDraw(h_err, "e2same0", lcolor = 335, lwidth = 1, msize = 0, fcolor = ROOT.kBlack, fstyle = 3004,) 
CMS.cmsDraw(data, "E1X0", mcolor=ROOT.kBlack)
#h1.Draw()
# Lower pad

dicanv.cd(2)
leg_ratio = CMS.cmsLeg(
    0.17, 0.97 - 0.05 * 5, 0.35, 0.97, textSize=0.05, columns=2
)
ratio = data.Clone("ratio")
ratio.Divide(totbkg.Clone())
for i in range(1,ratio.GetNbinsX()+1):
    if(ratio.GetBinContent(i)):
        ratio.SetBinError(i, math.sqrt(data.GetBinContent(i))/data.Clone().GetBinContent(i))
    else:
        ratio.SetBinError(i, 10^(-99))
        
yerr_root = ROOT.TGraphAsymmErrors()
yerr_root.Divide(data.Clone(), totbkg.Clone(),"pois")
for i in range(0,yerr_root.GetN()+1):
    yerr_root.SetPointY(i,1)
CMS.cmsDraw(yerr_root, "e2same0", lwidth = 100, msize = 0, fcolor = ROOT.kBlack, fstyle = 3004)  
CMS.cmsDraw(ratio, "E1X0", mcolor=ROOT.kBlack)
ref_line = ROOT.TLine(-2, 1, 2, 1)
CMS.cmsDrawLine(ref_line, lcolor=ROOT.kBlack, lstyle=ROOT.kDotted)
dicanv.Draw()
dicanv.SaveAs("cmsstyle_stackplot.pdf")
