import ROOT
ROOT.EnableImplicitMT()

fname= "../VLL2018_Test_VLL_M150_sample.root"

df = ROOT.RDataFrame("events_1L2J",fname)

print("No of events =",df.Count().GetValue())

df = df.Define("weight","1.0")

nominal_df = df.Vary("weight","ROOT::RVecD{weight*0.9, weight*1.1}",["down","up"])\
               .Histo1D(("event_ST","ST",1000,0,1000),"event_ST","weight")


histo = ROOT.RDF.Experimental.VariationsFor(nominal_df)
print(histo.GetKeys())


##---------------------------------------------------------------------------------------------------
def SetOverflowBin(histo):
    nbins = histo.GetNbinsX()
    histo.SetBinContent(nbins, histo.GetBinContent(nbins) + histo.GetBinContent(nbins+1)); ## Overflow
    histo.SetBinContent(1, histo.GetBinContent(1)+ histo.GetBinContent(0));                ## Underflow

##SetOverflow
SetOverflowBin(histo['nominal'])
SetOverflowBin(histo['weight:up'])
SetOverflowBin(histo['weight:down'])
##Rebin
rebin = 25
histo["nominal"].Rebin(rebin)
histo["weight:up"].Rebin(rebin)
histo["weight:down"].Rebin(rebin)
##SetColor
histo["weight:down"].SetLineColor(ROOT.kRed)
histo["weight:up"].SetLineColor(ROOT.kGreen+1)

##ratioHisto
ratioHist_up = histo["weight:up"].Clone()
ratioHist_up.Divide(histo["nominal"])
ratioHist_down = histo["weight:down"].Clone()
ratioHist_down.Divide(histo["nominal"])

#Legend
legend = ROOT.TLegend(0.95,0.50,0.80,0.86)
legend.AddEntry(histo["nominal"],"nominal","l")
legend.AddEntry(histo["weight:up"],"up","l")
legend.AddEntry(histo["weight:down"],"down","l")

##plotting
def PadStyling(pad,rpad):
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.20)
    pad.SetTopMargin(0.09)
    pad.SetBottomMargin(0.01)
    pad.SetTickx(1)
    pad.SetTicky(1)

    rpad.SetLeftMargin(0.15)
    rpad.SetRightMargin(0.20)
    rpad.SetTopMargin(0.02)
    rpad.SetBottomMargin(0.40)
    rpad.SetTickx(1)
    rpad.SetTicky(1)
    rpad.SetGrid(1)

def SetLegendStyle(legend):
    legend.SetTextFont(62)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.024)


##StartPlotting
canvas = ROOT.TCanvas("c1","systvar",650,600)
ROOT.gStyle.SetOptStat(0)

ratioPadSize = 0.3
mainPad  = ROOT.TPad("pad","pad",0,ratioPadSize,1,1)
ratioPad = ROOT.TPad("pad2","pad2",0,0,1.0,ratioPadSize)
PadStyling(mainPad,ratioPad)
mainPad.Draw()
ratioPad.Draw()

mainPad.cd()
mainPad.SetLogy(1)
histo["nominal"].Draw("HIST")
histo["weight:down"].Draw("HIST SAME")
histo["weight:up"].Draw("HIST SAME")
histo["nominal"].GetXaxis().SetLabelSize(1)
SetLegendStyle(legend)
legend.Draw()
mainPad.SetTickx(1)

ratioPad.cd()
ratioHist_up.Draw("ep")
ratioHist_down.Draw("epsame")
ratioHist_up.GetYaxis().SetRangeUser(0,2)
ratioHist_up.GetYaxis().SetLabelFont(43)
ratioHist_up.GetYaxis().SetLabelSize(12)
ratioHist_up.SetTitle('')
ratioHist_up.GetYaxis().SetNdivisions(503)
ratioHist_up.GetXaxis().SetNdivisions(513)
ratioHist_up.GetXaxis().SetTitle("event_ST")
ratioHist_up.GetXaxis().CenterTitle()
ratioHist_up.GetXaxis().SetTitleFont(43)
ratioHist_up.GetXaxis().SetTitleSize(20)
ratioHist_up.GetXaxis().SetTitleOffset(1.2)
ratioHist_up.GetXaxis().SetLabelFont(43)
ratioHist_up.GetXaxis().SetLabelSize(12)

canvas.Draw()
canvas.SaveAs("systvar.pdf")
