import ROOT
from optparse import OptionParser, OptionGroup, OptionValueError
import argparse
import os
import sys
import yaml
from array import array
import math

RF = ROOT.RooFit

ROOT.gROOT.LoadMacro("/afs/cern.ch/work/a/altaylor/public/ATLASStyle/AtlasStyle.C")
ROOT.gROOT.LoadMacro("/afs/cern.ch/work/a/altaylor/public/ATLASStyle/AtlasLabels.C")
ROOT.gROOT.LoadMacro("/afs/cern.ch/work/a/altaylor/public/ATLASStyle/AtlasUtils.C")
ROOT.gSystem.Load('RooTwoSidedCBShape/RooTwoSidedCBShape_cxx')

from ROOT import RooTwoSidedCBShape

## suppress error messages below warning
ROOT.RooMsgService.instance().setGlobalKillBelow(RF.WARNING) 

ROOT.SetAtlasStyle()
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# use a custom palette for the residual study 
NRGBs = 5
NCont = 255
stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
red = [ 0.00, 0.00, 0.87, 1.00, 0.51 ]
green = [ 0.00, 0.81, 1.00, 0.20, 0.00 ]
blue = [ 0.51, 1.00, 0.12, 0.00, 0.00 ]
stopsArray = array('d', stops)
redArray   = array('d', red)
greenArray = array('d', green)
blueArray  = array('d', blue)
ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, NCont)
ROOT.gStyle.SetNumberContours(NCont)

## decide which functions to test, either E1 or E2 for just now but could be extended in the YAML file
parser = argparse.ArgumentParser()

parser.add_argument('-myy', action='store', required=True, help='Function choice for myy')
parser.add_argument('-mjj', action='store', required=True, help='Function choice for mjj')
parser.add_argument('-catName', action='store', required=True, help='Give the name of the category')
parser.add_argument('-yamlFile', action='store', required=True, help='Give the name of the YAML file')
args = parser.parse_args()


def PlotSubtraction(obsVariable, nBins,subtractMax, originHisto):
    
    medianHist = ROOT.TH1F("median_"+combWS.var(obsVariable).GetName(), "median_"+combWS.var(obsVariable).GetName(), nBins, combWS.var(obsVariable).getMin(), combWS.var(obsVariable).getMax() )
    result = ROOT.TGraphErrors()

    for x in range(1,nBins+1):
        medianHist.SetBinContent(x,0.0)

    medianHist.SetLineColor(ROOT.kRed)
    medianHist.SetLineWidth(2)
    medianHist.GetXaxis().SetTitle(combWS.var(obsVariable).GetName()+" [GeV]")
    medianHist.GetYaxis().SetTitle("MC - Fit")

    medianHist.GetXaxis().SetTitleOffset(0.95)
    medianHist.GetYaxis().SetTitleOffset(0.7)
    medianHist.GetXaxis().SetTitleSize(0.1)
    medianHist.GetYaxis().SetTitleSize(0.1)
    medianHist.GetXaxis().SetLabelSize(0.1)
    medianHist.GetYaxis().SetLabelSize(0.1)
    medianHist.GetYaxis().SetRangeUser(-subtractMax, +subtractMax)

    pdfName = "bkg_"+combWS.var(obsVariable).GetName()

    line = ROOT.TLine()
    line.SetLineStyle(1)
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(1)
    line.SetLineStyle(2)

    m_ratioMin = -2.0
    m_ratioMaxi = 2.0

    line.DrawLine(m_ratioMin, 0.0, m_ratioMaxi, 0.0)

    n_events = data.sumEntries(obsVariable + " > " + str(combWS.var(obsVariable).getMin()) +  " && " + obsVariable + " < " + str(combWS.var(obsVariable).getMax()) )
    result = ROOT.TGraphErrors()
    increment = ( combWS.var(obsVariable).getMax() - combWS.var(obsVariable).getMin() ) / nBins
    combWS.var(obsVariable).setRange("fullRange", combWS.var(obsVariable).getMin(), combWS.var(obsVariable).getMax())

    int_tot = combWS.pdf(pdfName).createIntegral(  ROOT.RooArgSet( combWS.var(obsVariable) ), RF.NormSet( ROOT.RooArgSet( combWS.var(obsVariable) ) ), RF.Range("fullRange")  )
    val_tot = int_tot.getVal()

    pointIndex = 0
    pointIndexNonZero = 0

    result = ROOT.TGraphErrors()
    i_m = combWS.var(obsVariable).getMin()

    # need original histogram for the data - fit plot
    while i_m < combWS.var(obsVariable).getMax():

        combWS.var(obsVariable).setRange("range_"+ str(i_m), i_m, ( i_m + increment ) )

        int_curr = combWS.pdf(pdfName).createIntegral( ROOT.RooArgSet( combWS.var(obsVariable) ), RF.NormSet( ROOT.RooArgSet( combWS.var(obsVariable) ) ), RF.Range("range_"+ str(i_m))  )

        val_curr = int_curr.getVal()
        currMass = i_m + 0.5*increment
        # compute the central value of the PDF in that bin 
        curr_pdf_weight = n_events * ( val_curr / val_tot )
        var_name = combWS.var(obsVariable).GetName()
        curr_data_weight = data.sumEntries( obsVariable + " > " + str(i_m) + " && " + obsVariable +" < " + str(i_m + increment) )
        # do the subtraction 
        curr_weight = curr_data_weight - curr_pdf_weight
        result.SetPoint(pointIndex,currMass, curr_weight)
        curr_error = originHisto.GetBinError(pointIndex+1)
        result.SetPointError(pointIndex, 0.0, curr_error)
        
        pointIndex += 1
        i_m += increment

    histoAndTGraph = [ ]
    histoAndTGraph.append( medianHist )
    histoAndTGraph.append( result )

    return histoAndTGraph


## create a directory for the results
if not os.path.exists("Plots"):
    print 'making Plots directory...'
    os.makedirs("Plots")

## and create a sub-directory for the category studied..

if not os.path.exists("Plots/"+str(args.catName)):
    print 'making directory for category : ', args.catName
    os.makedirs("Plots/"+str(args.catName))


## open YAML file which has the configuration
configFile = yaml.safe_load(open(str(args.yamlFile)))

## get the background template
inputFile = ROOT.TFile(configFile["fileName"])
inputHisto = inputFile.Get(configFile["histoName"])

## check there is a consistency between the number of bins in the YAML file and in the histogram provided.

if inputHisto.ProjectionX().GetXaxis().GetNbins() == configFile["n_bins_yy"]:
    if inputHisto.ProjectionY().GetXaxis().GetNbins() == configFile["n_bins_jj"]:
        print('Good, number of bins in the provided histogram match the number of bins in the YAML file')
    else:
        sys.exit('Exiting, number of bins in histogram and number of bins in the YAML file dont match!' )





## get the signal model for mjj

mjjModelFile = ROOT.TFile(configFile["mjjWS"])
mjjWS = mjjModelFile.Get(configFile["mjjWSname"])
mjjPDF = mjjWS.pdf(configFile["mjjPDF"])

## get the signal model for myy

myyModelFile = ROOT.TFile(configFile["myyWS"])
myyWS = myyModelFile.Get(configFile["myyWSname"])
myyPDF = myyWS.pdf(configFile["myyPDF"])

## set up the fit
# 1 gev bins for myy
myy_min = configFile["myy_min"]
myy_max = configFile["myy_max"]
n_bins_yy = configFile["n_bins_yy"]

# 5 gev bins for mbb
mjj_min = configFile["mjj_min"]
mjj_max = configFile["mjj_max"]
n_bins_jj = configFile["n_bins_jj"]

varX = "myy"
varY = "mjj"

myy = ROOT.RooRealVar(varX, varX, myy_min, myy_max)
varset = ROOT.RooArgSet(myy)

mjj = ROOT.RooRealVar(varY, varY, mjj_min, mjj_max)
varset_mjj = ROOT.RooArgSet(mjj)

data = ROOT.RooDataHist('data', 'data', ROOT.RooArgList(myy,mjj), inputHisto)

## create workspace to do the fit 
combWS = ROOT.RooWorkspace("combWS","combWS")

getattr(combWS,'import')(data)
getattr(combWS,'import')(myyPDF)
## add suffix onto all variables with the except of mjj, the dataset name
getattr(combWS,'import')(mjjPDF, RF.RenameAllNodes(varY), RF.RenameAllVariablesExcept(varY, varY))

## set all the parameters in the diphoton and di b-jet PDFs to constant. 
myy_params = myyPDF.getParameters(ROOT.RooArgSet())
iter = myy_params.createIterator()
myy_variable = iter.Next()
while myy_variable:
    combWS.var(myy_variable.GetName()).setConstant(True)
    myy_variable = iter.Next()

mjj_params = mjjPDF.getParameters(ROOT.RooArgSet())
iter = mjj_params.createIterator()
mjj_variable = iter.Next()
while mjj_variable:
    if mjj_variable.GetName() != varY:
        combWS.var(mjj_variable.GetName()+"_"+varY).setConstant(True)
    mjj_variable = iter.Next()

combWS.factory("PROD::signalModel("+configFile["myyPDF"]+", "+configFile["mjjPDF"]+"_"+varY+")")

## background model
backgroundFuncMyy = configFile["bkgMyyPDF"][args.myy]
backgroundFuncMjj = configFile["bkgMjjPDF"][args.mjj]

combWS.factory(backgroundFuncMyy)
combWS.factory(backgroundFuncMjj)
combWS.factory("PROD::bkg_only(bkg_myy,bkg_mjj)")

## do a background only fit 

pdfname = "bkg_only"
bkg_model = combWS.pdf(pdfname)
bkg_model.fitTo(data,RF.Save(), RF.SumW2Error(True), RF.Offset(True), RF.Minimizer("Minuit2", "migrad"))

## do some background visuals
c1 = ROOT.TCanvas("can", "can", 800, 800)
pad1 = ROOT.TPad( "pad1", "pad1", 0.00, 0.33, 1.00, 1.00 )
pad2 = ROOT.TPad( "pad2", "pad2", 0.00, 0.00, 1.00, 0.33 )
pad1.SetBottomMargin(0.00001)
pad1.SetBorderMode(0)
pad2.SetTopMargin(0.00001)
pad2.SetBottomMargin(0.4)
pad2.SetBorderMode(0)
c1.cd()
pad1.Draw()
pad2.Draw()
pad1.cd()

x_frame = myy.frame(myy_min , myy_max , n_bins_yy)
data.plotOn(x_frame, RF.XErrorSize(0) )
bkg_model.plotOn(x_frame, RF.LineColor(ROOT.kRed) )

x_frame.Draw()
pad2.cd()

h_myy = inputHisto.ProjectionX()

plotsMyy = PlotSubtraction(varX,n_bins_yy,configFile["myySubtractionRange"],h_myy)
plotsMyy[0].Draw()
plotsMyy[1].Draw("EPSAME")

c1.SaveAs("Plots/"+str(args.catName)+"/"+args.catName+"_"+args.myy+"_"+args.mjj+"_background_myy.png")

pad1.cd()

# create plot of mjj now 
y_frame = mjj.frame(mjj_min, mjj_max, n_bins_jj)
data.plotOn(y_frame, RF.XErrorSize(0) )
bkg_model.plotOn(y_frame, RF.LineColor(ROOT.kRed) )
y_frame.Draw()

pad2.cd()

# plot the subtraction below
h_mjj = inputHisto.ProjectionY()
plotsMjj = PlotSubtraction(varY,n_bins_jj,configFile["mjjSubtractionRange"],h_mjj)
plotsMjj[0].Draw()
plotsMjj[1].Draw("EPSAME")
c1.SaveAs("Plots/"+str(args.catName)+"/"+args.catName+"_"+args.myy+"_"+args.mjj+"_background_mjj.png")

### make plots in 2D
### showing the 2D MC histogram, 2D PDF histogram and the residual between the two 

histo2D = bkg_model.createHistogram("bkg_histo",combWS.var(varX),RF.Binning(n_bins_yy), RF.YVar(combWS.var(varY) ,RF.Binning(n_bins_jj)))
histo2D.Scale( inputHisto.Integral() / histo2D.Integral() )

chi2 = 0.0
nbins = 0.0

h_residual = ROOT.TH2F("h_residual", "h_residual", n_bins_yy, myy_min, myy_max, n_bins_jj, mjj_min, mjj_max )
# reduce the binning for the residual plot
h_residual.RebinX(5)
histo2D.RebinX(5)
inputHisto.RebinX(5)

for i in range(1,histo2D.GetXaxis().GetNbins()+1):
    for j in range(1,histo2D.GetYaxis().GetNbins()+1):

        err_val = inputHisto.GetBinError(i,j)
        tmp_val = ( inputHisto.GetBinContent(i,j) - histo2D.GetBinContent(i,j) ) / err_val
        chi2 += tmp_val*tmp_val
        nbins += 1.0
        h_residual.SetBinContent(i,j,tmp_val)

chi2_final = chi2 / nbins

print('Reduced chi is :: ' , round(chi2_final,2))

## make 2D histogram plots
## clear the canvas
c1 = ROOT.TCanvas()

h_residual.Draw("COLZ")
c1.SetRightMargin(0.18)
c1.SetLeftMargin(0.15)

h_residual.GetYaxis().SetTitle("m_{bb} [GeV]")
h_residual.GetXaxis().SetTitle("m_{\gamma\gamma} [GeV]")
h_residual.GetZaxis().SetTitle("")
c1.SaveAs("Plots/"+str(args.catName)+"/"+args.catName+"_"+args.myy+"_"+args.mjj+"_residual.png")

inputHisto.Draw("COLZ")
c1.SetRightMargin(0.18)
c1.SetLeftMargin(0.15)

inputHisto.GetYaxis().SetTitle("m_{bb} [GeV]")
inputHisto.GetXaxis().SetTitle("m_{\gamma\gamma} [GeV]")
inputHisto.GetZaxis().SetTitle("")
c1.SaveAs("Plots/"+str(args.catName)+"/"+args.catName+"_"+args.myy+"_"+args.mjj+"_MC.png")

histo2D.Draw("COLZ")
c1.SetRightMargin(0.18)
c1.SetLeftMargin(0.15)

histo2D.GetYaxis().SetTitle("m_{bb} [GeV]")
histo2D.GetXaxis().SetTitle("m_{\gamma\gamma} [GeV]")
histo2D.GetZaxis().SetTitle("")
c1.SaveAs("Plots/"+str(args.catName)+"/"+args.catName+"_"+args.myy+"_"+args.mjj+"_PDF.png")

## create signal + background model

combWS.factory("SUM::model( nsig[0.0,-20.0,20.0]*signalModel, b[20.0,0.0,20000000.0]*bkg_only )")
total_model = combWS.pdf("model")

spur_signal_results = [ ]
peakPosition = 121.0

while peakPosition < 129.5:

    combWS.var(configFile["peakPositionName"]).setVal(peakPosition)

    total_model.fitTo(data,RF.Save(), RF.SumW2Error(True), RF.Offset(True), RF.Minimizer("Minuit2", "migrad"), RF.Extended())
    
    spur_signal = combWS.var("nsig").getVal()
    spur_signal_error = combWS.var("b").getError()
    z_spur = ( abs(spur_signal) / spur_signal_error )*100.0

    atuple = ( abs(spur_signal), z_spur, chi2_final)
    spur_signal_results.append(atuple)

    ## reset the number of signal and background

    combWS.var("nsig").setVal(0.0)
    combWS.var("b").setVal(0.0)

    peakPosition += 1.0

spur_signal_results.sort(key=lambda tup:tup[0],reverse=True)

## if no file exists, create new file 
if not os.path.exists(args.catName+".txt"):
    output_file = open(str(args.catName)+".txt", "w")
    output_file.write("spurious signal, Z_[spur] (%), chi2 / nDof"+'\n')
    

output_file = open(str(args.catName)+".txt", "a")
output_file.write(str(args.myy) + str(args.mjj) + ' & '  +''.join(str(round(item,3)) + " & " for item in spur_signal_results[0]) + '\n')
output_file.close()

del combWS
