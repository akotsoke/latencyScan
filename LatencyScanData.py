#!/bin/env python

# -*- coding: utf-8 -*-
"""

@author: Anastasia and Cameron Bravo (c.bravo@cern.ch)

"""
import ROOT as r
import numpy as np
from array import array
from optparse import OptionParser
from ROOT import TF1
#Define and parse options
parser = OptionParser()

parser.add_option("-i", "--infilename", type="string", dest="filename", default="LatencyData_Trimmed.root",
                  help="Specify Input Filename", metavar="filename")
parser.add_option("-o", "--outfilename", type="string", dest="outfilename", default="latencyAna.root",
                  help="Specify Output Filename", metavar="outfilename")

(options, args) = parser.parse_args()

r.gROOT.SetBatch()
f = r.TFile(options.filename)
t = f.latTree

Nhits_hs = {}
Nev_hs = {}

vfats= array( 'f' )
par= array( 'f' )
errY= array( 'f' )
errX=array('f')

par_NoPeaks= array( 'f' )
errY_NoPeaks= array( 'f' )
errX_NoPeaks=array('f')


par_Peak= array( 'f' )
errY_Peak= array( 'f' )
errX_Peak=array('f')

N_hits = [] 

def fline(x,par):
	if(x[0]>20 and x[0]<50):
		TF1.RejectPoint()
		return 0;

	return par[0];




for vfat in range(24):
    Nhits_hs[vfat] = r.TH1D("Nhits%i_h"%vfat,"Nhits%i_h"%vfat,256,-0.5,255.5)
    Nev_hs[vfat] = r.TH1D("Nev%i_h"%vfat,"Nev%i_h"%vfat,256,-0.5,255.5)
    N_hits.append([])
    pass

for evt in t:
    Nhits_hs[int(evt.vfatN)].AddBinContent(evt.lat + 1,evt.Nhits)
    Nev_hs[int(evt.vfatN)].AddBinContent(evt.lat + 1,evt.Nev)
    N_hits[int(evt.vfatN)].append(evt.Nhits)
    pass

N_hits= np.array(N_hits)
print np.amax(N_hits[1])   



canv = r.TCanvas("canv","canv",1000,1000)
canv.cd()
outF = r.TFile(options.outfilename,"RECREATE")

lat_ga = {}
fit_f = {}
fit_f_NoPeaks = {}
fit_f_Peak = {}

for vfat in range(0,24):
    fit_f[vfat] = r.TF1("fit%i_f"%vfat,"[0]",0,255)  #Fit function for the whole range
    fit_f_NoPeaks[vfat] = r.TF1("fit%i_f1"%vfat,fline,0,255,1) #Fit function for the range w/out the peaks
    fit_f_Peak[vfat]= r.TF1("fit%i_f_peak"%vfat,"[0]",30,40)  #Fit function for the peak range

    lat_ga[vfat] = r.TGraphAsymmErrors(Nhits_hs[vfat],Nev_hs[vfat])
    lat_ga[vfat].SetName("lat%i_ga"%vfat)
    vfats.append(vfat)
    lat_ga[vfat].Draw("AP")
    

    #Fit the whole range
    fit_f[vfat].SetParameter(0, 0.0005)
    lat_ga[vfat].Fit(fit_f[vfat],"Q")
    parameters = fit_f[vfat].GetParameters()
    par.append(parameters[0])
    errY.append(fit_f[vfat].GetParError(0))
    errX.append(0.0)
    lat_ga[vfat].Write()



    #Fit the range w/out the peaks
    fit_f_NoPeaks[vfat].SetParameter(0, 0.0005)
    lat_ga[vfat].Fit(fit_f_NoPeaks[vfat],"Q")
    parameters1 = fit_f_NoPeaks[vfat].GetParameters()
    par_NoPeaks.append(parameters1[0])
    errY_NoPeaks.append(fit_f_NoPeaks[vfat].GetParError(0))
    errX_NoPeaks.append(0.0)

   

    #Fit the Peaks
    fit_f_Peak[vfat].SetParameter(0, 0.03)
    lat_ga[vfat].Fit(fit_f_Peak[vfat],"R")
    parameters2 = fit_f_Peak[vfat].GetParameters()
    par_Peak.append(parameters2[0])
    errY_Peak.append(fit_f_NoPeaks[vfat].GetParError(0))
    errX_Peak.append(0.0)

    canv.SaveAs("VFAT%i_Nhits_vs_Latency.png" %vfat)


   


tgr= r.TGraphErrors(24,vfats,par,errX,errY)
tgr.SetName("Average_Nhits_per_Vfat")
tgr.SetMarkerStyle(20)
tgr.SetMarkerSize(0.7)
tgr.SetMarkerColor(r.kRed+1)
tgr.SetLineColor(r.kRed+1)
tgr.GetXaxis().SetLimits(0,24)
tgr.GetXaxis().SetTitle("VFAT Position")
tgr.GetYaxis().SetTitle("Average Nhits")
tgr.GetYaxis().SetTitleOffset(1.6)

tgr_NoPeaks= r.TGraphErrors(24,vfats,par_NoPeaks,errX_NoPeaks,errY_NoPeaks)
tgr_NoPeaks.SetMarkerStyle(20)
tgr_NoPeaks.SetMarkerSize(0.7)
tgr_NoPeaks.SetMarkerColor(r.kGreen)
tgr_NoPeaks.SetLineColor(r.kGreen)




tgr.Draw("AP")
tgr_NoPeaks.Draw("PSAME")
tgr.SetTitle("Average Nhits per Vfat")
canv.SaveAs("Average_Nhits_per_vfat.png")
tgr.Write()

par_fraction=array( 'f' )


par_Peak= np.array(par_Peak)
par_NoPeaks=np.array(par_NoPeaks)
par_fraction = np.array(par_Peak/par_NoPeaks)
vfats=np.array(vfats)
tgr_fraction= r.TGraph(24,vfats,par_fraction)
tgr_fraction.SetMarkerStyle(20)
tgr_fraction.SetMarkerSize(0.9)
tgr_fraction.SetMarkerColor(r.kBlue)
tgr_fraction.SetLineColor(r.kBlue)
tgr_fraction.Draw("AP")
tgr_fraction.GetXaxis().SetTitle("VFAT Position")
tgr_fraction.GetYaxis().SetTitle("Average_{InPeak}/Average_{OutPeak}")
tgr_fraction.GetYaxis().SetTitleOffset(1.3)
canv.SaveAs("ratio_InPeak_OutOfPeak.png")



outF.Close()

