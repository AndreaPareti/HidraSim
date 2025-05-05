import pandas as pd
import numpy as np
import ROOT 
import uproot


datafolder = "/home/apareti/hidra/TB24/ElectronEnergyScan/"
#datafolder = "/home/apareti/hidra/TB24/PionEnergyScan/"
mcfolder = "/home/apareti/hidra/HidraSim/analysis/"
treename = "Ftree"

# do not show plots when running through ssh
ROOT.gROOT.SetBatch(True)





def GetDF(run, Cut, filename, energy): 
    print("Using file: ", filename, energy)
    root_file = uproot.open(filename)
    tree = root_file[treename]
    data = tree.arrays(cut=Cut, library="pd")

    # define Asymmetry variable 
    data["AsymS"] = (data["TS11"] - data["TS15"]) / (data["TS11"] + data["TS15"] )
    data["AsymC"] = (data["TC11"] - data["TC15"]) / (data["TC11"] + data["TC15"] )

    # define partial barycenter variable 
    data["BaryS"] = (data["TS00"]-28.3*data["TS11"]+28.3*data["TS15"])/(data["TS00"]+data["TS11"]+data["TS15"])
    data["BaryC"] = (data["TC00"]-28.3*data["TC11"]+28.3*data["TC15"])/(data["TC00"]+data["TC11"]+data["TC15"])

    # input truth energy to dataframe
    data["TruthE"] = energy

    return data


def DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, markerStyle1=20, markerColor1=ROOT.kBlue, markerStyle2=20, markerColor2=ROOT.kRed+1, labelY="Counts", log=0):
    cCanva = ROOT.TCanvas(ctitle, "{0};{1};{2}".format(labelPlot, labelX1, labelY), 1400, 1200)
    cCanva.SetLogy(log)
    myHist1 = ROOT.TH1D(ctitle+"1", "{0};{1};{2}".format(labelPlot+"1", labelX1, labelY), nbinX, xmin, xmax)
    myHist2 = ROOT.TH1D(ctitle+"2", "{0};{1};{2}".format(labelPlot+"2", labelX1, labelY), nbinX, xmin, xmax)
    for xval1 in varX1:
        myHist1.Fill(xval1)
    for xval2 in varX2:
        myHist2.Fill(xval2)
        
    int1 = myHist1.Integral()
    scale1 = 1/int1; myHist1.Scale(scale1)
    int2 = myHist2.Integral()
    scale2 = 1/int2; myHist2.Scale(scale2)
    

    myMax = max(myHist1.GetMaximum(), myHist2.GetMaximum())
    scalemax = 1.1
    if(log>0):
        scalemax = 30
    myHist1.SetMaximum(myMax*scalemax); myHist2.SetMaximum(myMax*scalemax)

    myHist1.SetMarkerStyle(markerStyle1); myHist1.SetMarkerColor(markerColor1); myHist1.SetLineWidth(2); myHist1.SetLineColor(markerColor1)
    myHist2.SetMarkerStyle(markerStyle2); myHist2.SetMarkerColor(markerColor2); myHist2.SetLineWidth(2); myHist2.SetLineColor(markerColor2)

    myHist1.SetTitle(labelPlot)
    myHist1.Draw("")
    myHist2.Draw("HIST same")
    leg = ROOT.TLegend(0.70, 0.70, 0.87, 0.87)
    leg.AddEntry(myHist1, leg1, "l")
    leg.AddEntry(myHist2, leg2, "l")
    leg.SetLineWidth(0)
    leg.Draw()
    cCanva.SetLeftMargin(0.14)
    cCanva.SetRightMargin(0.09)

    cCanva.SaveAs(outfile)




def main():
    print("Hello there")
    dataruns = ["0786", "0766", "0772", "0774", "0775", "0778", "0779", "0792"]
    mcruns = ["0", "1", "2", "3", "4", "5", "6", "7"]
    energies = [10, 20, 30, 40, 60, 80, 100, 120]

    #dataruns = ["0766"]
    #mcruns = ["1"]
    #energies = [20]


    myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (MCounter<200) & (TailC<300) & (C2>160) & ( (totLeakage - L20)<5000 )"

    varProf = "YDWC2"
    cut_x_min = [-19.83, -16.74, -16.22, -15.95, -15.60, -16.12, -16.07, -15.50]
    cut_x_max = [23.90, 22.19, 23.27, 23.44, 24.27, 23.79, 23.63, 24.12]
    #cut_y_min = [-26.54, -25.80, -26.15, -26.15, -26.39, -25.63, -25.63, -26.03]
    #cut_y_max = [13.38, 10.89, 9.72, 9.50, 9.86, 10.89, 10.54, 10.17]


    cut_y_min = [-20, -20, -20, -20, -20, -20, -20, -20]
    cut_y_max = [4, 4, 4, 4, 4, 4, 4, 4]



    #dataruns = ["1000", "0967", "0966", "0965", "0963", "0962"]
    #ncruns = ["1", "3", "4", "5", "6", "7"]
    #energies = [20, 40, 60, 80, 100, 120]
    #myCut = "(abs(XDWC2 - XDWC1) < 5) & (abs(YDWC2 - YDWC1)<5) & (totPMTSene>0) & (PShower<450) & (TailC<400) & (TailC>170) & (totLeakage<6500) & (MCounter<160) & (PShower>350) & (YDWC2>-20) & (YDWC2<5) & (XDWC2>-20) & (XDWC2<20)"




    # Store ntuples as pandas dataframes
    #df10, df20, df30, df40, df60, df80, df100, df120
    dfs = []

    # array of dataframe with corrected energies
    dfCorrected_array = []

    # Arrays to store reco energy parameters
    MeanVec_S=[]; RmsVec_S=[]; MeanErrVec_S=[]; RmsErrVec_S=[] 
    MeanVec_C=[]; RmsVec_C=[]; MeanErrVec_C=[]; RmsErrVec_C=[] 
    MeanVec_Comb=[]; RmsVec_Comb=[]; MeanErrVec_Comb=[]; RmsErrVec_Comb=[] 
    ROOT.gStyle.SetOptStat(0)
    
    for index, (datarun, mcrun, energy) in enumerate(zip(dataruns, mcruns, energies)):
        datafilename = datafolder+"physics_sps2024_run" + datarun + ".root"
        mcfilename = mcfolder+"physics_HidraSimout_Run" + mcrun + ".root"

        ROOT.gStyle.SetOptFit(0)
        print("Current cuts on DWCs: ", cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index])
        dwcCut = " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 
        #CurrentCut = myCut + " & (XDWC2 > {0}) & (XDWC2 < {1}) & (YDWC2 > {2}) & (YDWC2 < {3})".format(cut_x_min[index], cut_x_max[index], cut_y_min[index], cut_y_max[index]) 
        CurrentCut = myCut + dwcCut

        dfdata = GetDF(datarun, CurrentCut, datafilename, energy)
        dfmc = GetDF(mcrun, "(totPMTSene>0) "+dwcCut, mcfilename, energy)

        var = "totPMTSene"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 100; xmin = energy-0.5*energy; xmax = energy+0.5*energy
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts")

        var = "totPMTCene"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 100; xmin = energy-0.5*energy; xmax = energy+0.5*energy
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts")

        var = "XDWC2"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 100; xmin = -30; xmax = 30
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts")

        var = "YDWC2"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 100; xmin = -30; xmax = 30
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts")

        #var = "PShower"
        #dfdata["PShower"] = dfdata["PShower"]-300
        #dfdata = dfdata[dfdata["PShower"]>0]

        var = "TS00"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 60; xmin = 0.; xmax = energy*1.2
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts", log=1)

        var = "TS11"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 60; xmin = 0.; xmax = energy
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts", log=1)

        var = "TS15"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 60; xmin = 0.; xmax = energy
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts", log=1)

        var = "TC00"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 60; xmin = 0.; xmax = energy*1.2
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts", log=1)

        var = "TC11"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 60; xmin = 0.; xmax = energy
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts", log=1)

        var = "TC15"
        outfile = "comparison_"+var+"_{0}GeV.png".format(energy); labelPlot= var+" {0} GeV".format(energy)
        varX1 = dfdata[var]; varX2 = dfmc[var]; ctitle = var; labelX1 = var+" [GeV]"; labelY = "Normalised counts"; labelX2 = labelX1
        leg1 = "TB data"; leg2 = "Simulation"; nbinX = 60; xmin = 0.; xmax = energy
        DrawTwoHist(outfile, ctitle, varX1, varX2, labelPlot, labelX1, labelX2, nbinX, xmin, xmax, leg1, leg2, labelY="Normalised counts", log=1)




        ##################
        ## Profile Plots of Sci energy vs Y coordinate
        #################

        var = "YDWC2"
        #ts00 = "TS00"; ts11 = "TS11"; ts15="TS15"; tc00 = "TC00"; tc11 = "TC11"; ts15="TC15"
        outfile = "DataEoverYProfileSci_{0}GeV.png".format(energy); labelPlot=var+" {0} GeV".format(energy)
        dataProfT00 = ROOT.TProfile("TS00dataprof_{0}GeV".format(energy), "TS00 data profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        dataProfT11 = ROOT.TProfile("TS11dataprof_{0}GeV".format(energy), "TS11 data profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        dataProfT15 = ROOT.TProfile("TS15dataprof_{0}GeV".format(energy), "TS15 data profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        dataProftotene = ROOT.TProfile("totSenedataprof_{0}GeV".format(energy), "TotSene data profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)

        for ydwc, ts00, ts11, ts15, totS in zip(dfdata["YDWC2"], dfdata["TS00"], dfdata["TS11"], dfdata["TS15"], dfdata["totPMTSene"]):
            dataProfT00.Fill(ydwc, ts00)
            dataProfT11.Fill(ydwc, ts11)
            dataProfT15.Fill(ydwc, ts15)
            dataProftotene.Fill(ydwc, totS)

        mcProfT00 = ROOT.TProfile("TS00mcprof_{0}GeV".format(energy), "TS00 mc profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        mcProfT11 = ROOT.TProfile("TS11mcprof_{0}GeV".format(energy), "TS11 mc profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        mcProfT15 = ROOT.TProfile("TS15mcprof_{0}GeV".format(energy), "TS15 mc profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        mcProftotene = ROOT.TProfile("totSenemcprof_{0}GeV".format(energy), "TotSene mc profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)

        for ydwc, ts00, ts11, ts15, totS in zip(dfmc["YDWC2"], dfmc["TS00"], dfmc["TS11"], dfmc["TS15"], dfmc["totPMTSene"]):
            mcProfT00.Fill(ydwc, ts00)
            mcProfT11.Fill(ydwc, ts11)
            mcProfT15.Fill(ydwc, ts15)
            mcProftotene.Fill(ydwc, totS)

        cProfS = ROOT.TCanvas("cProfS_{0}GeV".format(energy),"cProfS_{0}GeV".format(energy), 1400, 1200)

        maximum = max(dataProftotene.GetMaximum(), mcProftotene.GetMaximum())
        dataProftotene.SetMaximum(maximum*1.15); mcProftotene.SetMaximum(maximum*1.15)
        dataProftotene.SetTitle("Scintillation channel energy profile over Y coordinate, {0} GeV".format(energy))
        dataProftotene.GetXaxis().SetTitle("YDWC2 [mm]")
        dataProftotene.GetYaxis().SetTitle("Energy [GeV]")



        dataProftotene.SetMarkerSize(1); dataProftotene.SetMarkerStyle(20); dataProftotene.SetMarkerColor(ROOT.kAzure+1);  dataProftotene.SetLineWidth(2); dataProftotene.SetLineColor(ROOT.kAzure+1); dataProftotene.Draw()    
        dataProfT00.SetMarkerSize(1); dataProfT00.SetMarkerStyle(20); dataProfT00.SetMarkerColor(ROOT.kBlack);  dataProfT00.SetLineWidth(2); dataProfT00.SetLineColor(ROOT.kBlack); dataProfT00.Draw("same")
        dataProfT11.SetMarkerSize(1); dataProfT11.SetMarkerStyle(20); dataProfT11.SetMarkerColor(ROOT.kRed); dataProfT11.SetLineWidth(2); dataProfT11.SetLineColor(ROOT.kRed); dataProfT11.Draw("same")
        dataProfT15.SetMarkerSize(1); dataProfT15.SetMarkerStyle(20); dataProfT15.SetMarkerColor(ROOT.kGreen+2);  dataProfT15.SetLineWidth(2); dataProfT15.SetLineColor(ROOT.kGreen+2); dataProfT15.Draw("same")

        mcProftotene.SetMarkerSize(2); mcProftotene.SetMarkerStyle(51); mcProftotene.SetMarkerColor(ROOT.kAzure+1);  mcProftotene.SetLineWidth(2); mcProftotene.SetLineColor(ROOT.kAzure+1); mcProftotene.Draw("same")    
        mcProfT00.SetMarkerSize(2); mcProfT00.SetMarkerStyle(51); mcProfT00.SetMarkerColor(ROOT.kBlack);  mcProfT00.SetLineWidth(2); mcProfT00.SetLineColor(ROOT.kBlack); mcProfT00.Draw("same")
        mcProfT11.SetMarkerSize(2); mcProfT11.SetMarkerStyle(51); mcProfT11.SetMarkerColor(ROOT.kRed); mcProfT11.SetLineWidth(2); mcProfT11.SetLineColor(ROOT.kRed); mcProfT11.Draw("same")
        mcProfT15.SetMarkerSize(2); mcProfT15.SetMarkerStyle(51); mcProfT15.SetMarkerColor(ROOT.kGreen+2);  mcProfT15.SetLineWidth(2); mcProfT15.SetLineColor(ROOT.kGreen+2); mcProfT15.Draw("same")


        leg = ROOT.TLegend(0.77, 0.77, 0.95, 0.95)
        leg.AddEntry(dataProftotene, "totPMTSene (Data)", "pl")
        leg.AddEntry(dataProfT00, "TS00 (Data)", "pl")
        leg.AddEntry(dataProfT15, "TS15 (Data)", "pl")
        leg.AddEntry(dataProfT11, "TS11 (Data)", "pl")
        leg.AddEntry(mcProftotene, "totPMTSene (Sim)", "pl")
        leg.AddEntry(mcProfT00, "TS00 (Sim)", "pl")
        leg.AddEntry(mcProfT15, "TS15 (Sim)", "pl")
        leg.AddEntry(mcProfT11, "TS11 (Sim)", "pl")

        leg.Draw()
        cProfS.SaveAs("TowerSciProfileY{0}GeV.png".format(energy))




        ##################
        ## Profile Plotc of Cer energy vs Y coordinate
        #################

        var = "YDWC2"
        #tc00 = "TC00"; tc11 = "TC11"; tc15="TC15"; tc00 = "TC00"; tc11 = "TC11"; tc15="TC15"
        outfile = "DataEoverYProfileCer_{0}GeV.png".format(energy); labelPlot=var+" {0} GeV".format(energy)
        dataProfT00 = ROOT.TProfile("TC00dataprof_{0}GeV".format(energy), "TC00 data profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        dataProfT11 = ROOT.TProfile("TC11dataprof_{0}GeV".format(energy), "TC11 data profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        dataProfT15 = ROOT.TProfile("TC15dataprof_{0}GeV".format(energy), "TC15 data profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        dataProftotene = ROOT.TProfile("totCenedataprof_{0}GeV".format(energy), "TotCene data profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)

        for ydwc, tc00, tc11, tc15, totS in zip(dfdata["YDWC2"], dfdata["TC00"], dfdata["TC11"], dfdata["TC15"], dfdata["totPMTCene"]):
            dataProfT00.Fill(ydwc, tc00)
            dataProfT11.Fill(ydwc, tc11)
            dataProfT15.Fill(ydwc, tc15)
            dataProftotene.Fill(ydwc, totS)

        mcProfT00 = ROOT.TProfile("TC00mcprof_{0}GeV".format(energy), "TC00 mc profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        mcProfT11 = ROOT.TProfile("TC11mcprof_{0}GeV".format(energy), "TC11 mc profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        mcProfT15 = ROOT.TProfile("TC15mcprof_{0}GeV".format(energy), "TC15 mc profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)
        mcProftotene = ROOT.TProfile("totCenemcprof_{0}GeV".format(energy), "TotCene mc profile over YDWC2, {0} GeV".format(energy), 100, -22, 7)

        for ydwc, tc00, tc11, tc15, totC in zip(dfmc["YDWC2"], dfmc["TC00"], dfmc["TC11"], dfmc["TC15"], dfmc["totPMTCene"]):
            mcProfT00.Fill(ydwc, tc00)
            mcProfT11.Fill(ydwc, tc11)
            mcProfT15.Fill(ydwc, tc15)
            mcProftotene.Fill(ydwc, totC)



        cProfC = ROOT.TCanvas("cProfC_{0}GeV".format(energy),"cProfC_{0}GeV".format(energy), 1400, 1200)

        maximum = max(dataProftotene.GetMaximum(), mcProftotene.GetMaximum())
        dataProftotene.SetMaximum(maximum*1.15); mcProftotene.SetMaximum(maximum*1.15)
        dataProftotene.SetTitle("Cerenkov channel energy profile over Y coordinate, {0} GeV".format(energy))
        dataProftotene.GetXaxis().SetTitle("YDWC2 [mm]")
        dataProftotene.GetYaxis().SetTitle("Energy [GeV]")



        dataProftotene.SetMarkerSize(1); dataProftotene.SetMarkerStyle(20); dataProftotene.SetMarkerColor(ROOT.kAzure+1);  dataProftotene.SetLineWidth(2); dataProftotene.SetLineColor(ROOT.kAzure+1); dataProftotene.Draw()    
        dataProfT00.SetMarkerSize(1); dataProfT00.SetMarkerStyle(20); dataProfT00.SetMarkerColor(ROOT.kBlack);  dataProfT00.SetLineWidth(2); dataProfT00.SetLineColor(ROOT.kBlack); dataProfT00.Draw("same")
        dataProfT11.SetMarkerSize(1); dataProfT11.SetMarkerStyle(20); dataProfT11.SetMarkerColor(ROOT.kRed); dataProfT11.SetLineWidth(2); dataProfT11.SetLineColor(ROOT.kRed); dataProfT11.Draw("same")
        dataProfT15.SetMarkerSize(1); dataProfT15.SetMarkerStyle(20); dataProfT15.SetMarkerColor(ROOT.kGreen+2);  dataProfT15.SetLineWidth(2); dataProfT15.SetLineColor(ROOT.kGreen+2); dataProfT15.Draw("same")

        mcProftotene.SetMarkerSize(2); mcProftotene.SetMarkerStyle(51); mcProftotene.SetMarkerColor(ROOT.kAzure+1);  mcProftotene.SetLineWidth(2); mcProftotene.SetLineColor(ROOT.kAzure+1); mcProftotene.Draw("same")    
        mcProfT00.SetMarkerSize(2); mcProfT00.SetMarkerStyle(51); mcProfT00.SetMarkerColor(ROOT.kBlack);  mcProfT00.SetLineWidth(2); mcProfT00.SetLineColor(ROOT.kBlack); mcProfT00.Draw("same")
        mcProfT11.SetMarkerSize(2); mcProfT11.SetMarkerStyle(51); mcProfT11.SetMarkerColor(ROOT.kRed); mcProfT11.SetLineWidth(2); mcProfT11.SetLineColor(ROOT.kRed); mcProfT11.Draw("same")
        mcProfT15.SetMarkerSize(2); mcProfT15.SetMarkerStyle(51); mcProfT15.SetMarkerColor(ROOT.kGreen+2);  mcProfT15.SetLineWidth(2); mcProfT15.SetLineColor(ROOT.kGreen+2); mcProfT15.Draw("same")


        leg = ROOT.TLegend(0.77, 0.77, 0.95, 0.95)
        leg.AddEntry(dataProftotene, "totPMTCene (Data)", "pl")
        leg.AddEntry(dataProfT00, "TC00 (Data)", "pl")
        leg.AddEntry(dataProfT15, "TC15 (Data)", "pl")
        leg.AddEntry(dataProfT11, "TC11 (Data)", "pl")
        leg.AddEntry(mcProftotene, "totPMTCene (Sim)", "pl")
        leg.AddEntry(mcProfT00, "TC00 (Sim)", "pl")
        leg.AddEntry(mcProfT15, "TC15 (Sim)", "pl")
        leg.AddEntry(mcProfT11, "TC11 (Sim)", "pl")

        leg.Draw()
        cProfC.SaveAs("TowerCerProfileY{0}GeV.png".format(energy))







if __name__ == "__main__":
    main()
