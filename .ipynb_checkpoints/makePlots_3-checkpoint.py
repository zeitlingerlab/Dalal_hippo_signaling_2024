#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


def readSasaVals(dirname, file1, file2, fileComb, numPoints):
    fp1 = open(dirname + file1, "r")
    fp2 = open(dirname + file2, "r")
    fpC = open(dirname + fileComb, "r")
    ret = np.zeros((numPoints,))
    for i in range(numPoints):
        l1 = float(fp1.readline())
        l2 = float(fp2.readline())
        lC = float(fpC.readline())
        combSasa = l1 + l2
        Δsasa = (combSasa - lC) / 2
        ret[i] = Δsasa
    return ret



sasaP1W = readSasaVals("/n/projects/cm2363/tead4Simulation/amber_endsRestrained/production/sasas/", "proteinA.dat", "nucleic.dat", "proteinAnucleic.dat", 50000)
sasaP2W = readSasaVals("/n/projects/cm2363/tead4Simulation/amber_endsRestrained/production/sasas/", "proteinB.dat", "nucleic.dat", "proteinBnucleic.dat", 50000)
sasaPPW = readSasaVals("/n/projects/cm2363/tead4Simulation/amber_endsRestrained/production/sasas/","proteinA.dat", "proteinB.dat", "proteinprotein.dat", 50000)
sasaP1D = readSasaVals("/n/projects/cm2363/tead4Simulation/amber_mutation3/production/sasas/", "lefthalfprotein.dat", "nucleic.dat", "lefthalf.dat", 10000)
sasaP2D = readSasaVals("/n/projects/cm2363/tead4Simulation/amber_mutation3/production/sasas/", "righthalfprotein.dat", "nucleic.dat", "righthalf.dat", 10000)
sasaPPD = readSasaVals("/n/projects/cm2363/tead4Simulation/amber_mutation3/production/sasas/","lefthalf.dat", "righthalf.dat", "proteinprotein.dat", 10000)


def makeHists(data1, data2, color1, color2, ymin, ymax, name):


    fig = plt.figure(figsize=(3,2), dpi=600)
    axHist = fig.add_axes([0.3, 0.3, 0.6, 0.6])
    axHist.set_xlim((ymin, ymax))
    axHist.set_ylim((0, 2e-2))
    axHist.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
    def addPlot(dats, color, axHist, fill):
        
        histOfVals = np.histogram(dats[dats>1], bins=200)#, range=(ymin, ymax))
        histX = histOfVals[1][:-1]
        histogram = histOfVals[0]
        histogram = histogram / np.sum(histogram)
        if(fill):
            axHist.fill_between(histX, histogram, color = color, alpha=0.5)
        else:
            axHist.plot(histX, histogram, color = color, linestyle='solid')
        medValue = np.median(dats)
        axHist.plot( [medValue, medValue], [np.max(histogram) / 8,0],color=color, linestyle='solid' )

    addPlot(data1, color1, axHist, True)
    addPlot(data2, color2, axHist, False)
    #addPlot(sasaP1, "protein A-DNA", '#069AA8', axHist)
    #addPlot(sasaP2, "protein B-DNA", '#014D6D', axHist)
    fig.savefig(name + ".png")
    fig.savefig(name + ".pdf")

makeHists(sasaPPW, sasaPPD, '#CC6677', '#006677', 0, 400, "proteinProtein")
makeHists(sasaP1W, sasaP1D, '#CC6677', '#006677', 400, 1400, "protein1")
makeHists(sasaP2W, sasaP2D, '#CC6677', '#006677', 400, 1400, "protein2")

#addPlot(sasaPP, "protein-protein", '#CC6677',axHist)
#addPlot(sasaLP, "protein-left site", '#46DAC8',axHist)
#addPlot(sasaRP, "protein-right site", '#266DBD',axHist)
#axHist.set_xlabel("Buried surface area (Å²)")
#axHist.set_ylabel("Frequency")
#plt.savefig(outputName + ".png")
#plt.savefig(outputName + ".pdf")
