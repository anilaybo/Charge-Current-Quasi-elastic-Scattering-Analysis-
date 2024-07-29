#!/usr/bin/env python3
import ROOT as r
import rootUtils as ut
import os, sys
#import numpy as np

from physlib import *
from ROOT import TChain, TSelector, TTree, TClass

input_file = "~/Analysis/CCQE.root"
geo_file = "~/Analysis/geofile_full.conical.Genie-TGeant4.root"

work_dir= '~/'


f = r.TFile(input_file)
t = f.cret
h = r.TClass('FairMCEventHeader')
TClass = r.FairMCEventHeader


def configure(input_file):
    f = r.TFile(input_file)
    t = f.cret
    t.MakeClass()
    r.gROOT.ProcessLine('.L cret.C')
    cret_ch = r.TChain('cret', 'cret')
    cret_ch.Add(input_file)
    cret = r.cret(cret_ch)
    return cret

def finish():
    os.system('rm cret.C cret.h')

nEnt = t.GetEntries()

g = r.TFile(geo_file)
sGeo = g.FAIRGeom       #intall the FairShip's geometry for running the fGeo
fGeo = r.gGeoManager    #the module that being able to use Geometric classes and functions



inGeo, tGS, tLS = 0., 0., 0.    #Selection cut

Hadron = [-130, -211, -321, -2212, 130, 211, 321, 2212]  
Lepton = [11, 12, 13, 14, 15, 16, 17, 18]
Chargeless = [-14, 22, 111, 130, 421, 2112]
print ('okay')

def GeometrySelection(Pos):
    gCheck = []
    fGeo.SetCurrentPoint(Pos[0], Pos[1], Pos[2])
    init = fGeo.FindNode()
    for i in [-1,1]:
        for j in [-1,1]:
            fGeo.SetCurrentDirection(i,j,0)         #tracking in the emulsion
            fGeo.FindNextBoundary()
            if fGeo.GetStep() > 0.1:
                gCheck.append(True)
    fGeo.SetCurrentPoint(Pos[0], Pos[1], Pos[2]+0.5)        # z dimension +0.5
    finl = fGeo.FindNode()
    if init.GetMotherVolume() == finl.GetMotherVolume():      # mother volume is brick
        gCheck.append(True)
    if sum(gCheck) != 5: return False
    return True

def LocationSelection(Mom, Pdg):
    LCheck = []
    SlopeX = Slope(Mom[0], Mom[2])
    SlopeY = Slope(Mom[1], Mom[2])
    if Mom[3] >1.:
        LCheck.append(True)
    if Pdg not in Chargeless:
        LCheck.append(True)
    if (SlopeX < 1. or SlopeY < 1.):
        LCheck.append(True)
    if sum(LCheck) != 3: return False
    return True


def Slope(Pi, Pz):                      #Calculation for x,y,z axis
    Slope = r.TMath.ATan(Pi/Pz)
    return Slope



def ScatteringAngle(S01,S02,S11,S12):
    SpaceAngle = (S01-S11)**2+(S02-S12)**2 	#Space angle square
    SpAngle = r.TMath.Sqrt(SpaceAngle)
    return SpAngle



def cut(Ee,Sa):     #input Ee and Sa for functions
    cCheck = [] 
    if Ee > 1.0: 
       cCheck.append(True)
    else:
       cCheck.append(False)
    if Sa**2 < (2*511.0*10**-6)/Ee:
       cCheck.append(True)
    else: 
       cCheck.append(False)
    if sum(cCheck) != 2: return False 
    return True
  


def MomentumCut(P):
    mCheck = []
    for i in range(len(P)):
        if P[i] > 1.0:
           mCheck.append(True)
        else:
           mCheck.append(False)
    if True not in mCheck: return False
    return True


def counter(x):
    SelectEvent = []
    for Event in xrange(x):
        while x < Event:
            x += 1
        SelectEvent.append(x)
        CounterEvent = (SelectEvent/x)
    return x

def Multiplicity(PdgArr, RArr):
    Mult = 0
    for i in range(len(PdgArr)):
        if PdgArr[i] not in RArr:
            Mult += 1
    return Mult


def Bjorken(mom4_nu, mom4_nucl, mom4_lept):
    BjorkX = 2*(mom4_nu*mom4_lept)/(2*(mom4_nucl*(mom4_nu-mom4_lept)))
    return BjorkX


def Inelasticity(mom4_nu, mom4_nucl, mom4_lept):
    InelY = (mom4_nucl*(mom4_nu-mom4_lept))/(mom4_nucl*mom4_nu)
    return InelY


print ('its okay')

'''
#Electron Scattering Angle and Electron Energy Plot
c = r.TCanvas()
h = r.TH2D("htemp", 'Correlation between the electron scattering angle and electron energy for neutrino CCQE process', 40,0,50, 40,0,1)
'''


h = {}
ut.bookHist(h, 'EE', 'CCQE-Electron Energy', 20,0,250)
ut.bookHist(h, 'ESA', 'CCQE-Electron Scattering Angle', 40,0,4)
ut.bookHist(h, 'MCP', 'CCQE-Multiplicity- Number of Charge Particle at Neutrino Vertex', 20,0,20)
ut.bookHist(h, 'E1S', 'Correlation between the electron scattering and electron energy for neutrino CCQE process', 50,0,50, 60,0,1)
ut.bookHist(h, 'BjorkX', 'Bjorken X Distribution', 100,0,1)
ut.bookHist(h, 'BjorkXS', 'Bjorken X Distribution (Selected)', 100,0,1)
ut.bookHist(h, 'InelY', 'Inelasticity Y Distribution', 100,0,1)
ut.bookHist(h, 'InelYS', 'Inelasticity Y Distribution (Selected)', 100,0,1)




def makePlot(TClass):
    #Electron Energy Histogram
    ut.bookCanvas(h, key='Electron Energy', title='Electron Energy', nx=1920, ny=1080, cx=1, cy=1)
    c = h['Electron Energy'].cd(1)
    r.gStyle.SetOptStat('men')
    h['EE'].Draw()
    h['EE'].SetFillStyle(3345)
    h['EE'].SetFillColor(2)
    h['EE'].SetXTitle('Electron Energy')
    h['EE'].SetYTitle('Number of Event')
    h['Electron Energy'].Print(work_dir+'/Analysis/histogram/CCQE-EEnergy.png')
    #Electron Scattering Angle
    ut.bookCanvas(h, key='Electron Scattering Angle', title='Electron Scattering Angle', nx=1920, ny=1080, cx=1, cy=1)
    #r.gPad.SetLogy(1)
    r.gStyle.SetOptStat('en')
    h['ESA'].Draw()
    h['ESA'].SetFillStyle(3345)
    h['ESA'].SetFillColor(4)
    h['ESA'].SetXTitle('Electron Scattering Angle(mrad)')
    h['ESA'].SetYTitle('Number of Event')
    h['Electron Scattering Angle'].Print(work_dir+'/Analysis/histogram/CCQE-EScatAng.png')
    #Multiplicity (Number charged of particle at neutrino vertex)
    ut.bookCanvas(h, key='Multiplicity', title='Multiplicity (Number of charged particle at neutrino vertex)', nx=1920, ny=1080, cx=1, cy=1)
    c = h['Multiplicity'].cd(1)
    r.gStyle.SetOptStat('en')
    h['MCP'].Draw()
    h['MCP'].SetFillStyle(3345)
    h['MCP'].SetFillColor(6)
    h['MCP'].SetXTitle('Multiplicitiy')
    h['MCP'].SetYTitle('Number of Event')
    h['Multiplicity'].Print(work_dir+'/Analysis/histogram/CCQE-EMultiplicity.png')
    #Electron scattering angle and electron energy plot
    ut.bookCanvas(h, key='Electron Scat and Energy', title='Correlation between the electron scattering and electron energy for neutrino CCQE process', nx=1920, ny=1080, cx=1, cy=1)
    c = h['Electron Scat and Energy'].cd(1)
    r.gStyle.SetOptStat('menr')
    h['E1S'].Draw('Box')
    h['E1S'].SetXTitle('Energy')
    h['E1S'].SetYTitle('Scattering Angle')
    h['Electron Scat and Energy'].Print(work_dir+'/Analysis/histogram/CCQE-EScatAngandEng.png')
    #Bjorken X Distribution
    ut.bookCanvas(h, key='Bjork', title='Bjorken X Distribution', nx=1920, ny=1080, cx=1, cy=1)
    c = h['BjorkX'].Draw()
    h['BjorkXS'].SetFillStyle(3345)
    h['BjorkXS'].SetFillColor(2)
    h['BjorkXS'].Draw('same')
    h['Bjork'].Print(work_dir+'/Analysis/histogram/CCQE-BjorkenX.png')
    #Inelasticity Y Distribution
    ut.bookCanvas(h, key='Inel', title='Inelasticity Y Distribution', nx=1920, ny=1080, cx=1, cy=1)
    c= h['InelY'].Draw()
    h['InelYS'].SetFillStyle(3345)
    h['InelYS'].SetFillColor(2)
    h['InelYS'].Draw('same')
    h['Inel'].Print(work_dir+'/Analysis/histogram/CCQE-Inelasticity.png')





selected = 0.
allx = 0.

for Event in range(nEnt):
    t.GetEntry(Event)
    if (t.IntInGeo.at(0)):      #True=0, False=1
        EDauPdg = []
        EDauPos_i, EDauPos_j, EDauPos_k = [], [], []
        Mom_i, Mom_j, Mom_k, Mom_l = [], [], [], []
        PVPdg = []       #Primary Vertex Pdg
        GS, LS = [], []
        allx +=1.0
        P = []
        Pdg = []
        mom4_nu, mom4_nucl, mom4_lept = r.TLorentzVector(0., 0., 0., 0.,), r.TLorentzVector(0., 0., 0., 0.,), r.TLorentzVector(0., 0., 0., 0.,)
        nuEnergy= 0.
        for Vertex in range(t.VertexInfo.size()):
                if t.VertexInfo.at(Vertex) ==0:
                       # print ('nu', t.Energy.at(Vertex))
                        S01 = Slope(t.Px.at(Vertex), t.Pz.at(Vertex))
                        S02 = Slope(t.Py.at(Vertex), t.Pz.at(Vertex))
                        nuEnergy = t.Energy.at(Vertex)
                        mom4_nu += r.TLorentzVector(t.Px.at(Vertex), t.Py.at(Vertex), t.Pz.at(Vertex), nuEnergy)
                        mom4_nucl += r.TLorentzVector(0., 0., 0., (0.9383+0.9396)/2)
                if t.VertexInfo.at(Vertex) ==1:
                        P.append(t.P.at(Vertex))
                        Pos = []
                        Pos.append(t.StartX.at(Vertex))
                        Pos.append(t.StartY.at(Vertex))
                        Pos.append(t.StartZ.at(Vertex))
                        Mom_i.append(t.Px.at(Vertex))
                        Mom_j.append(t.Py.at(Vertex))
                        Mom_k.append(t.Pz.at(Vertex))
                        Mom_l.append(t.P.at(Vertex))
                        PVPdg.append(t.PdgCode.at(Vertex))
                        Pdg.append(t.PdgCode.at(Vertex))
                        if PVPdg[-1] in Lepton:
                            mom4_lept += r.TLorentzVector(t.Px.at(Vertex), t.Py.at(Vertex), t.Pz.at(Vertex), t.Energy.at(Vertex))
                        if t.PdgCode.at(Vertex) == 11:
                            #print ('pe', t.Energy.at(Vertex))
                            S11 = Slope(t.Px.at(Vertex), t.Pz.at(Vertex))
                            S12 = Slope(t.Py.at(Vertex), t.Pz.at(Vertex))
                            E1 = t.Energy.at(Vertex)
                            EPos = Pos
                            EDauPos_i.append(t.StartX.at(Vertex))
                            EDauPos_j.append(t.StartY.at(Vertex))
                            EDauPos_k.append(t.StartZ.at(Vertex))
                            EDauPdg.append(t.PdgCode.at(Vertex))
                            E1 = t.Energy.at(Vertex)
                            h['EE'].Fill(E1)

        if 11 in PVPdg: 
            S= ScatteringAngle(S01,S02,S11,S12)
            M = Multiplicity(Pdg,Chargeless)
            EDauPos = [EDauPos_i[0], EDauPos_j[0], EDauPos_k[0]]


            Bx = Bjorken(mom4_nu,mom4_nucl,mom4_lept)
            Iy = Inelasticity(mom4_nu,mom4_nucl,mom4_lept)
 
            D = cut(E1,S)

            #print ('Scattering Angle', S)
            #print ('Energy', E1)
            #print ('Direction cut', D)


            h['ESA'].Fill(S)
            h['MCP'].Fill(M)

            h['BjorkX'].Fill(Bx)
            h['BjorkXS'].Fill(Bx)
            h['InelY'].Fill(Iy)
            h['InelYS'].Fill(Iy)


            if M ==1:
               h['E1S'].Fill(E1,S)

            #print 'Multiplicity', M

            if MomentumCut(P) and cut(E1,S):
               selected += 1.0

            #h.Fill(E1,S)
      
            #Geometry Selection Check // Checked at Each Vertex
            if GeometrySelection(Pos):
               GS.append(True)
            else: GS.append(False)
            if GeometrySelection(EDauPos):
               GS.append(True)
            else: GS.append(False)


            #Location Selection Check // Checked at Neutrino Vertex
            for b in range(len(Mom_i)):
                Mom = [Mom_i[b],Mom_j[b],Mom_k[b],Mom_l[b]]
            if LocationSelection(Mom, PVPdg[b]):
                LS.append(True)
            else: LS.append(False)

      
#h.Draw('Scat')
#h.SetXTitle('Energy')
#h.SetYTitle('Scattering Angle')
#c.Print('CCQE.pdf')


print ('selected ratio', selected/allx)
print ('selected', selected)
print ('allx', allx)
print ('nEnt', nEnt)



makePlot(TClass)
