#! /usr/bin/env python

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--infile', type='string', action='store',
    dest='files',
    help='Input files')

parser.add_option('--outfile', type='string', action='store',
    default='outplots.root',
    dest='outname',
    help='Name of output file')

parser.add_option('--verbose', action='store_true',
    default=False,
    dest='verbose',
    help='Print debugging info')

parser.add_option('--maxEvts', type='int', action='store',
    default=-1,
    dest='maxEvts',
    help='Number of events to run. -1 is all events')

(options, args) = parser.parse_args()
argv = []

import ROOT, sys, copy
from DataFormats.FWLite import Events, Handle
from operator import itemgetter
import collections

infiles = [options.files]
printGen = True
events = Events (infiles)

jetAK4PtHandle, jetAK4Label, jetAK4PtProduct = Handle ("std::vector<float>"), ("jetsAK4"), ("jetAK4Pt")
jetAK4EtaHandle, jetAK4Label, jetAK4EtaProduct = Handle ("std::vector<float>"), ("jetsAK4"), ("jetAK4Eta")
jetAK4PhiHandle, jetAK4Label, jetAK4PhiProduct = Handle ("std::vector<float>"), ("jetsAK4"), ("jetAK4Phi")
jetAK4EHandle, jetAK4Label, jetAK4EProduct = Handle ("std::vector<float>"), ("jetsAK4"), ("jetAK4E")
jetAK4MassHandle, jetAK4Label, jetAK4MassProduct = Handle ("std::vector<float>"), ("jetsAK4"), ("jetAK4Mass")

genPtHandle, genLabel, genPtProd = Handle("std::vector<float>"), ("genPart"), ("genPartPt")
genEtaHandle, genLabel, genEtaProd = Handle("std::vector<float>"), ("genPart"), ("genPartEta")
genPhiHandle, genLabel, genPhiProd = Handle("std::vector<float>"), ("genPart"), ("genPartPhi")
genEHandle, genLabel, genEProd = Handle("std::vector<float>"), ("genPart"), ("genPartE")
genIDHandle, genLabel, genIDProd = Handle("std::vector<float>"), ("genPart"), ("genPartID")
genStatusHandle, genLabel, genStatusProd = Handle("std::vector<float>"), ("genPart"), ("genPartStatus")
genMomIDHandle, genLabel, genMomIDProd = Handle("std::vector<float>"), ("genPart"), ("genPartMomID")

f = ROOT.TFile(options.outname, "RECREATE")
f.cd()

hzetadist = ROOT.TH1D("hzetadist",";#zeta (#deltaR);;",90,0,4.5)
hzetaHdist = ROOT.TH1D("hzetaHdist","#zeta;;",90,0,4.5)
hzetaWdist = ROOT.TH1D("hzetaWdist","#zeta;;",90,0,4.5)
hzetaBackdist = ROOT.TH1D("hzetaBackdist","#zeta;;",90,0,4.5)
i = 1
for event in events:
     if (i>options.maxEvts):
          break
     i=i+1
     event.getByLabel (jetAK4Label, jetAK4PtProduct, jetAK4PtHandle)
     event.getByLabel (jetAK4Label, jetAK4EtaProduct, jetAK4EtaHandle)
     event.getByLabel (jetAK4Label, jetAK4PhiProduct, jetAK4PhiHandle)
     event.getByLabel (jetAK4Label, jetAK4EProduct, jetAK4EHandle)
     event.getByLabel (jetAK4Label, jetAK4MassProduct, jetAK4MassHandle)
     event.getByLabel (genLabel, genPtProd, genPtHandle)
     event.getByLabel (genLabel, genEtaProd, genEtaHandle)
     event.getByLabel (genLabel, genPhiProd, genPhiHandle)
     event.getByLabel (genLabel, genEProd, genEHandle)
     event.getByLabel (genLabel, genIDProd, genIDHandle)
     event.getByLabel (genLabel, genStatusProd, genStatusHandle)
     event.getByLabel (genLabel, genMomIDProd, genMomIDHandle)

     jetAK4Pt = jetAK4PtHandle.product()
     jetAK4Eta = jetAK4EtaHandle.product()
     jetAK4Phi = jetAK4PhiHandle.product()
     jetAK4E = jetAK4EHandle.product()
     jetAK4Mass = jetAK4MassHandle.product()

     genPt = genPtHandle.product()
     genEta = genEtaHandle.product()
     genPhi = genPhiHandle.product()
     genE = genEHandle.product()
     genID = genIDHandle.product()
     genStatus = genStatusHandle.product()
     genMomID = genMomIDHandle.product()

     if jetAK4Pt.size() <= 1:
  	  continue
     iW = -1
     ih = -1

     ####determine if higgs or W in event
     for igen in range(0,genPt.size()):
          if abs(genID.at(igen)) == 25 and abs(genMomID.at(igen))==6000006:
               ih = igen
          if abs(genID.at(igen)) == 24 and abs(genMomID.at(igen)) == 6:
               iW = igen

     ####Creating ordered dictionary to find maximum mass at specific Pt index
     named_map = collections.namedtuple('named_map','mass index')
     map = {}

     for iak4 in range(0, jetAK4Pt.size()):
          map[iak4] = jetAK4Mass.at(iak4)

     #print map.keys()
     #print map.values()
     map_sorted = sorted([named_map(v,k) for (k,v) in map.items()], reverse=True)
     maxMass = map_sorted[0]
     maxMass2 = map_sorted[1]
     #print maxMass.index,maxMass.mass
     #print maxMass2.index,maxMass2.mass
     maxak4 = maxMass.index
     max2ak4 = maxMass2.index
	
     ak4p4 = ROOT.TLorentzVector()
     ak4p4_2 = ROOT.TLorentzVector()
     
     ####Set Lorentz vectors for two highest mass jets
     ak4p4.SetPtEtaPhiE(jetAK4Pt.at(maxak4), jetAK4Eta.at(maxak4), jetAK4Phi.at(maxak4), jetAK4E.at(maxak4))
     ak4p4_2.SetPtEtaPhiE(jetAK4Pt.at(max2ak4), jetAK4Eta.at(max2ak4), jetAK4Phi.at(max2ak4), jetAK4E.at(max2ak4))
     
     ####Set Lorentz vector for invariant dijet mass
     invMassVec = ak4p4 + ak4p4_2
     
     ####scaled mass drop variable
     zeta = jetAK4Mass.at(maxak4) / (invMassVec.M()) * ak4p4.DeltaR(ak4p4_2)
     
     hzetaBackdist.Fill(zeta)
     if ih > -1 or iW > -1:

          hzetadist.Fill(zeta)
          if ih > -1:
               hzetaHdist.Fill(zeta)
          if iW > -1:
               hzetaWdist.Fill(zeta)
f.cd()
#hzetadist.Rebin(1)
f.Write()	
