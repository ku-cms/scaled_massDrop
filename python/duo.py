#! /usr/bin/env python
#Configuration

from optparse import OptionParser
parser = OptionParser()


parser.add_option('--infileSig', type='string', action='store',
			default='outplots.root',
			dest='infile1',
			help='signal input file')

parser.add_option('--infileBack', type='string',action='store',
			default='QCD_1.root',
			dest='infile2',
			help='background input file')

parser.add_option('--outfile', type='string', action='store',
			default='zetadist.root',
			dest='outfile1',
			help='Output file')

(options, args) = parser.parse_args()
argv = []

from ROOT import*

infile1 = TFile(options.infile1)
infile2 = TFile(options.infile2)
outfile = TFile(options.outfile1, 'RECREATE')
gStyle.SetOptStat(0)
hist1 = infile1.Get('hzetadist')
hist2 = infile2.Get('hzetadist')
hist1.Scale(1/hist1.Integral())
hist2.Scale(1/hist2.Integral())
c1 = TCanvas()
hist1.SetLineColor(kRed)
hist2.Draw()
hist1.Draw('SAME')

leg = TLegend(0.6, 0.6, 0.89, 0.89)
leg.SetFillColor(kWhite)
leg.AddEntry(hist1, 'Signal','l')
leg.AddEntry(hist2, 'Background','l')
leg.Draw('SAME')
c1.cd()
c1.Update()
c1.Write()
c1.SaveAs('zetadist_QCD.png')
outfile.Close()

