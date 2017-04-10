import ROOT
import glob


def computeEfficiencies(name, tree):
    print ">> Computing efficiencies for", name
    etaTot = ROOT.TH1F("{}_etaTot".format(name), "etaTot", 30, -3., 3. )
    etaCut = ROOT.TH1F("{}_etaCut".format(name), "etaCut", 30, -3., 3. )
    ptTot = ROOT.TH1F("{}_ptTot".format(name), "ptTot", 9, 5., 50. )
    ptCut = ROOT.TH1F("{}_ptCut".format(name), "ptCut", 9, 5., 50. )
    npvTot = ROOT.TH1F("{}_npvTot".format(name), "npvTot", 50, 0., 50. )
    npvCut = ROOT.TH1F("{}_npvCut".format(name), "npvCut", 50, 0., 50. )

    tree.Draw("genEta>>{}".format(etaTot.GetName()), "abs(genEta)<2.5 && genPt>30 && genPt<40 && nVtx<300")
    tree.Draw("genEta>>{}".format(etaCut.GetName()), "abs(genEta)<2.5 && genPt>30 && genPt<40 && nVtx<300 && isMatched==1")
    tree.Draw("genPt>>{}".format(ptTot.GetName()), "abs(genEta)<2.5   && genPt>30 && genPt<40 && nVtx<300")
    tree.Draw("genPt>>{}".format(ptCut.GetName()), "abs(genEta)<2.5   && genPt>30 && genPt<40 && nVtx<300 && isMatched==1")
    tree.Draw("nVtx>>{}".format(npvTot.GetName()), "abs(genEta)<2.5   && genPt>30 && genPt<40 ")
    tree.Draw("nVtx>>{}".format(npvCut.GetName()), "abs(genEta)<2.5   && genPt>30 && genPt<40  && isMatched==1")

    etaEff =  ROOT.TGraphAsymmErrors()
    etaEff.SetName("{}_etaEff".format(name))
    etaEff.BayesDivide(etaCut, etaTot)
    ptEff  =  ROOT.TGraphAsymmErrors()
    ptEff.SetName("{}_ptEff".format(name))
    ptEff.BayesDivide(ptCut, ptTot)
    npvEff  =  ROOT.TGraphAsymmErrors()
    npvEff.SetName("{}_npvEff".format(name))
    npvEff.BayesDivide(npvCut, npvTot)

    return etaEff, ptEff, npvEff


files = {}
files["PU0_74_RefMatching"] = ROOT.TFile.Open("./test.root")
files["PU0_73_RefMatching"] = ROOT.TFile.Open("../../../../../CMSSW_7_3_6/src/Regression/RegressionTrees/test/test_pu0_referenceMatching.root")
files["RelVal_RefMatching"] = ROOT.TFile.Open("../../../../../CMSSW_7_3_6/src/Regression/RegressionTrees/test/test_relval_referenceMatching.root")

output = ROOT.TFile.Open("efficiencies.root", "RECREATE")
histos = {}
for name,f in files.items():
    tree = f.Get("gedGsfElectronTree/RegressionTree")
    tree.__class__ = ROOT.TTree
    etaEff, ptEff, npvEff = computeEfficiencies(name, tree)
    etaEff.Write()
    ptEff.Write()
    npvEff.Write()
    histos[name] = (etaEff,ptEff,npvEff)


canvas = []
legend = []
canvas.append(ROOT.TCanvas("RelvalVsPU0", "RelvalVsPU0", 800, 800))
histos["RelVal_RefMatching"][0].SetMarkerStyle(20)
histos["RelVal_RefMatching"][0].SetMarkerColor(ROOT.kRed)
histos["RelVal_RefMatching"][0].SetLineColor(ROOT.kRed)
histos["PU0_73_RefMatching"][0].SetMarkerStyle(24)
histos["PU0_73_RefMatching"][0].SetMarkerColor(ROOT.kBlack)
histos["PU0_74_RefMatching"][0].SetMarkerStyle(20)
histos["PU0_74_RefMatching"][0].SetMarkerColor(ROOT.kBlack)

histos["RelVal_RefMatching"][0].GetXaxis().SetTitle("#eta^{gen}")
histos["RelVal_RefMatching"][0].GetYaxis().SetTitle("Efficiency")
histos["RelVal_RefMatching"][0].Draw("ap")
histos["PU0_73_RefMatching"][0].Draw("p same")
histos["PU0_74_RefMatching"][0].Draw("p same")
legend.append(ROOT.TLegend(0.4, 0.15, 0.6, 0.3))
legend[-1].SetFillColor(0)
legend[-1].AddEntry(histos["RelVal_RefMatching"][0], "RelVal", "p")
legend[-1].AddEntry(histos["PU0_73_RefMatching"][0], "73", "p")
legend[-1].AddEntry(histos["PU0_74_RefMatching"][0], "74", "p")
legend[-1].Draw()



for c in canvas:
    c.Print("plots/{}.png".format(c.GetName()))


output.Close()
