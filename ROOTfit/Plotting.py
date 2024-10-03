import ROOT
import math
import numpy
ROOT.gROOT.SetBatch(True)

CANVASXSIZE = CANVASYSIZE = 1200
CANVASMARGINS_TLBR = (0.04, 0.12, 0.12, 0.04)
CANVASMARGINS2D_TLBR = (0.14, 0.12, 0.14, 0.16)

MAXDELTANLLSHOWN = 100
TITLEOFFSETX = 1.0
TITLEOFFSETY = 1.2
TITLEOFFSETZ = None
LABELSIZE = 0.035
TITLESIZE = 0.04
TEXTSIZE = 0.04

WILSON_MAPPING = {
    'cHWB': 'c_{HWB}',
    'cHDD': 'c_{HDD}',
    'cHQ1': 'c_{HQ}^{(1)}',
    'cHQ3': 'c_{HQ}^{(3)}',
    'cHbq': 'c_{Hb}',
    'cHd': 'c_{Hd}',
    'cHe11': 'c_{He,11}',
    'cHe33': 'c_{He,22}',
    'cHe33': 'c_{He,33}',
    'cHj1': 'c_{Hq}^{(1)}',
    'cHj3': 'c_{Hq}^{(3)}',
    'cHl111': 'c_{Hl,11}^{(1)}',
    'cHl122': 'c_{Hl,22}^{(1)}',
    'cHl133': 'c_{Hl,33}^{(1)}',
    'cHl311': 'c_{Hl,11}^{(3)}',
    'cHl322': 'c_{Hl,22}^{(3)}',
    'cHl333': 'c_{Hl,33}^{(3)}',
    'cHu': 'c_{Hu}',
    'cll1221': 'c_{ll,1221}',
}


def graphProfile1D(name, result):
    result = sorted(result)
    graph = ROOT.TGraph()
    result = [r for r in result if not math.isnan(
        r[1]) and not math.isinf(r[1])]
    for i, r in enumerate(result):
        y = r[1]
        # The graphs shows 2DeltaNLL
        y *= 2
        graph.SetPoint(i, r[0], y)
    graph.SetName(name)
    return graph


def graphProfile2D(name, result):
    graph = ROOT.TGraph2D()
    result = [r for r in result if not math.isnan(
        r[2]) and not math.isinf(r[2])]
    for i, r in enumerate(result):
        # The graphs shows 2DeltaNLL
        z = r[2]
        z *= 2
        graph.SetPoint(i, r[0], r[1], z)
    graph.SetName(name)
    return graph


def plotProfile1D(result, name, outfile, xtitle, ytitle='#minus 2#Delta log(L)', linecolor=ROOT.kRed+1, linewidth=2, linestyle=1, conts2d=False, pdfout=None):

    graph = graphProfile1D('graph_'+name, result)
    outfile.cd()
    graph.Write()

    if xtitle in WILSON_MAPPING:
        xtitle = WILSON_MAPPING[xtitle]

    if conts2d:
        cl1 = 2.297
        cl2 = 5.991
    else:
        cl1 = 1
        cl2 = 3.841
    c = ROOT.TCanvas('canvas_'+name, name, CANVASXSIZE, CANVASYSIZE)
    c.SetTopMargin(CANVASMARGINS_TLBR[0])
    c.SetLeftMargin(CANVASMARGINS_TLBR[1])
    c.SetBottomMargin(CANVASMARGINS_TLBR[2])
    c.SetRightMargin(CANVASMARGINS_TLBR[3])
    c.SetTicks(1, 1)
    xs = numpy.array(graph.GetX())
    ys = numpy.array(graph.GetY())
    imin = -1
    ymin = 1e9
    ymax = -1
    imax = -1
    x0 = 0.
    x1 = 0.
    for i, y in enumerate(ys):
        if y < ymin:
            ymin = y
            imin = i
            x0 = xs[i]
        if y > ymax:
            ymax = y
            imax = i
            x1 = xs[i]
    print(imin, ys[:imin], xs[:imin])
    g_inverse1 = ROOT.TGraph(imin, ys[:imin], xs[:imin])
    g_inverse2 = ROOT.TGraph(graph.GetN()-imin, ys[imin:], xs[imin:])
    ci1 = g_inverse1.Eval(cl1), g_inverse2.Eval(cl1)
    ci2 = g_inverse1.Eval(cl2), g_inverse2.Eval(cl2)
    graph.SetLineColor(linecolor)
    graph.SetLineStyle(linestyle)
    graph.SetLineWidth(linewidth)
    graph.Draw('AL')
    graph.GetHistogram().GetYaxis().SetRangeUser(
        0., min(MAXDELTANLLSHOWN, graph.GetHistogram().GetYaxis().GetXmax()))
    graph.GetHistogram().GetXaxis().SetTitleSize(TITLESIZE)
    graph.GetHistogram().GetYaxis().SetTitleSize(TITLESIZE)
    graph.GetHistogram().GetXaxis().SetLabelSize(LABELSIZE)
    graph.GetHistogram().GetYaxis().SetLabelSize(LABELSIZE)
    if TITLEOFFSETX is not None:
        graph.GetHistogram().GetXaxis().SetTitleOffset(TITLEOFFSETX)
    if TITLEOFFSETY is not None:
        graph.GetHistogram().GetYaxis().SetTitleOffset(TITLEOFFSETY)
    graph.GetHistogram().SetLineColor(ROOT.kBlack)
    graph.GetHistogram().GetXaxis().SetTitle(xtitle)
    graph.GetHistogram().GetYaxis().SetTitle(ytitle)

    xmax = graph.GetHistogram().GetXaxis().GetXmax()
    xmin = graph.GetHistogram().GetXaxis().GetXmin()
    cl1line = ROOT.TLine(xmin, cl1, xmax, cl1)
    cl2line = ROOT.TLine(xmin, cl2, xmax, cl2)
    cl1line.Draw()
    cl2line.Draw()
    ci2line1 = ROOT.TLine(ci2[0], 0, ci2[0], cl2)
    ci2line1.SetLineStyle(2)
    ci2line2 = ROOT.TLine(ci2[1], 0, ci2[1], cl2)
    ci2line2.SetLineStyle(2)
    ci1line1 = ROOT.TLine(ci1[0], 0, ci1[0], cl1)
    ci1line1.SetLineStyle(2)
    ci1line2 = ROOT.TLine(ci1[1], 0, ci1[1], cl1)
    ci1line2.SetLineStyle(2)
    err_low = abs(ci1[0]-x0)
    err_high = abs(ci1[1]-x0)
    rnd = 2-int(math.log10(err_low/2+err_high/2))
    central = round(x0, rnd)
    ci1text = ROOT.TLatex(0.5, 0.85, 'Best fit: {}_{{ #minus{} }}^{{ #plus{} }}'.format(
        central, round(err_low, rnd), round(err_high, rnd)))
    ci1text.SetNDC()
    ci1text.SetTextSize(TEXTSIZE)
    ci1text.SetTextFont(42)
    ci1text.SetTextAlign(21)
    ci1text.Draw()
    ci2text = ROOT.TLatex(0.5, 0.79, '95% CL: [{},{}]'.format(
        round(ci2[0], rnd), round(ci2[1], rnd)))
    ci2text.SetNDC()
    ci2text.SetTextSize(TEXTSIZE)
    ci2text.SetTextFont(42)
    ci2text.SetTextAlign(21)
    ci1line1.Draw()
    ci1line2.Draw()
    ci2line1.Draw()
    ci2line2.Draw()
    ci2text.Draw()
    printAndSave(c, outfile, pdfout)


def plotProfile2D(result, name, outfile, xtitle, ytitle, ztitle='#minus 2#Delta log(L)', linecolor=ROOT.kBlack, outfolder='.', isNLL=True, conts1d=False, toy_result=None, linecolor_toys=ROOT.kRed+1, pdfout=None, colors=True):
    ROOT.gStyle.SetPalette(ROOT.kThermometer)
    ROOT.TColor.InvertPalette()

    graph = graphProfile2D('graph_'+name, result)

    if xtitle in WILSON_MAPPING:
        xtitle = WILSON_MAPPING[xtitle]
    if ytitle in WILSON_MAPPING:
        ytitle = WILSON_MAPPING[ytitle]

    if conts1d:
        cl1 = 1
        cl2 = 3.841
    else:
        cl1 = 2.297
        cl2 = 5.991
    c = ROOT.TCanvas('canvas_'+name, name, 1024, 1024)
    c.SetTopMargin(CANVASMARGINS2D_TLBR[0])
    c.SetLeftMargin(CANVASMARGINS2D_TLBR[1])
    c.SetBottomMargin(CANVASMARGINS2D_TLBR[2])
    c.SetRightMargin(CANVASMARGINS2D_TLBR[3])

    c.SetTicks(1, 1)

    graph.GetHistogram().SetTitle('')
    graph.GetHistogram().GetXaxis().SetTitle(xtitle)
    graph.GetHistogram().GetYaxis().SetTitle(ytitle)
    graph.GetHistogram().GetZaxis().SetTitle(ztitle)
    graph.GetHistogram().GetXaxis().SetTitleSize(TITLESIZE)
    graph.GetHistogram().GetYaxis().SetTitleSize(TITLESIZE)
    graph.GetHistogram().GetZaxis().SetTitleSize(TITLESIZE)
    graph.GetHistogram().GetXaxis().SetLabelSize(LABELSIZE)
    graph.GetHistogram().GetYaxis().SetLabelSize(LABELSIZE)
    graph.GetHistogram().GetZaxis().SetLabelSize(LABELSIZE)
    if TITLEOFFSETX is not None:
        graph.GetHistogram().GetXaxis().SetTitleOffset(TITLEOFFSETX)
    if TITLEOFFSETY is not None:
        graph.GetHistogram().GetYaxis().SetTitleOffset(TITLEOFFSETY)
    if TITLEOFFSETZ is not None:
        graph.GetHistogram().GetZaxis().SetTitleOffset(TITLEOFFSETZ)
    if colors:
        red = numpy.array([0.1, 1.])
        green = numpy.array([0.1, 1.])
        blue = numpy.array([1, 1.])
        stops = numpy.array([0.0, 1.00])
        ROOT.TColor.CreateGradientColorTable(2, stops, red, green, blue, 1001)

        pad1 = ROOT.TPad("pad1", "", 0, 0, 1, 1)
        pad2 = ROOT.TPad("pad2", "", 0, 0, 1, 1)
        for p in pad1, pad2:
            p.SetRightMargin(0.16)
            p.SetTopMargin(0.12)
            p.SetLeftMargin(0.1)
            p.SetBottomMargin(0.12)
        pad2.SetFillStyle(0)
        pad2.SetFillColor(0)
        pad2.SetFrameFillStyle(0)
        pad1.Draw()
        pad1.cd()
        graph.GetHistogram().GetZaxis().SetRangeUser(0, 10)
        h0 = graph.GetHistogram().Clone()
        h0.SetContour(20)
        h0.SetLineColor(linecolor)
        h0.SetLineWidth(2)
        h1 = h0.Clone()
        h0.SetContour(10)
        h0.Draw('cont4z')
        pad1.Update()
        pad1.Modified()
        c.cd()
        pad2.Draw()
        pad2.cd()
    else:
        h1 = graph.GetHistogram().Clone()
        c.cd()
    h1.SetContour(1)
    h1.SetContourLevel(0, cl1)
    h1.SetLineColor(linecolor)
    h1.SetLineWidth(3)
    h1.SetLineStyle(2)
    h1.Draw('same cont3')
    h2 = graph.GetHistogram().Clone('h2')
    h2.SetContour(1)
    h2.SetContourLevel(0, cl2)
    h2.SetLineColor(linecolor)
    h2.SetLineWidth(3)
    h2.Draw('cont3 same')
    c.Draw()
    printAndSave(c, outfile, pdfout)


def printAndSave(c, outfile, pdfout):
    outfile.cd()
    c.Write()
    if pdfout is not None:
        c.SaveAs(pdfout)
