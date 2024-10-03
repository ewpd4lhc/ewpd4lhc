#!/usr/bin/env python3

import yaml
import argparse
from ROOT import RooArgSet, RooArgList, RooRealVar, RooFormulaVar, TMatrixDSym, RooDataSet, RooWorkspace, RooGaussian, RooMultiVarGaussian, RooProdPdf, RooFit, TFile
from ROOT.RooStats import ModelConfig


def make_prediction(smprediction, lin, quad):
    param_names = set()
    formula = str(smprediction)
    if not lin is None:
        for c in lin:
            if isinstance(lin[c], tuple) and lin[c][0] != 0. or lin[c] != 0.:
                if isinstance(lin[c], tuple):
                    formula += "+{}*({}-{})".format(lin[c][0], c, lin[c][1])
                else:
                    formula += "+{}*{}".format(lin[c], c)
                param_names.add(c)
    if not quad is None:
        for cs in quad:
            if quad[cs] != 0.:
                name1, name2 = cs.split('*')
                formula += "+{}*{}*{}".format(quad[cs], name1, name2)
                param_names.add(name1)
                param_names.add(name2)
    return sorted(list(param_names)), formula


def make_workspace(outfn, yamldata, rescaling=False):
    # if covariance too small, rescale measurements, to avoid numerical issues
    covMatrix = TMatrixDSym(len(yamldata))
    for i1, o1 in enumerate(yamldata):
        for i2, o2 in enumerate(yamldata):
            covMatrix[i1, i2] = yamldata[o1]['correlation'][o2] * \
                yamldata[o1]['error']*yamldata[o2]['error']

    if rescaling:
        det = covMatrix.Determinant()
        rescale = 1
        for n in range(1, 10):
            if (10**n)**(2*len(yamldata))*det > 10e10:
                break
            else:
                rescale = 10**n
        # rescaling
        for o in yamldata:
            yamldata[o]['smprediction'] *= rescale
            yamldata[o]['measurement'] *= rescale
            yamldata[o]['error'] *= rescale
            if 'theoerr' in yamldata[o]:
                for x in yamldata[o]['theoerr']:
                    yamldata[o]['theoerr'][x] *= rescale
            if 'd6lin' in yamldata[o]:
                for c in yamldata[o]['d6lin']:
                    yamldata[o]['d6lin'][c] *= rescale
            if 'd8lin' in yamldata[o]:
                for c in yamldata[o]['d8lin']:
                    yamldata[o]['d8lin'][c] *= rescale
            if 'd6quad' in yamldata[o]:
                for c in yamldata[o]['d6quad']:
                    yamldata[o]['d6quad'][c] *= rescale

        covMatrix = TMatrixDSym(len(yamldata))
        for i1, o1 in enumerate(yamldata):
            for i2, o2 in enumerate(yamldata):
                if 'correlation' in yamldata[o1] and o2 in yamldata[o1]['correlation']:
                    covMatrix[i1, i2] = yamldata[o1]['correlation'][o2] * \
                        yamldata[o1]['error']*yamldata[o2]['error']
                    covMatrix[i2, i1] = yamldata[o1]['correlation'][o2] * \
                        yamldata[o1]['error']*yamldata[o2]['error']

    predictions = []
    muVec = RooArgList()
    xVec = RooArgList()
    rooparameters = {}
    coefficients = []
    predictions = []
    measurements = []
    for o in yamldata:
        lin = dict(yamldata[o]['d6lin']) if 'd6lin' in yamldata[o] else dict()
        if 'd8lin' in yamldata[o]:
            lin.update(yamldata[o]['d8lin'])
            for c in yamldata[o]['d8lin']:
                if not c in coefficients:
                    coefficients.append(c)
        quad = yamldata[o]['d6quad'] if 'd6quad' in yamldata[o] else dict()
        for c in list(quad.keys())+list(lin.keys()):
            if not c in coefficients:
                coefficients.append(c)

        isSMinput = False
        if 'smdependence' in yamldata[o]:
            for x in yamldata[o]['smdependence']:
                if x == o:
                    if len(yamldata[o]['smdependence']) != 1 or yamldata[o]['smdependence'][x] != 1:
                        raise Exception(
                            'Observable {} depends on itself but not like expected for input parameter'.format(o))
                    isSMinput = True
                else:
                    lin[x] = (yamldata[o]['smdependence'][x],
                              yamldata[x]['smprediction'])
        if 'theoerr' in yamldata[o]:
            lin.update(yamldata[o]['theoerr'])

        if isSMinput:
            if len(lin) != 0 or len(quad) != 0:
                raise Exception(
                    'Observable {} is assumed to be input put has dependence on other parameters'.format(o))
            formula = o
            params = RooArgList()
            if o in rooparameters:
                prediction = rooparameters[o]
                params.add(rooparameters[o])
            else:
                central = yamldata[o]['smprediction']
                up = yamldata[o]['smprediction']+yamldata[o]['error']*100
                down = yamldata[o]['smprediction']-yamldata[o]['error']*100
                prediction = RooRealVar(o, o, central, down, up)
                rooparameters[o] = prediction
                params.add(rooparameters[o])

        else:
            param_names, formula = make_prediction(
                yamldata[o]['smprediction'], lin, quad)
            params = RooArgList()
            for p in param_names:
                if p in rooparameters:
                    params.add(rooparameters[p])
                else:
                    central = 0.
                    up = 100.
                    down = -100.
                    if p in yamldata:
                        central = yamldata[p]['smprediction']
                        up = yamldata[p]['smprediction'] + \
                            yamldata[p]['error']*100
                        down = yamldata[p]['smprediction'] - \
                            yamldata[p]['error']*100
                    c = RooRealVar(p, p, central, down, up)
                    rooparameters[p] = c
                    params.add(c)
            prediction = RooFormulaVar(
                'predicted_{}'.format(o), formula, params)

        predictions.append(prediction)
        muVec.add(prediction)
        x = RooRealVar('measured_{}'.format(
            o), 'measured_{}'.format(o), yamldata[o]['measurement'])
        x.setConstant()
        xVec.add(x)
        measurements.append(x)

    constraints = []
    global_obs = []
    for p in rooparameters:
        if 'theoerr' in p:
            go = RooRealVar('centralvalue_'+p, 'centralvalue_'+p, 0, -10, 10)
            go.setConstant()
            global_obs.append(go)
            constraint = RooGaussian(
                "constraint_"+p, "constraint_"+p, rooparameters[p], go, RooFit.RooConst(1.))
            constraints.append(constraint)

    if len(constraints) == 0:
        pdf = RooMultiVarGaussian('model', 'model', xVec, muVec, covMatrix)
    else:
        mvg = RooMultiVarGaussian(
            'gaussian', 'gaussian', xVec, muVec, covMatrix)
        pdf = RooProdPdf('model', 'model', [mvg]+constraints)
    data = RooArgSet()
    data.add(xVec)
    dataset = RooDataSet('data', 'data', data)
    dataset.add(RooArgSet(data))

    workspace = RooWorkspace('w')
    getattr(workspace, 'import')(pdf)
    getattr(workspace, 'import')(dataset)
    modelconfig = ModelConfig(workspace)
    modelconfig.SetName('modelconfig')
    modelconfig.SetPdf(pdf)

    pois = RooArgSet()
    nps = RooArgSet()
    for c in rooparameters:
        if c in coefficients:
            pois.add(rooparameters[c])
        else:
            nps.add(rooparameters[c])

    modelconfig.SetParametersOfInterest(pois)
    modelconfig.SetObservables(RooArgSet(data))
    modelconfig.SetNuisanceParameters(nps)
    modelconfig.SetGlobalObservables(global_obs)

    outf = TFile(outfn, 'RECREATE')
    workspace.Write()
    modelconfig.Write()
    outf.Close()


def main():
    parser = argparse.ArgumentParser(
        description='Creates RooFit model from yaml file.')
    parser.add_argument('--input', help='path to yaml file',
                        default='ewpo_out.yml')
    parser.add_argument('--output', help='output file name',
                        default='ewpo_out.root')
    parser.add_argument('--rescale', action='store_true',
                        help='rescale all variables to improve numerical stability', default=False)
    args = parser.parse_args()
    with open(args.input) as f:
        yamldata = yaml.safe_load(f)
        make_workspace(args.output, yamldata, rescaling=args.rescale)


if __name__ == '__main__':
    main()
