#!/usr/bin/env python3

import ROOT
import argparse
import itertools
import math
import os
import yaml
import Plotting


def main(infilename, outfolder, actual_poi_names, float_par_names, scan_range, npoints, workspace_name, modelconfig_name, data_name, doplots):
    """
    Main function to fit a model to data and optionally scan over a range of POIs (parameters of interest).
    The results can be saved in a YAML file and plotted if desired.

    Args:
        infilename (str): Path to the input ROOT file containing the workspace.
        outfolder (str): Path to the output folder where results will be saved.
        actual_poi_names (list): List of names of POIs (parameters of interest) to fit.
        float_par_names (list): List of parameters to float during the fit.
        scan_range (list): List of ranges to scan for each POI.
        npoints (int): Number of points to scan in each POI range.
        workspace_name (str): Name of the workspace inside the ROOT file (search for single workspace if None).
        modelconfig_name (str): Name of the model configuration inside the ROOT file (search for single workspace if None).
        data_name (str): Name of the dataset inside the workspace.
        doplots (bool): Whether to create plots from the scan results.

    """

    # Load model, data, and parameters from the input file
    pdf, data, pois, nps = load_data(
        infilename, workspace_name, modelconfig_name, data_name)

    # Combine POIs with float parameters, and set unused parameters as constant
    if float_par_names != 'all' and actual_poi_names != 'all':
        float_par_names = list(set(float_par_names + actual_poi_names))
        for p in pois:
            if not p.GetName() in float_par_names:
                p.setConstant()
                
    # Perform fit in case of user-defined POIs or if there are nuisance parameters
    if len(actual_poi_names) > 0 or len(nps)>0:
        actual_pois = ROOT.RooArgSet()
        for pn in actual_poi_names:
            for p in pois:
                if p.GetName() == pn:
                    actual_pois.add(p)
        assert len(actual_poi_names) == len(actual_pois)

        fitresult = pdf.fitTo(data, Save=True, PrintLevel=0,
                              Strategy=2, Minos=actual_pois)
        if fitresult.status() != 0:
            print('WARNING: fit did not converge (possibly degenerate?)')
            return

    # Perform scan over scan ranges if specified
    if len(scan_range) > 0:
        if len(actual_poi_names) != len(scan_range):
            raise Exception("Need to give a scan range for each (or no) POI!")

        res = scan(pdf, data, actual_pois, scan_range, npoints)

        # Save scan results to a YAML file
        outname = 'scan_2DeltaNLL_' + '_'.join(actual_poi_names)
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        with open(os.path.join(outfolder, outname + '.yml'), 'w') as f:
            print('Writing scan results to', os.path.join(
                outfolder, outname + '.yml'))
            res_dictlist = []
            for r in res:
                rdict = {}
                for ia, a in enumerate(actual_poi_names):
                    rdict[a] = r[ia]
                rdict['nll'] = r[-1]
                res_dictlist.append(rdict)
            yaml.dump(res_dictlist, f)

        # Create plots if required
        if doplots:
            rootoutfile = ROOT.TFile(os.path.join(
                outfolder, 'plots.root'), 'recreate')
            pdfout = os.path.join(outfolder, 'scan_2DeltaNLL_{}{}'.format(
                '_'.join(actual_poi_names), '.pdf'))
            if len(actual_poi_names) == 1:
                Plotting.plotProfile1D(
                    res, outname, rootoutfile, actual_poi_names[0], pdfout=pdfout)
            elif len(actual_poi_names) == 2:
                Plotting.plotProfile2D(
                    res, outname, rootoutfile, actual_poi_names[0], actual_poi_names[1], pdfout=pdfout)
            print('Saving canvas as', pdfout)
            rootoutfile.Close()


def load_data(infilename, workspace_name, modelconfig_name, data_name):
    """
    Load the workspace, model configuration, and data from the input ROOT file.

    Args:
        infilename (str): Path to the input ROOT file.
        workspace_name (str): Name of the workspace (None to autodetect).
        modelconfig_name (str): Name of the model configuration (None to autodetect).
        data_name (str): Name of the dataset (None to autodetect).

    Returns:
        tuple: Returns the pdf (probability density function), data, parameters of interest (pois), and nuisance parameters (nps).

    Raises:
        Exception: If multiple or no datasets are found.
    """
    f = ROOT.TFile(infilename)
    workspace = modelconfig = None

    # Load model configuration and workspace
    if modelconfig_name is None:
        modelconfig = getObjFromFile(f, ROOT.RooStats.ModelConfig)
    else:
        modelconfig = f.Get(modelconfig_name)

    if workspace_name is None:
        workspace = getObjFromFile(f, ROOT.RooWorkspace)
    else:
        workspace = f.Get(workspace_name)

    workspace.Print()

    # Load dataset
    if data_name is None:
        datasets = workspace.allData()
        if len(datasets) == 1:
            data = datasets.back()
        elif len(datasets) == 0:
            raise Exception("No datasets found in file")
        elif len(datasets) > 1:
            raise Exception(
                "More than one datasets found in file, please give dataset name")
    else:
        data = workspace[data_name]

    # Extract PDF and parameters
    pdf = modelconfig.GetPdf()
    nps = modelconfig.GetNuisanceParameters()
    pois = modelconfig.GetParametersOfInterest()
    f.Close()
    return pdf, data, pois, nps


def getObjFromFile(f, objtype):
    """
    Retrieve an object of a given type from a ROOT file.
    """
    obj = None
    for key in f.GetListOfKeys():
        if isinstance(key.ReadObj(), objtype):
            if obj is not None:
                raise Exception(
                    "More than one object of type found in file, please specify the name")
            obj = key.ReadObj()
    return obj


def split_arg(arg):
    """
    Parse a comma-separated string argument into a list.
    Secondary separation by ":" will create list of lists.
    'all' will be return as 'all'
    """
    if arg == '':
        return []
    elif arg == 'all':
        return 'all'
    return [[y.strip() for y in x.strip().split(':')] if ':' in x else x.strip() for x in arg.split(',')]


def scan(pdf, data, pois, scan_ranges, nsteps, inclOffset=False, print_level=-1):
    """
    Perform a scan over a set of parameters of interest (POIs) and compute the negative log-likelihood (NLL) at each point.

    Args:
        pdf (ROOT.RooAbsPdf): The probability density function to fit.
        data (ROOT.RooAbsData): The dataset to fit.
        pois (ROOT.RooArgSet): The parameters of interest to scan.
        scan_ranges (list): List of scan ranges for each POI.
        nsteps (int): Number of points to scan for each POI.
        inclOffset (bool): Whether to include offsetting in the NLL.
        print_level (int): Level of verbosity for printing scan information.

    Returns:
        list: List of tuples containing the scan points and their corresponding NLL values.

    """
    nsteps = [nsteps] * len(pois)
    assert len(scan_ranges) == len(pois)
    assert len(nsteps) == len(pois)

    for sr in scan_ranges:
        assert len(sr) == 2

    scan_ranges = [(float(sr[0]), float(sr[1])) for sr in scan_ranges]
    print("Performing scan of ranges:")
    for poi, sr in zip(pois, scan_ranges):
        print(poi.GetName(), sr)

    scan_points = []
    for sr, ns in zip(scan_ranges, nsteps):
        scan_points.append(
            [sr[0] + i * float(sr[1] - sr[0]) / (ns - 1) for i in range(ns)])

    grid = list(itertools.product(*scan_points))
    result = []

    # Check if profiling is necessary
    do_profiling = False
    for par in pdf.getParameters(data):
        if not par.isConstant() and not par in pois:
            do_profiling = True
            break

    # Create NLL (negative log-likelihood) function
    nll = pdf.createNLL(data, ROOT.RooFit.Verbose(False))
    minimizer = ROOT.RooMinimizer(nll)
    minimizer.setPrintLevel(-1)
    minimizer.setStrategy(2)
    minimizer.setErrorLevel(0.5)
    minimizer.setOffsetting(True)

    # Perform an initial minimization to find the global minimum of the NLL
    status = minimizer.migrad()
    if status != 0:
        raise Exception("Bad initial fit!")

    global_minNLL = nll.getVal()
    global_minNLL_point = tuple([p.getVal() for p in pois])

    # Scan the grid and evaluate NLL at each point
    for ipoint, gridpoint in enumerate(grid):
        for poi, val in zip(pois, gridpoint):
            poi.setVal(val)
            poi.setConstant()

        # Minimize NLL wrt to floating pars
        if do_profiling:
            minimizer.setPrintLevel(print_level)
            status = minimizer.migrad()
            if status != 0:
                print(
                    "WARNING: bad scan at {}, skipping and resetting likelihood".format(gridpoint))
                for poi, val in zip(pois, gridpoint):
                    poi.setVal(0)
                    poi.setConstant(False)
                    minimizer.migrad()
                continue

        minNLL = nll.getVal()
        if not inclOffset:
            minNLL = minNLL - global_minNLL
            if minNLL < -0.2:
                raise Exception(
                    "Bad scan, found global minimum better than best fit at " + ','.join([str(x) for x in gridpoint]))

        result.append(gridpoint + (minNLL,))
        if print_level >= 0:
            print(gridpoint + (minNLL,))

    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Fit and scan parameters in a ROOT workspace')
    parser.add_argument(
        '-i', '--input', help='Path to ROOT file containing workspace')
    parser.add_argument(
        '-p', '--pois', help='Parameters to fit, comma separated', default='')
    parser.add_argument(
        '-f', '--float', help='Parameters to float in the fit, comma separated.', default='')
    parser.add_argument('-w', '--workspacename',
                        help='Name of workspace.', default=None)
    parser.add_argument('-d', '--dataname',
                        help='Name of dataset.', default=None)
    parser.add_argument('-m', '--modelconfigname',
                        help='Name of modelconfig.', default=None)
    parser.add_argument(
        '-s', '--scan', help='Parameter ranges to scan for POIs, comma separated (e.g -3:3,1:2)', default='')
    parser.add_argument('-n', '--npoints',
                        help='Number of points in scan', default=20)
    parser.add_argument('-o', '--outfolder',
                        help='Output folder', default='rootfit')
    parser.add_argument('--plots', help='Create plots?',
                        action='store_true', default=False)

    args = parser.parse_args()

    # Call the main function with the parsed arguments
    main(
        infilename=args.input,
        outfolder=args.outfolder,
        actual_poi_names=split_arg(args.pois),
        float_par_names=split_arg(args.float),
        scan_range=split_arg(args.scan),
        npoints=int(args.npoints),
        workspace_name=args.workspacename,
        modelconfig_name=args.modelconfigname,
        data_name=args.dataname,
        doplots=args.plots
    )
