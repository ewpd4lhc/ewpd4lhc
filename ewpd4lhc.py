#!/usr/bin/env python3

import configparser
import argparse
import os
import yaml
from math import floor, log10
import SMcalculator
import SMEFTlikelihood
import Utils
import LINAfit
import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# Allowed uncertainties (add to covariance, introduce nuisance parameter, ignore)
SUPPORTED_ERROR_TYPES = ['covariance', 'nuispar', 'off']
SIN2THETAEEFF_DERIVED=['Ae','AFBe','sin2thetaeeff']
SIN2THETALEFF_DERIVED=['Al','AFBl','sin2thetaleff']


def main():
    """
    Main function to create a  SMEFT EWPD likelihood in a YAML file and optionally a ROOFit workspace.

    This function handles command-line arguments, reads configuration files, and initializes the necessary parameters for the likelihood calculation.
    """
    parser = argparse.ArgumentParser(
        description='Creates description of SMEFT EWPD likelihood in yaml file.')
    parser.add_argument('--config', help='path to config file',
                        default='config/ewpd4lhc.cfg')
    parser.add_argument('--data', help='path to data file', default=None)
    parser.add_argument(
        '--output', help='path to output yaml file', default='ewpd_out.yml')
    parser.add_argument(
        '--root_output', help='path to output root file', default=None)
    parser.add_argument('--lin', action='store_true',
                        help='only linear dim6 EFT parametrization', default=False)
    parser.add_argument('--smonly', action='store_true',
                        help='no EFT parametrization', default=False)
    parser.add_argument(
        '--decimals', help='Number of decimals stored in yaml file', default=12)
    parser.add_argument('--norescale_root', action='store_true',
                        help='By default, measurements and predictions in root workspace are rescaled such that absolute value likelihood becomes larger, improving numerics. This turns this off.', default=False)

    # Parse command-line arguments
    args = parser.parse_args()
    outfn = args.output
    smonly = args.smonly
    linonly = args.lin
    root_output = args.root_output
    ndigits = int(args.decimals)
    rescale_root = not args.norescale_root

    # Parse additional command line arguments and read config file
    scheme, symmetry, observables, sminputs, coefficients, Lambda, d6lin_type, d6quad_type, d8lin_type, theo_err_treatment, para_err_treatment, datapath = readCfg(
        args)

    if 'sin2theta' in scheme:
        lfu = symmetry in SMEFTlikelihood.LFU_SYMMETRIES
        measurement_names, measurement_central, measurement_covariance, observable_mapping, fixed_prediction = Utils.load_data(datapath + '.yml')
        sin2thetaleff,sin2thetaleff_err = determine_sin2thetaleff(measurement_names, measurement_central, measurement_covariance, observable_mapping, scheme, lfu)
        if lfu:
            observables=[o for o in observables if not observable_mapping[o] in SIN2THETALEFF_DERIVED]
        else:
            observables=[o for o in observables if not observable_mapping[o] in SIN2THETAEEFF_DERIVED]
        dataupdates={'sin2thetaleff':{'central':sin2thetaleff,'error':sin2thetaleff_err,'prediction':'sin2thetaleff'}}
    else:
        dataupdates=None
    

    # Initialize the EWPD likelihood object
    ewpd = SMEFTlikelihood.EWPDlikelihood(
        scheme=scheme,
        symmetry=symmetry,
        observables=observables,
        Lambda=Lambda,
        d6linSource=d6lin_type,
        d6quadSource=d6quad_type,
        d8linSource=d8lin_type,
        datapath=datapath + '.yml',
        parapath=datapath + '_SMEFTpara/',
        dataupdates=dataupdates
    )

    # Return normalized observables for better numerics (not tested)
    normalize = False

    # Get measurements and predictions
    meas = ewpd.measured(incl_inputs=True, normalized=normalize)
    names = list(meas.keys())
    pred = ewpd.predicted(incl_inputs=True, normalized=normalize)

    # Get covariance and correlation matrices of measurements alone
    meas_cov = ewpd.covariance(
        incl_inputs=True, with_para_err=False, with_theo_err=False, normalized=normalize)
    # Get SM parametrization (and include theory uncertainties in parameter
    smparas, smpara = ewpd.linSMpara(
        incl_inputs=True, with_theo_err=theo_err_treatment.lower() != 'off', normalized=normalize)

    # Full error and covariance calculations
    full_err = ewpd.error(
        incl_inputs=para_err_treatment == 'nuispar',
        with_para_err=para_err_treatment == 'covariance',
        with_theo_err=theo_err_treatment == 'covariance',
        normalized=normalize
    )
    full_cor = ewpd.correlation(
        incl_inputs=para_err_treatment == 'nuispar',
        with_para_err=para_err_treatment == 'covariance',
        with_theo_err=theo_err_treatment == 'covariance'
    )
    full_cov = ewpd.covariance(
        incl_inputs=para_err_treatment == 'nuispar',
        with_para_err=para_err_treatment == 'covariance',
        with_theo_err=theo_err_treatment == 'covariance',
        normalized=normalize
    )
    # Filter Wilson coefficients based on user input
    coefficients_d6 = [c for c in ewpd.wilson_coefficients_d6() if coefficients is None or c in coefficients]
    coefficients_d8 = [c for c in ewpd.wilson_coefficients_d8() if coefficients is None or c in coefficients]

    # d6 linear, quadratic, and d8 contributions if specified
    d6lin = ewpd.d6lin(normalized=normalize) if d6lin_type != 'none' else None
    d6quad = ewpd.d6quad(
        normalized=normalize) if d6quad_type != 'none' else None
    d8lin = ewpd.d8lin(normalized=normalize) if d8lin_type != 'none' else None

    # Prepare all observables for output
    all_observables = [o for o in observables if o not in sminputs] + sminputs

    # Initialize dictionaries to store theoretical error and parameter dependence
    theoerr_para = {}
    smdependence_para = {}
    for o in all_observables:
        theoerr_para[o] = {}
        smdependence_para[o] = {}
        for x in smpara[o]:
            if 'theoerr' in x and theo_err_treatment == 'nuispar':
                theoerr_para[o][x] = smpara[o][x]
            if 'theoerr' not in x and para_err_treatment == 'nuispar':
                smdependence_para[o][x] = smpara[o][x]

    # Prepare output data
    out = {}
    for o in all_observables:
        if o in sminputs and (para_err_treatment == 'covariance' or para_err_treatment == 'off'):
            continue
        out[o] = {}
        out[o]['measurement'] = round(meas[o], ndigits)
        out[o]['smprediction'] = round(pred[o], ndigits)
        out[o]['error'] = round(full_err[o], ndigits)

        # Include d6lin, d6quad, d8lin contributions if available
        if d6lin_type is not None and o in d6lin:
            out[o]['d6lin'] = {k: round(v, ndigits) for k, v in d6lin[o].items(
            ) if coefficients is None or k in coefficients}
        if d6quad_type is not None and o in d6quad:
            out[o]['d6quad'] = {k: round(v, ndigits) for k, v in d6quad[o].items(
            ) if coefficients is None or k.split('*')[0] in coefficients and k.split('*')[1] in coefficients}
        if d8lin_type is not None and o in d8lin:
            out[o]['d8lin'] = {k: round(v, ndigits) for k, v in d8lin[o].items(
            ) if coefficients is None or k in coefficients}

        # Round correlations for output
        corr = full_cor[o]
        for o1 in corr:
            corr[o1] = round(corr[o1], ndigits)
        out[o]['correlation'] = corr

        # Include theoretical error and parameter dependence if applicable
        if theo_err_treatment == 'nuispar':
            out[o]['theoerr'] = {}
        if para_err_treatment == 'nuispar':
            out[o]['smdependence'] = {}

        for p in theoerr_para[o]:
            out[o]['theoerr'][p] = round(theoerr_para[o][p], ndigits)
        for p in smdependence_para[o]:
            out[o]['smdependence'][p] = round(smdependence_para[o][p], ndigits)

    runSMfit(names, sminputs, meas, meas_cov, pred, smparas, smpara)

    if not args.smonly and d6lin_type is not None:
        eft_para_fit={}
        for o in d6lin:
            eft_para_fit[o]={}
            for c in coefficients_d6:
                eft_para_fit[o][c] = d6lin[o][c] if c in d6lin[o] else 0
        runSMEFTfit(names, meas, meas_cov, full_cov, pred, 
                    smdependence_para, theoerr_para, coefficients_d6, eft_para_fit)

    printModel(out, theo_err_treatment, para_err_treatment,
               sminputs, d6lin_type, d6quad_type, d8lin_type)

    # Write output to YAML file
    if '/' not in outfn:
        path = os.path.dirname(__file__)
    else:
        path = ''
    if not outfn.lower().endswith('.yml') and not outfn.lower().endswith('.yaml'):
        outfn += '.yml'

    print('Writing output to', os.path.join(path, outfn))
    with open(os.path.join(path, outfn), 'w') as f:
        yaml.dump(out, f)

    if root_output is not None:
        print('Writing output to: '+os.path.join(path, root_output))
        import yaml2root
        yaml2root.make_workspace(root_output, out, rescale_root)


def readCfg(args):
    """
    Reads the configuration file and command-line arguments for setting up the EWPD fit.
    """
    # Set up data path (in the directory of this file if not specified in command-line)
    path = os.path.dirname(__file__)
    datapath = os.path.join(path, 'data') if args.data is None else args.data

    cfg = args.config

    # Read configuration file
    config = configparser.ConfigParser()
    if os.path.exists(cfg):
        config.read(cfg)
    elif os.path.exists(os.path.join(path, cfg)):
        config.read(os.path.join(path, cfg))
    else:
        raise Exception(f"Did not find config path {cfg}")

    # arguments from config file
    observables = config['General']['observables'].split()
    Lambda = float(config['SMEFT']['Lambda'])
    coefficients = config['SMEFT'].get('coefficients', None)

    # Determine EFT parameterization sources (e.g. SMEFTsim or various papers potentially including higher order corrections)
    if args.smonly:
        d6lin_type = d6quad_type = d8lin_type = None
    elif args.lin:
        d6lin_type = config['SMEFT'].get(
            'd6lin', None) if config['SMEFT']['d6lin'].lower() != 'none' else None
        d6quad_type = d8lin_type = None
    else:
        d6lin_type = config['SMEFT'].get(
            'd6lin', None) if config['SMEFT']['d6lin'].lower() != 'none' else None
        d6quad_type = config['SMEFT'].get(
            'd6quad', None) if config['SMEFT']['d6quad'].lower() != 'none' else None
        d8lin_type = config['SMEFT'].get(
            'd8lin', None) if config['SMEFT']['d8lin'].lower() != 'none' else None

    # How to handle uncertainties
    para_err_treatment = config['Uncertainties']['parametric'].lower()
    theo_err_treatment = config['Uncertainties']['theory'].lower()
    if para_err_treatment not in SUPPORTED_ERROR_TYPES:
        raise Exception(f"Error {para_err_treatment} not supported")
    if theo_err_treatment not in SUPPORTED_ERROR_TYPES:
        raise Exception(f"Error {theo_err_treatment} not supported")

    # Determine the input scheme for SM calculations: MW,MZ,GF or alpha,MZ,GF
    scheme = config['General']['inputscheme']
    if not hasattr(SMcalculator.INPUTSCHEME, scheme):
        raise Exception(f"Unknown scheme {scheme}")
    sminputs = getattr(SMcalculator.INPUTSCHEME, scheme).value

    # Validate symmetry
    symmetry = config['SMEFT']['symmetry']
    if symmetry not in SMEFTlikelihood.SUPPORTED_SYMMETRIES:
        raise Exception(f"Symmetry {symmetry} not supported")

    # Adjust data path based on symmetry
    if symmetry in SMEFTlikelihood.LFU_SYMMETRIES:
        datapath = os.path.join(datapath, 'pole_observables_lfu')
    else:
        datapath = os.path.join(datapath, 'pole_observables_nolfu')

    return scheme, symmetry, observables, sminputs, coefficients, Lambda, d6lin_type, d6quad_type, d8lin_type, theo_err_treatment, para_err_treatment, datapath

def determine_sin2thetaleff(measurement_names, measurement_central, measurement_covariance, observable_mapping, scheme, lfu):
    if lfu:
        sin2thetaleff_obs=SIN2THETALEFF_DERIVED
    else:
        sin2thetaleff_obs=SIN2THETAEEFF_DERIVED
    s2t_key='s2t'

    sminputs=dict(measurement_central)
    sminputs['sin2thetaeeff']=SMcalculator.SIN2THETALEFF_DEFAULT
    aux_sm=SMcalculator.EWPOcalculator(scheme=scheme,input_dict=sminputs)

    sin2thetaleff_measurement_names = []
    sin2thetaleff_measurement_para = {}
    sin2thetaleff_measurement_prediction = {}
    for obs in measurement_names:
        if observable_mapping[obs] in sin2thetaleff_obs:
            sin2thetaleff_measurement_names.append(obs)
            sin2thetaleff_measurement_prediction[obs]=aux_sm.get(observable_mapping[obs])
            sin2thetaleff_measurement_para[obs]={}
            sin2thetaleff_measurement_para[obs][s2t_key]=aux_sm.derivative(observable_mapping[obs],'sin2thetaleff')

    res, cov, param_res, para_cov = LINAfit.lina_fit(sin2thetaleff_measurement_names, measurement_central, sin2thetaleff_measurement_prediction, measurement_covariance, [s2t_key], sin2thetaleff_measurement_para)
    sin2thetaleff=aux_sm.sin2thetaleff()+param_res[s2t_key]
    sin2thetaleff_err=para_cov[s2t_key][s2t_key]**0.5
    return sin2thetaleff,sin2thetaleff_err


def printModel(out, theo_err_treatment, para_err_treatment, sminputs, d6lin_type, d6quad_type, d8lin_type):
    """
    Prints the statistical model summary, including observed data, correlations, and error treatments.
    """
    column_width=20
    print('\n' + '*'*80)
    print('Statistical model:\n' + '='*80)

    # Header for the observables section
    header = 'Observable', 'Measurement', 'Prediction', 'Total Error'
    print(''.join([str(x).ljust(column_width, ' ') for x in header]))
    print('-'*80)

    # Print each observable's measurement, prediction, and error
    for o in out:
        print(o.ljust(column_width, ' ') + ''.join([str(x).ljust(column_width, ' ') for x in (
            out[o]['measurement'], out[o]['smprediction'], out[o]['error'])]))

    print('='*80)
    print('Correlation')
    print('-'*80)

    # Print the correlation matrix
    print(' '.join([x[:6].ljust(6, ' ') for x in [''] + list(out.keys())]))
    for oy in out:
        print(' '.join([oy[:6].ljust(
            6, ' ')] + ['{:6.3f}'.format(round(out[ox]['correlation'][oy], 3)) for ox in out]))

    print('='*80)

    # Theory uncertainties section
    print('Theory uncertainties:')
    print('-'*80)

    if theo_err_treatment == 'off':
        print('turned off')
    elif theo_err_treatment == 'covariance':
        print('included in total error and correlation printed above')
    elif theo_err_treatment == 'nuispar':
        print('modelled by the following nuisance parameters:')
        errs = []
        for o in out:
            if 'theoerr' in out[o]:
                for e in out[o]['theoerr']:
                    if e not in errs:
                        errs.append(e)
        for e in errs:
            print('{}, effect on {}'.format(e, ', '.join(
                [o + ': ' + str(out[o]['theoerr'][e]) for o in out if 'theoerr' in out[o] and e in out[o]['theoerr']])))

    print('='*80)

    # Input parameter dependence section
    print('Input parameter dependence:')
    print('-'*80)

    if para_err_treatment == 'off':
        print('turned off')
    elif para_err_treatment == 'covariance':
        print('included in error and correlation printed above')
    elif para_err_treatment == 'nuispar':
        print('included in parametrization:')
        for o in out:
            if 'smdependence' in out[o] and o not in sminputs:
                para = ''
                for x in out[o]['smdependence']:
                    a = out[o]['smdependence'][x]
                    if a > 0:
                        para += ' +'
                    else:
                        para += ' '
                    para += '{}*d{}'.format(a, x)
                print('{} = {}{}'.format(o, out[o]['smprediction'], para))

    print('='*80)

    # SMEFT parametrization section
    print('SMEFT parametrization:')
    print('-'*80)
    print('Linear dimension-six contributions: {}'.format(d6lin_type))
    print('Quadratic dimension-six contributions: {}'.format(d6quad_type))
    print('Dimension-eight contributions: {}'.format(d8lin_type))
    print('='*80)


def runSMfit(names, sminputs, meas, meas_cov, pred, smparas, smpara):
    """
    Runs the Standard Model (SM) fit and prints the results.
    """

    exp_err = {n: meas_cov[n][n]**0.5 for n in meas_cov}

    print('Performing SM fit...\n')

    # Fit the data with parametric errors
    fit_res, fit_cov, params_res, params_cov = LINAfit.lina_fit_with_nps(
        names, meas, meas_cov, pred, smparas, smpara,
        nps=[p for p in smparas if 'theoerr' in p]
    )

    # Perform indirect fit (excluding each observable one at a time)
    indirect = {}
    indirect_err = {}
    for n in names:
        fit_res_indirect, fit_cov_indirect, params_res_indirect, params_cov_indirect = LINAfit.lina_fit_with_nps(
            names, meas, meas_cov, pred, smparas, smpara,
            nps=[p for p in smparas if 'theoerr' in p],
            spectators=[n]
        )
        indirect[n] = fit_res_indirect[n]
        indirect_err[n] = fit_cov_indirect[n][n]**0.5

    width_const = 13

    # Print the fit results
    print('Analytic SM fit results:')
    header = ['Observable', 'Direct', '  +-',
              'Fit', '  +-', 'Indirect', '  +-', 'Pull']
    postFitTable(width_const,header,names,meas,exp_err,fit_res,fit_cov,indirect,indirect_err)

def postFitTable(width_const,header,names,meas,exp_err,fit_res,fit_cov,indirect,indirect_err):
    
    print('='*len(header)*width_const)
    print(''.join([' '+x.ljust(width_const-1, ' ') for x in header]))
    print('-'*len(header)*width_const)

    for n in names:
        exponent=5 if n=='Gmu' else 0
        if n in fit_res:
            rnd = -int(floor(log10(10**exponent*exp_err[n]))) + 1
            print(''.join([(n+(' [10^{}]'.format(-exponent) if exponent!=0 else '')).ljust(width_const, ' ')] +
                          [('{:'+str(width_const)+'.' + str(max(rnd, 0)) + 'f}').format(round(10**exponent*x, rnd)) if isinstance(x, float) else x.rjust(width_const, ' ')
                           for x in [meas[n], exp_err[n], fit_res[n], fit_cov[n][n]**0.5, indirect[n], indirect_err[n]]]),
                  '{:6.1f}'.format(round((fit_res[n] - meas[n]) / exp_err[n], 1)))

    print('='*len(header)*width_const)


def runSMEFTfit(names, meas, meas_cov, full_cov, pred, smdependence_para, theoerr_para, coefficients_d6, d6lin):
    """
    Runs the linear SMEFT fit for dimension-six operators and prints the results.
    """

    exp_err = {n: meas_cov[n][n]**0.5 for n in meas_cov}

    # parametrize deviations from SM expectations in terms of pulls (relative to exp error), for better numerics
    smdependence_para_pulls={}
    for o in smdependence_para:
        smdependence_para_pulls[o]={}
        for inp in smdependence_para[o]:
            smdependence_para_pulls[o][inp+'_pull']=smdependence_para[o][inp]*exp_err[inp]
            
    print('\nPerforming linear d6 EFT fits...\n')

    # One-at-a-time fit for each d6 coefficient
    print('One-at-a-time results:')
    for c in coefficients_d6:
        fit_res, fit_cov, params_res, params_cov, directions, blind_directions, unconstrained = LINAfit.eft_fit(
            [n for n in names if n in full_cov], meas, full_cov, pred,
            smdependence_para_pulls, theoerr_para, d6lin, [c]
        )
        err = params_cov[c][c]**0.5
        rnd = -int(floor(log10(err))) + 1
        print(c, round(params_res[c], rnd), '+-', round(err, rnd))

    # Multi-dimensional fit for all d6 coefficients
    fit_res, fit_cov, params_res, params_cov, directions, blind_directions, unconstrained = LINAfit.eft_fit(
        [n for n in names if n in full_cov], meas, full_cov, pred,
        smdependence_para_pulls, theoerr_para, d6lin, coefficients_d6, doevs=True
    )
    
    print('\nMulti-dimensional fit (central, err, pull):\n')
    for ic, c in enumerate(params_res):
        err = params_cov[c][c]**0.5
        rnd = -int(floor(log10(err))) + 1
        if directions is not None:
            print('Direction ', end='')
        print(f"{c}:", round(params_res[c], rnd), round(
            err, rnd), round(params_res[c]/err, 1), end='\t')
        if directions is None:
            print('')
        else:
            print('\n'+''.join(['{0:+}*{1}'.format(round(directions[ic][1][x], 3), x)
                               for x in directions[ic][1] if abs(directions[ic][1][x]) > 0.0005]))
            importance = {}
            for n in names:
                if n not in full_cov or n not in d6lin:
                    continue
                importance[n] = 0.
                for c in coefficients_d6:
                    importance[n] += d6lin[n][c] * \
                        directions[ic][1][c] * directions[ic][0]
                importance[n] /= full_cov[n][n]**0.5
            importance = dict(
                sorted(importance.items(), key=lambda item: abs(item[1]), reverse=True))
            print('Mainly constrained by: ' +
                  ', '.join([x for x in importance if abs(importance[x]) > 0.1**0.5]))
            print()
    # Handle blind directions if present
    if blind_directions is not None:
        for d in blind_directions:
            print('Blind direction: ' + ''.join(['{0:+}*{1}'.format(
                round(d[1][x], 2), x) for x in d[1] if abs(d[1][x]) > 0.1]))

    # Perform indirect fit for each observable
    indirect = {}
    indirect_err = {}
    for n in names:
        if n not in full_cov:
            continue
        fit_res_indirect, fit_cov_indirect, params_res_indirect, params_cov_indirect, directions, blind_directions, unconstrained = LINAfit.eft_fit(
            [n for n in names if n in full_cov], meas, full_cov, pred,
            smdependence_para_pulls, theoerr_para, d6lin, coefficients_d6,
            doevs=True, spectators=[n]
        )
        if n in unconstrained:
            indirect[n] = float('nan')
            indirect_err[n] = float('nan')
        else:
            indirect[n] = fit_res_indirect[n]
            indirect_err[n] = fit_cov_indirect[n][n]**0.5

    # Print the SMEFT fit results
    print('\nAnalytic SMEFT fit results:')
    width_const=13
    header = ['Observable', 'Direct', '  +-',
              'Fit', '  +-', 'Indirect', '  +-', 'Pull']
    postFitTable(width_const,header,names,meas,exp_err,fit_res,fit_cov,indirect,indirect_err)


if __name__ == '__main__':
    main()
