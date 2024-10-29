import numpy as np
import Utils


def lina_fit(names, measured, predicted, covariance, parameters, parametrization, spectators=[], doevs=False, evcutoff=1e5, unconstrainedcutoff=1e-2):
    """
    Solves a chi2-fit of a linear parametrization with linear algebra.

    Args:
    names (list): List of observable names.
    measured (dict): Measured values of the observables.
    predicted (dict): Predicted values of the observables.
    covariance (dict): Covariance matrix of the measured values.
    parameters (list): List of model parameters.
    parametrization (dict): Dictionary describing how the observables depend on the parameters.
    spectators (list): List of observables to be treated as spectators (ignored in fit). Default is [].
    doevs (bool): Whether to perform an eigenvalue (EV) analysis. Default is False. 
    evcutoff (float): Cutoff for eigenvector with lowest eigenvalue considered. Default is 1e5.
    unconstrainedcutoff (float): Cutoff for Observables considered unconstrained as they align with blind direction. Default is 1e-2.

    Returns:
    Tuple containing:
        - dres (dict): Post-fit values for each observable.
        - dcov (dict): Post-fit covariance matrix of observables.
        - dparam_res (dict): Fitted parameter values.
        - dpara_cov (dict): Covariance matrix of the fitted parameters.
        - sensitive_directions (list, optional): Sensitive directions from EV analysis (if doevs is True).
        - blind_directions (list, optional): Blind directions from EV analysis (if doevs is True).
        - unconstrained (list, optional): List of observables unconstrained by the fit (if doevs is True).
    """

    # to be safe, make new lists
    names = list(names)
    params = list(parameters)

    # Convert measured and predicted values to arrays
    meas = np.array([measured[n] for n in names])
    pred = np.array([predicted[n] for n in names])
    # Differeces between measured and predicted
    deltas = np.array([measured[n] - predicted[n] for n in names])

    # Convert covariance and parametrization dicts to numpy  matrices
    cov = np.array(Utils.toMatrix(names, names, covariance))
    para = np.array(Utils.toMatrix(names, params, parametrization))

    # Invert the covariance matrix
    covI = np.linalg.inv(cov)

    # Zero out covariance for spectators (observables to be ignored)
    for i, n in enumerate(names):
        for j, m in enumerate(names):
            if n in spectators or m in spectators:
                covI[i][j] = 0

    repara = None

    # Eigenvalue (EV) analysis to identify sensitive and blind directions
    if doevs:
        # Inverse covariance (Fisher information) in space of parametrization
        a = np.transpose(para) @ covI @ para

        # Its eigenvalues
        w, v = np.linalg.eigh(a)
        evs = list(zip(w, np.transpose(v)))

        # Sort by eigenvalue, most sensitive direction first
        evs.sort(reverse=True, key=lambda x: x[0])

        # Convert eigenvalue sensitivities and directions into readable format
        evs = [(1 / abs(r[0]) ** 0.5, dict((x, y) for y, x in sorted(zip(r[1].tolist(),
                params), key=lambda x: abs(x[0]), reverse=True))) for r in evs]
        sensitive_directions = []
        blind_directions = []

        # Sort into blind and sensitive directions
        for e, v in evs:
            if e > evcutoff * evs[0][0]:
                blind_directions.append((e, v))
            else:
                sensitive_directions.append((e, v))

        # In case there are blind directions, create a reparametrization matrix from parameter basis to EV basis
        if len(blind_directions) > 0:
            reparas = []
            for d in sensitive_directions:
                reparas.append(Utils.toList(params, d[1]))
            repara = np.array(reparas)

        # Identify unconstrained observables
        # They are unconstrained if their parametrization receives contributions from blind directions
        unconstrained = []
        for n in names:
            # Fit direction
            p = np.array(Utils.toList(params, parametrization[n]))
            pnorm = np.linalg.norm(p)
            # I forgot why I did this :(
            if pnorm < 1e-30:
                continue
            # Determine which observables get contributions from blind directions
            # This is the case if normalized_param dot normalized_blind_dir >
            maxproj = 0.0
            for direction in blind_directions:
                d = np.array(Utils.toList(params, direction[1]))
                dnorm = np.linalg.norm(d)
                maxproj = max(maxproj, abs(np.dot(p, d) / pnorm / dnorm))
            if maxproj > unconstrainedcutoff:
                unconstrained.append(n)
    # Reparametrize if sensitive directions have been identified
    # In that case eigenvectors are the parameters
    if repara is not None:
        para = para @ np.transpose(repara)
        params = range(len(sensitive_directions))

    # Solve for fitted parameters and compute covariance
    paraT = np.transpose(para)
    param_res = np.linalg.solve(paraT @ covI @ para, paraT @ covI @ deltas)
    covP = np.linalg.inv(paraT @ covI @ para)

    # Convert results to dictionaries
    dparam_res = {n: float(r) for n, r in zip(params, param_res)}
    dpara_cov = {n: {m: float(covP[i][j]) for j, m in enumerate(
        params)} for i, n in enumerate(params)}

    # Post-fit residuals for observables
    res = deltas - para @ param_res
    # And the covariance (linear error propagation)
    cov_postfit = para @ covP @ paraT

    # Post fit result (subtract residual)
    dres = {n: float(m - r) for n, r, m in zip(names, res, meas)}
    # Post fit cov (create a dict from the matrix)
    dcov = {n: {m: float(cov_postfit[i][j]) for j, m in enumerate(
        names)} for i, n in enumerate(names)}

    # If a reparametrization has been performed, information on sensitive directions is reported
    if repara is not None:
        return dres, dcov, dparam_res, dpara_cov, sensitive_directions, blind_directions, unconstrained
    else:
        return dres, dcov, dparam_res, dpara_cov


def lina_fit_with_nps(measurement_names, measurement_central, measurement_covariance, predictions, parameters, parametrization, nps, spectators=[], doevs=False):
    """
    A lina fit with nuisance parameters (nps). The nuisance parameter dependence should already be part of the parametrization.
    Gaussian constraints are added for each nuisance parameter.


    Args:
    measurement_names (list): Names of the measured observables.
    measurement_central (dict): Central values of the measurements.
    measurement_covariance (dict): Covariance matrix of the measurements.
    predictions (dict): Predicted values of the observables.
    parameters (list): List of model parameters.
    parametrization (dict): Parametrization of the observables in terms of the parameters.
    nps (list): List of nuisance parameters (NPs) to account for theoretical errors.
    spectators (list): List of observables to be treated as spectators (ignored in fit). Default is [].
    doevs (bool): Whether to perform an eigenvalue (EV) analysis. Default is False.

    Returns:
    Result of lina_fit.
    """

    # Initialize central values of measurement and NPs
    central_values = dict(measurement_central)

    # Add nuisance parameters to central values and predictions (add them to our Gaussian model as independent measurement)
    for p in nps:
        central_values[p] = 0.0
        predictions[p] = 0.0

    # Extend measurement names to include nuisance parameters
    measurment_plus_np_names = measurement_names + nps

    # Update parametrization to include nuisance parameters
    actual_para = dict(parametrization)
    for p in nps:
        actual_para[p] = {}
        actual_para[p][p] = 1

    # Extend covariance matrix to include nuisance parameters
    full_covariance = {}
    for n in measurement_names:
        full_covariance[n] = {}
        for m in measurement_names:
            full_covariance[n][m] = measurement_covariance[n][m]
    for p in nps:
        full_covariance[p] = {}
        full_covariance[p][p] = 1.0

    # Call lina_fit to perform the full fit
    return lina_fit(measurment_plus_np_names, central_values, predictions, full_covariance, parameters, actual_para, spectators, doevs=doevs)


def eft_fit(names, central_values, covariance, predictions, sm_para, theoerr_para, eft_para, eft_paras, spectators=[], doevs=False):
    """
    Fits the data using a linear EFT parametrization.

    Args:
    names (list): Names of the observables.
    central_values (dict): Central values of the measured observables.
    covariance (dict): Covariance matrix of the observables.
    predictions (dict): Predicted values of the observables.
    sm_para (dict): Parametrization of the observables in terms of Standard Model (SM) parameters.
    theoerr_para (dict): Theoretical error parametrization.
    eft_para (dict): Parametrization of the observables in terms of EFT parameters.
    eft_paras (list): List of EFT parameters (dimension-six coefficients).
    spectators (list): List of observables to be treated as spectators (ignored in fit). Default is [].
    doevs (bool): Whether to perform an eigenvalue (EV) analysis. Default is False.

    Returns:
    Result of lina_fit (but with fixed number of return values).
    """

    parameters = []
    nps = []
    para = {}

    # Construct parametrization for both SM and EFT parameters
    for n in names:
        para[n] = {}
        for p in sm_para[n]:
            if p not in parameters:
                parameters.append(p)
            para[n][p] = sm_para[n][p]
        if n in eft_para:
            for p in eft_paras:
                if p not in parameters:
                    parameters.append(p)
                para[n][p] = eft_para[n][p]
        for p in theoerr_para[n]:
            if p not in parameters:
                parameters.append(p)
                nps.append(p)
            para[n][p] = theoerr_para[n][p]

    # Perform the fit using nuisance parameters (NPs)
    res = lina_fit_with_nps(names, central_values, covariance,
                            predictions, parameters, para, nps, spectators, doevs)

    # If EV analysis is performed, unpack the results accordingly
    if len(res) > 4:
        fit_res, fit_cov, params_res, params_cov, directions, blind_directions, unconstrained = res
    else:
        fit_res, fit_cov, params_res, params_cov = res
        directions = None
        blind_directions = None
        unconstrained = []

    return fit_res, fit_cov, params_res, params_cov, directions, blind_directions, unconstrained
