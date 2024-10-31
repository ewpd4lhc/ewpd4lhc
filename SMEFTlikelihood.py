import Utils
import SMcalculator,LINAfit
import logging

# Supported symmetries
SUPPORTED_SYMMETRIES = ['general','top', 'topU3l', 'U35']
LFU_SYMMETRIES = ['topU3l', 'U35']

class EWPDlikelihood:
    """
    Class creating the Electroweak Precision Data (EWPD) likelihood.

    This class combines SM predictions from the EWPO calculator class
    with SMEFT parametrizations to model the EWPD SMEFT likelihood.
    """

    def __init__(self, scheme, symmetry, observables=None, Lambda=1, d6linSource='SMEFTsim', d6quadSource='SMEFTsim', d8linSource=None, datapath='data/pole_observables_lfu.yml', parapath='data/pole_observables_lfu_SMEFTpara', dataupdates=None):
        """
        Initialize the EWPDlikelihood object.

        Args:
            scheme (str): The scheme for SM calculations and SMEFT parametrization ('alpha','MW','sin2theta','alphaMW','alphasin2theta').
            symmetry (str): Symmetry used in the model ('general','top','topU3l','U35').
            observables (list): List of observables to consider (None for all available).
            Lambda (float): New physics energy scale in TeV(!).
            d6linSource (str): Source of dimension-6 linear terms.
            d6quadSource (str): Source of dimension-6 quadratic terms.
            d8linSource (str): Source of dimension-8 linear terms.
            datapath (str): Path to the file containing measurement data.
            parapath (str): Path to the SMEFT parametrization folder.
            dataupdates (dict): Optional updates to measurement data. Default is None.
        """
        # Ensure the scheme is valid
        if not hasattr(SMcalculator.INPUTSCHEME, scheme):
            raise Exception(f"Scheme {scheme} not supported")

        # Load measurement data and relevant mappings
        self._all_measurement_names, self._measurement_central, self._measurement_covariance, self._observable_mapping, self._fixed_prediction = Utils.load_data(
            datapath)

        # Update measurements wrt data files if any updates provided
        if dataupdates is not None:
            self._update_measurements(dataupdates)

        # Initialize SM predictions based on central measurement values
        self._smpredictions = SMcalculator.EWPOcalculator(
            input_dict=self._measurement_central, scheme=getattr(SMcalculator.INPUTSCHEME, scheme))

        # Select observables for the likelihood calculation
        if observables is None:
            self._observables = [
                o for o in self._all_measurement_names if not o in self._smpredictions.inputs()]
        else:
            self._observables = []
            for o in observables:
                if not o in self._all_measurement_names:
                    raise Exception(
                        "No measurement for observable {}".format(o))
                if self._observable_mapping[o] not in self._smpredictions.inputs():
                    self._observables.append(o)

        # Set the energy scale and load parameter data for dimension-6 and dimension-8 terms
        self._Lambda = Lambda
        if d6linSource is not None:
            self._d6linpara = Utils.load_para(
                parapath + '/{}_{}_{}.yml'.format(d6linSource, symmetry, scheme), 'd6lin')
        else:
            self._d6linpara = None
        if d6quadSource is not None:
            self._d6quadpara = Utils.load_para(
                parapath + '/{}_{}_{}.yml'.format(d6quadSource, symmetry, scheme), 'd6quad')
        else:
            self._d6quadpara = None
        if d8linSource is not None:
            self._d8linpara = Utils.load_para(
                parapath + '/{}_{}_{}.yml'.format(d8linSource, symmetry, scheme), 'd8lin')
        else:
            self._d8linpara = None

        # Validate the selected observables
        for o in self._observables:
            if not o in self._all_measurement_names:
                raise Exception("Observable {} not found in data file. Available observables: {}".format(
                    o, ", ".join(self._all_measurement_names)))
            if self._fixed_prediction[o] is None and not self._observable_mapping[o] in self._smpredictions.outputs():
                raise Exception(
                    "Observable {} not implemented as a SM output".format(o))
            for para in self._d6linpara, self._d6quadpara, self._d8linpara:
                if para is not None and not self._observable_mapping[o] in para:
                    raise Exception(
                        "Observable {} not implemented for chosen parametrizations".format(o))

    def observables(self, incl_inputs=False):
        """
        Get the list of observables.

        Args:
            incl_inputs (bool): Whether to include inputs in the observable list.

        Returns:
            list: List of observable names.
        """
        if incl_inputs:
            return list(self._smpredictions.inputs() + self._observables)
        return list(self._observables)

    def observables_to_d8(self):
        """
        Get the list of observables for wich dimension-8 terms exist.

        Returns:
            list: List of observables for which dimension-8 operators are implemented.
        """
        res = []
        for o in self.observables():
            if o in self._d8linpara:
                res.append(o)
        return res

    def measured(self, normalized=False, incl_inputs=False):
        """
        Get the measured values for the observables.

        Args:
            normalized (bool): Whether to return normalized values.
            incl_inputs (bool): Whether to include inputs.

        Returns:
            dict: Measured values for observables.
        """
        res = {}
        for o in self.observables(incl_inputs):
            res[o] = float(self._measurement_central[o])
        return self._normalize(res) if normalized else res

    def predicted(self, incl_inputs=False, normalized=False):
        """
        Get the predicted values for the observables.

        Args:
            incl_inputs (bool): Whether to include inputs.
            normalized (bool): Whether to return normalized values.

        Returns:
            dict: Predicted values for observables.
        """
        res = {}
        for o in self.observables(incl_inputs):
            res[o] = self._smpredictions.get(
                self._observable_mapping[o]) if self._fixed_prediction[o] is None else self._fixed_prediction[o]
        return self._normalize(res) if normalized else res

    def covariance(self, incl_inputs=False, with_para_err=True, with_theo_err=True, with_d8_err=False, normalized=False):
        """
        Compute the total covariance matrix for all observables. Includes measurement and, by default, parametric and theory uncertainties.

        Args:
            incl_inputs (bool): Whether to include inputs in covariance.
            with_para_err (bool): Whether to include parametric uncertainties.
            with_theo_err (bool): Whether to include theory uncertainties.
            with_d8_err (bool): Whether to include uncertainty due to unknown dim-8 coefficients (assuming they are O(1))
            normalized (bool): Whether to return normalized values.

        Returns:
            dict: Covariance matrix.
        """
        cov = self._smpredictions.covariance(self._measurement_covariance, self._observable_mapping,
                                             add_exp_err=True, add_para_err=with_para_err, add_theo_err=with_theo_err)

        res = {}
        for o1 in self.observables(incl_inputs):
            res[o1] = {}
            for o2 in self.observables(incl_inputs):
                res[o1][o2] = cov[o1][o2]
        if with_d8_err:
            covd8 = self.covariance_d8error()
            for o1 in self.observables(incl_inputs):
                for o2 in self.observables(incl_inputs):
                    res[o1][o2] += covd8[o1][o2]
                    
        return self._normalize_cov(res) if normalized else res

    def correlation(self, incl_inputs=False, with_para_err=True, with_theo_err=True, with_d8_err=False):
        """
        Compute the correlation matrix for the observables.

        Args:
            incl_inputs (bool): Whether to include inputs in covariance.
            with_para_err (bool): Whether to include parametric uncertainties.
            with_theo_err (bool): Whether to include theory uncertainties.
            with_d8_err (bool): Whether to include uncertainty due to unknown dim-8 coefficients (assuming they are O(1))

        Returns:
            dict: Correlation matrix.
        """
        res = self.covariance(incl_inputs, with_para_err,
                              with_theo_err, with_d8_err, False)
        for o1 in self.observables(incl_inputs):
            for o2 in self.observables(incl_inputs):
                if o1 != o2:
                    res[o1][o2] /= (res[o1][o1] ** 0.5 * res[o2][o2] ** 0.5)
        for o in self.observables(incl_inputs):
            res[o][o] = 1.
        return res

    def error(self, incl_inputs=False, with_para_err=True, with_theo_err=True, with_d8_err=False, normalized=False):
        """
        Compute the error for the observables.

        Args:
            incl_inputs (bool): Whether to include inputs in covariance.
            with_para_err (bool): Whether to include parametric uncertainties.
            with_theo_err (bool): Whether to include theory uncertainties.
            with_d8_err (bool): Whether to include uncertainty due to unknown dim-8 coefficients (assuming they are O(1))

        Returns:
            dict: Errors for each observable.
        """
        cov = self.covariance(incl_inputs, with_para_err,
                              with_theo_err, with_d8_err, normalized)
        res = {}
        for o in self.observables(incl_inputs):
            res[o] = cov[o][o] ** 0.5
        return res

    def prediction_error(self, incl_inputs=False, with_para_err=True, with_theo_err=True, normalized=False):
        """
        Compute the prediction error for all observables, due to parametric and theory uncertainties.

        Args:
            incl_inputs (bool): Whether to include inputs.
            with_para_err (bool): Whether to include parameter uncertainties.
            with_theo_err (bool): Whether to include theoretical uncertainties.
            normalized (bool): Whether to return normalized values.

        Returns:
            dict: Prediction errors for each observable.
        """
        cov = self._smpredictions.covariance(self._measurement_covariance, self._observable_mapping,
                                             add_exp_err=False, add_para_err=with_para_err, add_theo_err=with_theo_err)
        res = {}
        for o in self.observables(incl_inputs):
            res[o] = cov[o][o] ** 0.5
        return self._normalize(res) if normalized else res

    def covariance_d8error(self, normalized=False):
        """
        Compute the covariance matrix for dimension-8 operators.
        Each dim-8 coefficient introduces an O(1) uncertainty.

        Args:
            normalized (bool): Whether to return normalized values.

        Returns:
            dict: Covariance matrix for dimension-8 operators.
        """
        pred = self.predicted()
        res = {}
        for o1 in self.observables():
            res[o1] = {}
            for o2 in self.observables():
                res[o1][o2] = 0.
        para = self.d8lin()
        for c in self.wilson_coefficients_d8():
            for o1 in self.observables():
                for o2 in self.observables():
                    if c in para[o2] and c in para[o1]:
                        res[o1][o2] += (para[o1][c] * para[o2][c])
        return self._normalize_cov(res) if normalized else res

    def d6lin(self, normalized=False):
        """
        Get the dimension-6 linear contribution for each observable.

        Args:
            normalized (bool): Whether to return normalized values.

        Returns:
            dict: Dimension-6 linear contributions.
        """
        if self._d6linpara is None:
            return None
        res = {}
        for o in self.observables():
            res[o] = dict(self._d6linpara[self._observable_mapping[o]])
            for c in res[o]:
                res[o][c] /= self._Lambda**2
        return self._normalize(res) if normalized else res

    def d6quad(self, normalized=False):
        """
        Get the dimension-6 quadratic contribution for each observable.

        Args:
            normalized (bool): Whether to return normalized values.

        Returns:
            dict: Dimension-6 quadratic contributions.
        """
        if self._d6quadpara is None:
            return None
        res = {}
        for o in self.observables():
            if self._d6quadpara is not None:
                res[o] = dict(self._d6quadpara[self._observable_mapping[o]])
                for c in res[o]:
                    res[o][c] /= self._Lambda**4
        return self._normalize(res) if normalized else res

    def d8lin(self, normalized=False):
        """
        Get the dimension-8 linear contribution for each observable.

        Args:
            normalized (bool): Whether to return normalized values.

        Returns:
            dict: Dimension-8 linear contributions.
        """
        if self._d8linpara is None:
            return None
        res = {}
        for o in self.observables():
            res[o] = dict(self._d8linpara[self._observable_mapping[o]])
            for c in res[o]:
                res[o][c] /= self._Lambda**4
        return self._normalize(res) if normalized else res

    def linSMpara(self, incl_inputs=False, with_theo_err=True, normalized=False):
        """
        Get the linearized SM parameterization.

        Args:
            incl_inputs (bool): Whether to include SM input observables (with trivial dependence).
            with_theo_err (bool): Whether to include the dependence on theoretical uncertainties.
            normalized (bool): Whether to return normalized values.

        Returns:
            tuple: Inputs and linear parameters for each observable.
        """
        res = {}
        inputs, para = self._smpredictions.linear_para(
            incl_nuis_pars=with_theo_err)
        for o in self.observables(incl_inputs):
            if self._observable_mapping[o] in para:
                res[o] = para[self._observable_mapping[o]]
            elif o in self._fixed_prediction:
                res[o] = {}
            else:
                raise Exception(
                    f"Unknown SM dependence for {o} (if it's a fixed prediction add it to data file)")

        return list(inputs), res

    def wilson_coefficients_d6(self):
        """
        Get the list of dimension-6 Wilson coefficients with non-zero impact.

        Returns:
            list: Wilson coefficients for dimension-6 operators.
        """
        res = []
        for o in self.observables():
            if self._d6linpara is not None:
                for c in self._d6linpara[self._observable_mapping[o]]:
                    if not c in res and self._d6linpara[self._observable_mapping[o]][c]!=0:
                        res.append(c)
            if self._d6quadpara is not None:
                if self._observable_mapping[o] in self._d6quadpara:
                    for c2 in self._d6quadpara[self._observable_mapping[o]]:
                        if self._d6quadpara[self._observable_mapping[o]][c2]==0:
                            continue
                        for c in c2.split('*'):
                            if not c in res:
                                res.append(c)
        return res

    def wilson_coefficients_d8(self):
        """
        Get the list of dimension-8 Wilson coefficients with non-zero impact.

        Returns:
            list: Wilson coefficients for dimension-8 operators.
        """
        res = []
        if self._d8linpara is not None:
            for o in self.observables():
                for c in self._d8linpara[self._observable_mapping[o]]:
                    if not c in res and self._d8linpara[self._observable_mapping[o]][c]!=0:
                        res.append(c)
        return res
    
    def _update_measurements(self,dataupdates):
        for o in dataupdates:
            if not 'central' in dataupdates[o]:
                raise Exception(f"Need to give central value for data update {o}")
            if not 'error' in dataupdates[o]:
                raise Exception(f"Need to give error for data update {o}")
            if not o in self._all_measurement_names:
                self._all_measurement_names.append(o)
                self._measurement_covariance[o] = {}
                if not 'prediction' in dataupdates[o]:
                    raise Exception(f"Need to give prediction for data update {o}")
            self._measurement_central[o] = dataupdates[o]['central']
            self._measurement_covariance[o][o] = dataupdates[o]['error']**2
            if 'prediction' in dataupdates[o]:
                if isinstance(dataupdates[o],(float,int)):
                    self._fixed_prediction[o] = dataupdates[o]['prediction']
                    self._observable_mapping[o] = o
                else:
                    self._fixed_prediction[o] = None
                    self._observable_mapping[o] = dataupdates[o]['prediction']
            # Update covariance for related observables assuming no correlation
            for o1 in self._all_measurement_names:
                if not o == o1:
                    self._measurement_covariance[o][o1] = 0.
                    self._measurement_covariance[o1][o] = 0.

    def _normalize(self, d):
        """
        Normalize a dictionary (or dict of dicts) of values by the predicted values.

        Args:
            d (dict): Dictionary of observable values.

        Returns:
            dict: Normalized dictionary of values.
        """
        res = {}
        pred = self.predicted(incl_inputs=True)
        for o in d:
            if isinstance(d[o], dict):
                res[o] = {}
                for c in d[o]:
                    res[o][c] = d[o][c] / pred[o]
            else:
                res[o] = d[o] / pred[o]
        return res

    def _normalize_cov(self, d):
        """
        Normalize a covariance matrix by predicted values.

        Args:
            d (dict): Covariance matrix.

        Returns:
            dict: Normalized covariance matrix.
        """
        res = {}
        pred = self.predicted(incl_inputs=True)
        for o1 in d:
            res[o1] = {}
            for o2 in d[o1]:
                res[o1][o2] = d[o1][o2] / (pred[o1] * pred[o2])
        return res

