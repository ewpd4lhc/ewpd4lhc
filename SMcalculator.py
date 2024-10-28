from math import log
import Utils
import numpy as np
import logging

# SM inputs
SMINPUTS_ALPHASCHEME = ['MH', 'mt', 'alphas', 'MZ', 'Deltaalpha']
SMINPUTS_MWSCHEME = ['MH', 'mt', 'alphas', 'MZ', 'MW']
# The EWPOs predicted or used as input by this class
EWPOS = ['MH', 'mt', 'alphas', 'MZ', 'Deltaalpha', 'GammaZ', 'Rl', 'Re', 'Rmu', 'Rtau', 'Rc', 'Rb', 'sigmahad', 'Al', 'Ae', 'Amu', 'Atau',
         'AFBl', 'AFBe', 'AFBmu', 'AFBtau', 'AFBb', 'AFBc', 'Ab', 'Ac', 'MW', 'GammaW', 'BrWhad', 'BrWe', 'BrWmu', 'BrWtau', 'BrWl', 'sin2thetaleff']


class EWPOcalculator:
    """
    A class that calculates Standard Model predictions for electroweak precision observables (EWPOs)
    The class supports two schemes ('MW' and 'alpha') for parameter inputs .
    It provides functionality to compute various EWPOs, their derivatives, and theory and parametric covariance.
    """

    def __init__(self, MH=None, mt=None, alphas=None, MZ=None, Deltaalpha=None, MW=None, input_dict=None, scheme=None):
        """
        Initializes the EWPOcalculator with either a set of 5 input parameters (including one of MW or Deltaalpha)
        or a dictionary of observable values (can be more than 5 but most include at least 5 inputs) + the scheme chosen.
        Units are GeV.

        Example usage:
        sm = EWPOcalculator(MH=125.,mt=173.,alphas=0.118,MZ=91.2,MW=80.4)
        or
        sm = EWPOcalculator(input_dict=dict_of_input_par_values,scheme='MW')

        Args:
            MH (float): Higgs mass.
            mt (float): Top quark mass.
            alphas (float): Strong coupling constant.
            MZ (float): Z boson mass.
            Deltaalpha (float): Electroweak alpha parameter (optional).
            MW (float): W boson mass (optional).
            input_dict (dict): Dictionary containing SM input parameters.
            scheme (str): Specifies which scheme to use ('MW' or 'alpha').

        Raises:
            Exception: If the input values are inconsistent or missing required parameters.
        """
        if not input_dict is None:
            # Ensure input dictionary is used exclusively when provided
            for x in MH, mt, alphas, MZ, Deltaalpha, MW:
                if not x is None:
                    raise Exception(
                        "Set input with dict only or without dict please")
            if scheme is None or scheme not in ['MW', 'alpha']:
                raise Exception("Need to choose 'MW' or 'alpha' scheme")
            # Set values based on the MW or alpha scheme
            if scheme == 'MW':
                for x in SMINPUTS_MWSCHEME:
                    if x not in input_dict:
                        raise Exception("Need to set {}".format(x))
                MH = input_dict['MH']
                mt = input_dict['mt']
                alphas = input_dict['alphas']
                MZ = input_dict['MZ']
                MW = input_dict['MW']
            elif scheme == 'alpha':
                for x in SMINPUTS_ALPHASCHEME:
                    if x not in input_dict:
                        raise Exception("Need to set {}".format(x))
                MH = input_dict['MH']
                mt = input_dict['mt']
                alphas = input_dict['alphas']
                MZ = input_dict['MZ']
                Deltaalpha = input_dict['Deltaalpha']
        else:
            # Ensure all individual parameters are set if not using a dictionary
            for x in MH, mt, alphas, MZ:
                if x is None:
                    raise Exception(
                        "Need to give all 5 inputs (e.g.: EWPOcalculator(MH=125,mt=173,MW=80.4,MZ=91.2,alphas=0.118)) or dict with all 5 inputs")
        if MW is None and Deltaalpha is None:
            raise Exception("Need to set MW or Deltaalpha")

        # Set interpolation constants and initial values
        self._set_interpol_constants()

        self._MH_original = MH
        self._mt_original = mt
        self._alphas_original = alphas
        self._MZ_original = MZ
        self._Deltaalpha_original = Deltaalpha
        self._MW_original = MW

        # Determine whether MW or alpha scheme is used
        self._MWscheme = False
        if MW is not None:
            if Deltaalpha is not None:
                raise Exception("Only MW OR Deltaalpha supported as input")
            self._MWscheme = True
            self._SMinputs = SMINPUTS_MWSCHEME
            if scheme is not None and scheme == 'alpha':
                raise Exception(
                    "Given inputs for MW scheme but specified 'alpha'")
        else:
            self._SMinputs = SMINPUTS_ALPHASCHEME
            if scheme is not None and scheme == 'MW':
                raise Exception(
                    "Given inputs for alpha scheme but specified 'MW'")

        self._outputs = [x for x in EWPOS if x not in self._SMinputs]

        # Initialize with the provided inputs
        self.update(MH, mt, alphas, MZ, Deltaalpha, MW)

    def reset(self):
        """
        Resets the EWPOcalculator to its original input values.
        """
        self.update(MH=self._MH_original, mt=self._mt_original, alphas=self._alphas_original,
                    MZ=self._MZ_original, Deltaalpha=self._Deltaalpha_original, MW=self._MW_original)

    def update(self, MH=None, mt=None, alphas=None, MZ=None, Deltaalpha=None, MW=None, input_dict=None):
        """
        Updates the SM parameters and interpolates new values for observables based on input changes.

        Args:
            MH (float): Updated Higgs mass.
            mt (float): Updated Top quark mass.
            alphas (float): Updated strong coupling constant.
            MZ (float): Updated Z boson mass.
            Deltaalpha (float): Updated electroweak alpha parameter (optional).
            MW (float): Updated W boson mass (optional).
            input_dict (dict): Dictionary of updated SM input parameters.

        Raises:
            Exception: If the inputs are inconsistent or both MW and Deltaalpha are set.
        """
        if not input_dict is None:
            # Ensure input dictionary is used exclusively when provided
            for x in MH, mt, alphas, MZ, Deltaalpha, MW:
                if not x is None:
                    raise Exception(
                        "Set input with dict only or without dict please")
            # Update the values from the dictionary
            MH = input_dict.get('MH')
            mt = input_dict.get('mt')
            alphas = input_dict.get('alphas')
            MZ = input_dict.get('MZ')
            if self._MWscheme:
                MW = input_dict.get('MW')
            else:
                Deltaalpha = input_dict.get('Deltaalpha')

        if MW is not None and Deltaalpha is not None:
            raise Exception("Should only set MW OR Deltaalpha")

        # Update the specific parameters based on the provided inputs
        if MH is not None:
            self._update_MH(MH)
        if mt is not None:
            self._update_mt(mt)
        if alphas is not None:
            self._update_alphas(alphas)
        if MZ is not None:
            self._update_MZ(MZ)
        if Deltaalpha is not None:
            self._update_Deltaalpha(Deltaalpha)
        if MW is not None:
            self._update_MW(MW)

        if self._MWscheme:
            self._update_Deltaalpha_common(self._Deltaalpha_from_MW())

    def MW(self):
        """
        Returns the value of the W boson mass.

        Returns:
            float: MW value if MW scheme is used, otherwise MW derived from Deltaalpha.
        """
        if self._MWscheme:
            return self._MW
        else:
            return self._MW_from_Deltaalpha()

    def Deltaalpha(self):
        """
        Returns the value of the running alpha parameter Deltaalpha.

        Returns:
            float: Deltaalpha value.
        """
        if self._MWscheme:
            return self._Deltaalpha_from_MW()
        return self._Deltaalpha

    def MH(self):
        """
        Returns the value of the Higgs mass.

        Returns:
            float: Higgs mass (MH).
        """
        return self._MH

    def mt(self):
        """
        Returns the value of the top quark mass.

        Returns:
            float: Top quark mass (mt).
        """
        return self._mt

    def alphas(self):
        """
        Returns the value of the strong coupling constant (alphas).

        Returns:
            float: Strong coupling constant (alphas).
        """
        return self._alphas

    def MZ(self):
        """
        Returns the value of the Z boson mass.

        Returns:
            float: Z boson mass (MZ).
        """
        return self._MZ

    def GammaZ(self):
        """
        Returns the Z boson decay width (GammaZ).

        Returns:
            float: Z boson decay width in GeV.
        """
        return self._GammaObs(*self._A_GAMMA_Z) / 1e3

    def GammaW(self):
        """
        Returns the W boson decay width (GammaW).

        Returns:
            float: W boson decay width in GeV.
        """
        return self._GammaW()

    def BrWhad(self):
        """
        Returns the W boson hadronic branching ratio (BrWhad).

        Returns:
            float: W boson hadronic branching ratio.
        """
        return self._BrWhad()

    def BrWl(self):
        """
        Returns the W boson branching ratio to a specific lepton flavour.

        Returns:
            float: W boson branching ratio to leptons  (Br(W -> l)).
        """
        return self._BrWl()

    def BrWe(self):
        """
        Returns the W boson branching ratio to electrons.

        Returns:
            float: W boson branching ratio to electrons (Br(W -> e)).
        """
        return self._BrWl()

    def BrWmu(self):
        """
        Returns the W boson branching ratio to muons.

        Returns:
            float: W boson branching ratio to muons (Br(W -> mu)).
        """
        return self._BrWl()

    def BrWtau(self):
        """
        Returns the W boson branching ratio to taus.

        Returns:
            float: W boson branching ratio to taus (Br(W -> tau)).
        """
        return self._BrWl()

    def Rl(self):
        """
        Returns the ratio of Z boson decays to leptons over hadrons (Rl).

        Returns:
            float: Ratio of leptonic to hadronic Z decays (Rl).
        """
        return self._GammaObs(*self._A_R_L) / 1e3

    def Rmu(self):
        """
        Returns the ratio of Z boson decays to muons over hadrons (Rmu).

        Returns:
            float: Ratio of Z decays to muons (same as Rl).
        """
        return self.Rl()

    def Re(self):
        """
        Returns the ratio of Z boson decays to electrons over hadrons (Re).

        Returns:
            float: Ratio of Z decays to electrons (same as Rl).
        """
        return self.Rl()

    def Rtau(self):
        """
        Returns the ratio of Z boson decays to taus over hadrons (Rtau).

        Returns:
            float: Ratio of Z decays to taus. Slightly different from Rl due to mass difference.
        """
        return self.Rl() * 1.0023

    def Rc(self):
        """
        Returns the ratio of Z boson decays to charm quarks over hadrons (Rc).

        Returns:
            float: Ratio of Z decays to charm quarks.
        """
        return self._GammaObs(*self._A_R_C) / 1e5

    def Rb(self):
        """
        Returns the ratio of Z boson decays to bottom quarks over hadrons (Rb).

        Returns:
            float: Ratio of Z decays to bottom quarks.
        """
        return self._GammaObs(*self._A_R_B) / 1e5

    def sigmahad(self):
        """
        Returns the hadronic cross-section at the Z pole (sigma_had).

        Returns:
            float: Hadronic cross-section at the Z boson pole.
        """
        return self._GammaObs(*self._A_SIGMA_HAD)

    def Al(self):
        """
        Returns the leptonic asymmetry parameter (Al).

        Returns:
            float: Leptonic asymmetry parameter Al.
        """
        return self._Al()

    def Ae(self):
        """
        Returns the electron asymmetry parameter (Ae).

        Returns:
            float: Electron asymmetry parameter (same as Al).
        """
        return self.Al()

    def Amu(self):
        """
        Returns the muon asymmetry parameter (Amu).

        Returns:
            float: Muon asymmetry parameter (same as Al).
        """
        return self.Al()

    def Atau(self):
        """
        Returns the tau asymmetry parameter (Atau).

        Returns:
            float: Tau asymmetry parameter (same as Al).
        """
        return self.Al()

    def Ab(self):
        """
        Returns the bottom quark asymmetry parameter (Ab).

        Returns:
            float: Bottom quark asymmetry parameter.
        """
        return self._Ab()

    def Ac(self):
        """
        Returns the charm quark asymmetry parameter (Ac).

        Returns:
            float: Charm quark asymmetry parameter.
        """
        return self._Ac()

    def AFBl(self):
        """
        Returns the forward-backward asymmetry for leptons (AFBl).

        Returns:
            float: Forward-backward asymmetry for leptons.
        """
        return self._AFBl()

    def AFBe(self):
        """
        Returns the forward-backward asymmetry for electrons (AFBe).

        Returns:
            float: Forward-backward asymmetry for electrons (same as AFBl).
        """
        return self.AFBl()

    def AFBmu(self):
        """
        Returns the forward-backward asymmetry for muons (AFBmu).

        Returns:
            float: Forward-backward asymmetry for muons (same as AFBl).
        """
        return self.AFBl()

    def AFBtau(self):
        """
        Returns the forward-backward asymmetry for taus (AFBtau).

        Returns:
            float: Forward-backward asymmetry for taus (same as AFBl).
        """
        return self.AFBl()

    def AFBb(self):
        """
        Returns the forward-backward asymmetry for bottom quarks (AFBb).

        Returns:
            float: Forward-backward asymmetry for bottom quarks.
        """
        return self._AFBb()

    def AFBc(self):
        """
        Returns the forward-backward asymmetry for charm quarks (AFBc).

        Returns:
            float: Forward-backward asymmetry for charm quarks.
        """
        return self._AFBc()

    def sin2thetaleff(self):
        """
        Returns the effective weak mixing angle (sin^2(theta)_eff).

        Returns:
            float: Effective weak mixing angle sin^2(theta)_eff.
        """
        return self._sin2thetaleff()

    def inputs(self):
        """
        Returns the list of Standard Model (SM) input parameters used in the calculation.

        Returns:
            list: A list of SM input parameters.
        """
        return list(self._SMinputs)

    def outputs(self):
        """
        Returns the list of output observables computed by the SM calculator.

        Returns:
            list: A list of observables.
        """
        return list(self._outputs)

    def get(self, obs):
        """
        Retrieves the value of a specified observable.

        Args:
            obs (str): The observable to retrieve.

        Returns:
            float: The value of the specified observable.

        Raises:
            Exception: If the observable is not implemented.
        """
        if not obs in EWPOS:
            raise Exception(str(obs) + " is not implemented")
        else:
            return getattr(self, self._lfu_version(obs))()

    def getall(self):
        """
        Retrieves the values for all observables in the SM calculator.

        Returns:
            dict: A dictionary where keys are observable names and values are their corresponding values.
        """
        d = {}
        for obs in EWPOS:
            d[obs] = self.get(obs)
        return dict(d)

    def set(self, name, value):
        """
        Sets the value of a specific Standard Model input and updates the state accordingly.

        Args:
            name (str): The name of the SM input parameter.
            value (float): The value to set for the input parameter.

        Raises:
            Exception: If the input parameter name is not recognized.
        """
        if not name in self._SMinputs:
            raise Exception(str(name) + " is not a SM input")
        else:
            getattr(self, '_update_' + name)(value)

        if self._MWscheme:
            self._update_Deltaalpha_common(self._Deltaalpha_from_MW())

    def derivative(self, obs, inp):
        """
        Computes the derivative of an observable with respect to an input parameter.

        Args:
            obs (str): The observable whose derivative is to be computed.
            inp (str): The input parameter with respect to which the derivative is taken.

        Returns:
            float: The derivative value.

        Raises:
            Exception: If the observable or input parameter is not recognized or implemented.
        """
        obs = self._lfu_version(obs)

        if obs == inp:
            return 1.0

        if not obs in EWPOS:
            raise Exception(str(obs) + " is not implemented")

        if obs in self._SMinputs:
            raise Exception(str(obs) + " is an input")

        if not inp in self._SMinputs:
            raise Exception(str(inp) + " is not an input")

        if not self._MWscheme:
            if hasattr(self, '_d{}_d{}'.format(obs, inp)):
                return getattr(self, '_d{}_d{}'.format(obs, inp))()
            else:
                logging.debug('No {} derivative of {}'.format(inp, obs))
                return 0
        else:
            # Special handling for MW scheme and Deltaalpha
            dalpha_dmw = self._dDeltaalpha_dMW()
            if obs == 'Deltaalpha':
                dobs_dalpha = 1.0
            elif hasattr(self, '_d{}_d{}'.format(obs, 'Deltaalpha')):
                dobs_dalpha = getattr(
                    self, '_d{}_d{}'.format(obs, 'Deltaalpha'))()
            else:
                logging.debug('No {} derivative of {}'.format(inp, obs))
                dobs_dalpha = 0.0

            if inp == 'MW':
                return dobs_dalpha * dalpha_dmw
            else:
                dobs_dinp = getattr(self, '_d{}_d{}'.format(obs, inp))(
                ) if hasattr(self, '_d{}_d{}'.format(obs, inp)) else 0.0
                if obs == 'Deltaalpha':
                    return dobs_dinp
                else:
                    dmw_dinp = getattr(self, '_d{}_d{}'.format('MW', inp))()
                    return dobs_dinp - dobs_dalpha * dalpha_dmw * dmw_dinp

    def derivatives(self, obs):
        """
        Computes the derivatives of an observable with respect to all input parameters.

        Args:
            obs (str): The observable whose derivatives are to be computed.

        Returns:
            dict: A dictionary where the keys are input parameters and values are the derivative values.
        """
        obs = self._lfu_version(obs)
        derivatives_dict = {}

        # Calculate derivative for each SM input
        for inp in self._SMinputs:
            derivatives_dict[inp] = self.derivative(obs, inp)

        return derivatives_dict

    def linear_para(self, incl_nuis_pars=False, outputs=None):
        """
        Computes the linear parameterization of the model's observables with respect to inputs.

        Args:
            incl_nuis_pars (bool): Flag to include nuisance parameters in the parameterization.
            outputs (list, optional): List of outputs to include. Defaults to all outputs.

        Returns:
            tuple: A list of parameters and a dictionary of linear coefficients for each observable.
        """
        if outputs is None:
            outputs = self.outputs()

        parameterization = {}
        parameters = self.inputs()

        for o in EWPOS:
            parameterization[o] = {}
            if o in self.inputs():
                parameterization[o][o] = 1.0
                continue
            # Calculate derivatives for each observable and input
            for i in self.inputs():
                parameterization[o][i] = self.derivative(o, i)
            if incl_nuis_pars:
                if self.uncorr_theoerr(o) != 0:
                    errname = 'theoerr_{}'.format(o)
                    parameterization[o][errname] = self.uncorr_theoerr(o)
                    if errname not in parameters:
                        parameters.append(errname)
                if self.sin2theta_theoerr(o) != 0:
                    errname = 'theoerr_sin2theta'
                    parameterization[o][errname] = self.sin2theta_theoerr(o)
                    if errname not in parameters:
                        parameters.append(errname)

        return parameters, parameterization

    def uncorr_theoerr(self, obs):
        """
        Retrieves the uncorrelated theoretical error for a given observable.

        Args:
            obs (str): The observable for which the error is to be retrieved.

        Returns:
            float: The uncorrelated theoretical error.

        Raises:
            Exception: If the observable is not recognized or implemented.
        """
        obs = self._lfu_version(obs)

        if not obs in EWPOS:
            logging.debug(
                str(obs) + " is not implemented, no uncorr theory error")
            return 0.0

        if not hasattr(self, '_theoerr_{}'.format(obs)):
            logging.debug(
                'No uncorrelated theoretical error for {}'.format(obs))
            return 0.0

        return getattr(self, '_theoerr_{}'.format(obs))()

    def sin2theta_theoerr(self, obs):
        """
        Retrieves the theoretical error associated with the effective weak mixing angle (sin^2(theta)_eff).

        Args:
            obs (str): The observable for which the error is to be retrieved.

        Returns:
            float: The theoretical error for sin^2(theta)_eff.

        Raises:
            Exception: If the observable is not recognized or implemented.
        """
        obs = self._lfu_version(obs)

        if not obs in EWPOS:
            logging.debug(
                str(obs) + " is not implemented, no sin2theta theory error")
            return 0.0

        if not hasattr(self, '_s2terr_{}'.format(obs)):
            return 0.0

        return getattr(self, '_s2terr_{}'.format(obs))()

    def covariance(self, input_covariance, observable_mapping={}, add_exp_err=False, add_theo_err=True, add_para_err=True):
        """
        Computes the covariance matrix for observables based on input covariance of parameters.

        Args:
            input_covariance (dict): The covariance matrix of input parameters.
            observable_mapping (dict, optional): Mapping of observable names to input parameter names.
            add_exp_err (bool): Flag to include experimental errors in the covariance.
            add_theo_err (bool): Flag to include theoretical errors in the covariance.
            add_para_err (bool): Flag to include parameterization errors in the covariance.

        Returns:
            dict: The computed covariance matrix.

        Raises:
            Exception: If the input covariance matrix is improperly formatted or missing parameters.
        """
        for inp in self._SMinputs:
            if inp not in input_covariance:
                raise Exception(
                    "Need to provide {} in measurement covariance".format(inp))

        for n in input_covariance:
            for m in input_covariance:
                if m not in input_covariance[n]:
                    raise Exception(
                        "Bad covariance, missing {}:{}".format(n, m))

        names = input_covariance.keys()
        observable_map = {}

        # Map observables to inputs
        for n in names:
            if n in observable_mapping:
                observable_map[n] = observable_mapping[n]
            else:
                observable_map[n] = n

        covariance = {}
        for n in names:
            covariance[n] = {}
            for m in names:
                covariance[n][m] = input_covariance[n][m] if add_exp_err else 0.0

        if add_para_err:
            # Perform eigenvalue decomposition for input variations
            cov_matrix = Utils.toMatrix(
                self._SMinputs, self._SMinputs, input_covariance)
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
            shifts = []

            for value, vector in zip(eigenvalues, np.transpose(eigenvectors)):
                input_shifts = {}
                for x, name in zip(vector, self._SMinputs):
                    input_shifts[name] = value ** 0.5 * x
                output_shifts = {}
                for obs in self._outputs:
                    output_shifts[obs] = 0.0
                    for sminp in self._SMinputs:
                        output_shifts[obs] += self.derivative(
                            obs, sminp) * input_shifts[sminp]
                output_shifts.update(input_shifts)
                shifts.append(output_shifts)

            for shift in shifts:
                for n in names:
                    for m in names:
                        if observable_map[n] in shift and observable_map[m] in shift:
                            covariance[n][m] += float(shift[observable_map[n]]
                                                      * shift[observable_map[m]])

        if add_theo_err:
            for n in names:
                if n not in self.inputs():
                    covariance[n][n] += self.uncorr_theoerr(
                        observable_map[n]) ** 2
                for m in names:
                    if n not in self.inputs():
                        covariance[n][m] += self.sin2theta_theoerr(
                            observable_map[n]) * self.sin2theta_theoerr(observable_map[m])

        return dict(covariance)

    # returns the lepton flavour universal observable name (Al=Ae=Amu=Atau in the SM etc)
    def _lfu_version(self, obs):
        if obs in ['Ae', 'Amu', 'Atau']:
            return 'Al'
        if obs in ['AFBe', 'AFBmu', 'AFBtau']:
            return 'AFBl'
        if obs in ['Re', 'Rmu', 'Rtau']:
            return 'Rl'
        if obs in ['BrWe', 'BrWmu', 'BrWtau']:
            return 'BrWl'

        return obs

    # update Higgs mass and various related parameters used in interpolation formulas
    def _update_MH(self, MH):
        self._MH = MH
        self._LH = log(self._MH/self._MH_REF)
        self._dH = (self._MH/self._MH_REF)-1.
        self._dH_MW = log(self._MH/self._MH_REFMW)
        self._dh_MW = (self._MH/self._MH_REFMW)**2
        self._LH_old = log(self._MH/self._MH_REFMW)
        self._dH_old = (self._MH/self._MH_REFMW)**2

    # update top mass and various related parameters used in interpolation formulas
    def _update_mt(self, mt):
        self._mt = mt
        self._dt = (self._mt/self._MT_REF)**2-1
        self._dt_MW = (self._mt/self._MT_REFMW)**2-1
        self._dt_old = (self._mt/self._MT_REFMW)**2-1

    # update the strong coupling alpha_s and various related parameters used in interpolation formulas
    def _update_alphas(self, alphas):
        self._alphas = alphas
        self._dalphas = self._alphas/self._ALPHAS_REF-1
        self._dalphas_MW = self._alphas/self._ALPHAS_REFMW-1
        self._dalphas_old = self._alphas/self._ALPHAS_REFMW-1

    # update Z mass and various related parameters used in interpolation formulas
    def _update_MZ(self, MZ):
        self._MZ = MZ
        self._dZ = self._MZ/self._MZ_REF-1
        self._dZ_MW = (self._MZ/self._MZ_REFMW)-1
        self._dZ_old = (self._MZ/self._MZ_REFMW)-1

    # update W mass (only possible if MW taken as input)
    def _update_MW(self, MW):
        if self._MWscheme:
            if MW > 0:
                self._MW = MW
        else:
            raise Exception("Using alpha scheme, should not set MW")

    # update Delta alpha, the paraemter determining running alpha(0)->alpha(MZ) (only possible if MW not taken as input)
    def _update_Deltaalpha(self, Deltaalpha):
        if self._MWscheme:
            raise Exception("Using MW scheme, should not set alpha")
        elif Deltaalpha > 0:
            self._Deltaalpha = Deltaalpha
            self._update_Deltaalpha_common(Deltaalpha)

    # update various parameters related to Deltaalpha and used in interpolation formulas
    def _update_Deltaalpha_common(self, Deltaalpha):
        self._dDeltaalpha = Deltaalpha/self._DELTAALPHA_REF-1
        self._dDeltaalpha_MW = (Deltaalpha/self._DELTAALPHA_REFMW)-1
        self._dDeltaalpha_old = (Deltaalpha/self._DELTAALPHA_REFMW)-1

    # predict MW from Deltaalpha and other SM inputs
    def _MW_from_Deltaalpha(self):
        c = self._C_MW
        return c[0] - c[1]*self._dH_MW - c[2]*self._dH_MW**2 + c[3]*self._dH_MW**4 + c[4]*(self._dh_MW - 1) - c[5]*self._dDeltaalpha_MW + c[6]*self._dt_MW - c[7]*self._dt_MW**2 - c[8]*self._dH_MW*self._dt_MW + c[9]*self._dh_MW*self._dt_MW - c[10]*self._dalphas_MW + c[11]*self._dZ_MW

    # predict Deltaalpha from MW and other SM inputs
    def _Deltaalpha_from_MW(self):
        c = self._C_MW
        dDeltaalpha = (c[0] - self._MW - c[1]*self._dH_MW - c[2]*self._dH_MW**2 + c[3]*self._dH_MW**4 + c[4]*(self._dh_MW - 1) + c[6]*self._dt_MW -
                       c[7]*self._dt_MW**2 - c[8]*self._dH_MW*self._dt_MW + c[9]*self._dh_MW*self._dt_MW - c[10]*self._dalphas_MW + c[11]*self._dZ_MW)/c[5]
        return (1+dDeltaalpha)*self._DELTAALPHA_REFMW

    # Linearized SM input depedence of Deltaalpha (from MW interpolation formulas)
    def _dDeltaalpha_dMZ(self):
        c = self._C_MW
        return self._DELTAALPHA_REFMW/c[5] * c[11]/self._MZ_REFMW

    def _dDeltaalpha_dalphas(self):
        c = self._C_MW
        return self._DELTAALPHA_REFMW/c[5] * (- c[10])/self._ALPHAS_REFMW

    def _dDeltaalpha_dmt(self):
        c = self._C_MW
        return self._DELTAALPHA_REFMW/c[5] * (c[6] - c[7]*2*self._dt_MW - c[8]*self._dH_MW + c[9]*self._dh_MW) * 2*self._mt/self._MT_REFMW**2

    def _dDeltaalpha_dMH(self):
        c = self._C_MW
        return self._DELTAALPHA_REFMW/c[5] * (-c[1] - c[2]*2*self._dH_MW + c[3]*4*self._dH_MW**3 - c[8]*self._dt_MW + 2*c[4] + 2*c[9]*self._dt_MW)/self._MH

    def _dDeltaalpha_dMW(self):
        c = self._C_MW
        return -self._DELTAALPHA_REFMW/c[5]

    # Linearized SM input dependence of MW
    def _dMW_dDeltaalpha(self):
        c = self._C_MW
        return -c[5]/self._DELTAALPHA_REFMW

    def _dMW_dalphas(self):
        c = self._C_MW
        return -c[10]/self._ALPHAS_REFMW

    def _dMW_dmt(self):
        c = self._C_MW
        return (2*self._mt/self._MT_REFMW**2)*(c[6] - 2*c[7]*self._dt_MW - c[8]*self._dH_MW + c[9]*self._dh_MW)

    def _dMW_dMZ(self):
        c = self._C_MW
        return c[11]/self._MZ_REFMW

    def _dMW_dMH(self):
        c = self._C_MW
        return (-c[1] - 2*c[2]*self._dH_MW + 4*c[3]*self._dH_MW**3-c[8]*self._dt_MW + 2*c[4]*self._dH_MW + 2*c[9]*self._dH_MW*self._dt_MW)/self._MH

    # Linearized alphas dependence of EWPOs
    def _dGammaZ_dalphas(self):
        return self._dGammaObs_dalphas(*self._A_GAMMA_Z)/1e3

    def _dRl_dalphas(self):
        return self._dGammaObs_dalphas(*self._A_R_L)/1e3

    def _dRc_dalphas(self):
        return self._dGammaObs_dalphas(*self._A_R_C)/1e5

    def _dRb_dalphas(self):
        return self._dGammaObs_dalphas(*self._A_R_B)/1e5

    def _dsigmahad_dalphas(self):
        return self._dGammaObs_dalphas(*self._A_SIGMA_HAD)

    # Linearized Deltaalpha dependence of EWPOs
    def _dGammaZ_dDeltaalpha(self):
        return self._dGammaObs_dDeltaalpha(*self._A_GAMMA_Z)/1e3

    def _dRl_dDeltaalpha(self):
        return self._dGammaObs_dDeltaalpha(*self._A_R_L)/1e3

    def _dRc_dDeltaalpha(self):
        return self._dGammaObs_dDeltaalpha(*self._A_R_C)/1e5

    def _dRb_dDeltaalpha(self):
        return self._dGammaObs_dDeltaalpha(*self._A_R_B)/1e5

    def _dsigmahad_dDeltaalpha(self):
        return self._dGammaObs_dDeltaalpha(*self._A_SIGMA_HAD)

    # Linearized mt dependence of EWPOs
    def _dGammaZ_dmt(self):
        return self._dGammaObs_dmt(*self._A_GAMMA_Z)/1e3

    def _dRl_dmt(self):
        return self._dGammaObs_dmt(*self._A_R_L)/1e3

    def _dRc_dmt(self):
        return self._dGammaObs_dmt(*self._A_R_C)/1e5

    def _dRb_dmt(self):
        return self._dGammaObs_dmt(*self._A_R_B)/1e5

    def _dsigmahad_dmt(self):
        return self._dGammaObs_dmt(*self._A_SIGMA_HAD)

    # Linearized MZ dependence of EWPOs
    def _dGammaZ_dMZ(self):
        return self._dGammaObs_dMZ(*self._A_GAMMA_Z)/1e3

    def _dRl_dMZ(self):
        return self._dGammaObs_dMZ(*self._A_R_L)/1e3

    def _dRc_dMZ(self):
        return self._dGammaObs_dMZ(*self._A_R_C)/1e5

    def _dRb_dMZ(self):
        return self._dGammaObs_dMZ(*self._A_R_B)/1e5

    def _dsigmahad_dMZ(self):
        return self._dGammaObs_dMZ(*self._A_SIGMA_HAD)

    # Linearized MH dependence of EWPOs
    def _dGammaZ_dMH(self):
        return self._dGammaObs_dMH(*self._A_GAMMA_Z)/1e3

    def _dRl_dMH(self):
        return self._dGammaObs_dMH(*self._A_R_L)/1e3

    def _dRc_dMH(self):
        return self._dGammaObs_dMH(*self._A_R_C)/1e5

    def _dRb_dMH(self):
        return self._dGammaObs_dMH(*self._A_R_B)/1e5

    def _dsigmahad_dMH(self):
        return self._dGammaObs_dMH(*self._A_SIGMA_HAD)

    # various interpolation formulas -- see below for interpolation constants and references

    def _sin2theta_lbeff(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        res = s0+d1*self._LH+d2*self._LH**2+d3*self._LH**4+d4*self._dDeltaalpha+d5*self._dt+d6 * \
            self._dt**2+d7*self._dt*self._LH+d8*self._dalphas + \
            d9*self._dalphas*self._dt+d10*self._dZ
        return res*1e-4

    def _dsin2theta_lbeff_dDeltaalpha(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        res = d4/self._DELTAALPHA_REF
        return res*1e-4

    def _dsin2theta_lbeff_dalphas(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        res = (d8+d9*self._dt)/self._ALPHAS_REF
        return res*1e-4

    def _dsin2theta_lbeff_dmt(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        res = (2*self._mt/self._MT_REF**2) * \
            (d5+2*d6*self._dt+d7*self._LH+d9*self._dalphas)
        return res*1e-4

    def _dsin2theta_lbeff_dMZ(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        res = d10/self._MZ_REF
        return res*1e-4

    def _dsin2theta_lbeff_dMH(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        res = (d1+2*d2*self._LH+4*d3*self._LH**3+d7*self._dt)/self._MH
        return res*1e-4

    def _sin2theta_udeff(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return s0+d1*self._LH_old+d2*self._LH_old**2+d3*self._LH_old**4+d4*(self._dH_old**2-1)+d5*self._dDeltaalpha_old+d6*self._dt_old+d7*self._dt_old**2+d8*self._dt_old*(self._dH_old-1)+d9*self._dalphas_old+d10*self._dZ_old

    def _dsin2theta_udeff_dDeltaalpha(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return d5/self._DELTAALPHA_REFOLD

    def _dsin2theta_udeff_dalphas(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return d9/self._ALPHAS_REFOLD

    def _dsin2theta_udeff_dmt(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return (2*self._mt/self._MT_REFOLD**2)*(d6+2*d7*self._dt_old+d8*(self._dH_old-1))

    def _dsin2theta_udeff_dMZ(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return d10/self._MZ_REFOLD

    def _dsin2theta_udeff_dMH(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return (1./self._MH)*(d1+2*d2*self._LH_old+4*d3*self._LH_old**3)+(2./self._MH_REFOLD)*(2*d4*self._dH_old+d8*self._dt_old)

    def _sin2thetaleff(self):
        return self._sin2theta_lbeff(*self._D_SIN2THETA_LEFF)

    def _dsin2thetaleff_dDeltaalpha(self):
        return self._dsin2theta_lbeff_dDeltaalpha(*self._D_SIN2THETA_LEFF)

    def _dsin2thetaleff_dalphas(self):
        return self._dsin2theta_lbeff_dalphas(*self._D_SIN2THETA_LEFF)

    def _dsin2thetaleff_dmt(self):
        return self._dsin2theta_lbeff_dmt(*self._D_SIN2THETA_LEFF)

    def _dsin2thetaleff_dMZ(self):
        return self._dsin2theta_lbeff_dMZ(*self._D_SIN2THETA_LEFF)

    def _dsin2thetaleff_dMH(self):
        return self._dsin2theta_lbeff_dMH(*self._D_SIN2THETA_LEFF)

    def _sin2thetabeff(self):
        return self._sin2theta_lbeff(*self._D_SIN2THETA_BEFF)

    def _dsin2thetabeff_dmt(self):
        return self._dsin2theta_lbeff_dmt(*self._D_SIN2THETA_BEFF)

    def _dsin2thetabeff_dMZ(self):
        return self._dsin2theta_lbeff_dMZ(*self._D_SIN2THETA_BEFF)

    def _dsin2thetabeff_dMH(self):
        return self._dsin2theta_lbeff_dMH(*self._D_SIN2THETA_BEFF)

    def _dsin2thetabeff_dDeltaalpha(self):
        return self._dsin2theta_lbeff_dDeltaalpha(*self._D_SIN2THETA_BEFF)

    def _dsin2thetabeff_dalphas(self):
        return self._dsin2theta_lbeff_dalphas(*self._D_SIN2THETA_BEFF)

    def _sin2thetaueff(self):
        return self._sin2theta_udeff(*self._D_SIN2THETA_UEFF)

    def _dsin2thetaueff_dmt(self):
        return self._dsin2theta_udeff_dmt(*self._D_SIN2THETA_UEFF)

    def _dsin2thetaueff_dMZ(self):
        return self._dsin2theta_udeff_dMZ(*self._D_SIN2THETA_UEFF)

    def _dsin2thetaueff_dMH(self):
        return self._dsin2theta_udeff_dMH(*self._D_SIN2THETA_UEFF)

    def _dsin2thetaueff_dDeltaalpha(self):
        return self._dsin2theta_udeff_dDeltaalpha(*self._D_SIN2THETA_UEFF)

    def _dsin2thetaueff_dalphas(self):
        return self._dsin2theta_udeff_dalphas(*self._D_SIN2THETA_UEFF)

    def _sin2thetadeff(self):
        return self._sin2theta_udeff(*self._D_SIN2THETA_DEFF)

    def _dsin2thetadeff_dmt(self):
        return self._dsin2theta_udeff_dmt(*self._D_SIN2THETA_DEFF)

    def _dsin2thetadeff_dMZ(self):
        return self._dsin2theta_udeff_dMZ(*self._D_SIN2THETA_DEFF)

    def _dsin2thetadeff_dMH(self):
        return self._dsin2theta_udeff_dMH(*self._D_SIN2THETA_DEFF)

    def _dsin2thetadeff_dDeltaalpha(self):
        return self._dsin2theta_udeff_dDeltaalpha(*self._D_SIN2THETA_DEFF)

    def _dsin2thetadeff_dalphas(self):
        return self._dsin2theta_udeff_dalphas(*self._D_SIN2THETA_DEFF)

    # Asymmetries and their input parameter dependence
    # Calculated from effective sin2theta values

    def _Af(self, absQf, s2t):
        return (1-4*absQf*s2t)/(1-4*absQf*s2t+8*(absQf*s2t)**2)

    def _dAf_ds2t(self, absQf, s2t):
        return 16*absQf**2*s2t*(2*absQf*s2t-1)/(8*absQf**2*s2t**2-4*absQf*s2t+1)**2

    def _dAl_ds2t(self):
        absQf = 1
        s2t = self._sin2thetaleff()
        return self._dAf_ds2t(absQf, s2t)

    def _dAu_ds2t(self):
        absQf = 2./3.
        s2t = self._sin2thetaueff()
        return self._dAf_ds2t(absQf, s2t)

    def _dAd_ds2t(self):
        absQf = 1./3.
        s2t = self._sin2thetadeff()
        return self._dAf_ds2t(absQf, s2t)

    def _dAc_ds2t(self):
        absQf = 2./3.
        s2t = self._sin2thetaueff()
        return self._dAf_ds2t(absQf, s2t)

    def _dAs_ds2t(self):
        absQf = 1./3.
        s2t = self._sin2thetadeff()
        return self._dAf_ds2t(absQf, s2t)

    def _dAb_ds2t(self):
        absQf = 1./3.
        s2t = self._sin2thetabeff()
        return self._dAf_ds2t(absQf, s2t)

    def _Al(self):
        absQf = 1
        s2t = self._sin2thetaleff()
        return self._Af(absQf, s2t)

    def _dAl_dDeltaalpha(self):
        return self._dAl_ds2t()*self._dsin2thetaleff_dDeltaalpha()

    def _dAl_dalphas(self):
        return self._dAl_ds2t()*self._dsin2thetaleff_dalphas()

    def _dAl_dmt(self):
        return self._dAl_ds2t()*self._dsin2thetaleff_dmt()

    def _dAl_dMZ(self):
        return self._dAl_ds2t()*self._dsin2thetaleff_dMZ()

    def _dAl_dMH(self):
        return self._dAl_ds2t()*self._dsin2thetaleff_dMH()

    def _Ab(self):
        absQf = 1./3.
        s2t = self._sin2thetabeff()
        return self._Af(absQf, s2t)

    def _dAb_dDeltaalpha(self):
        return self._dAb_ds2t()*self._dsin2thetabeff_dDeltaalpha()

    def _dAb_dalphas(self):
        return self._dAb_ds2t()*self._dsin2thetabeff_dalphas()

    def _dAb_dmt(self):
        return self._dAb_ds2t()*self._dsin2thetabeff_dmt()

    def _dAb_dMZ(self):
        return self._dAb_ds2t()*self._dsin2thetabeff_dMZ()

    def _dAb_dMH(self):
        return self._dAb_ds2t()*self._dsin2thetabeff_dMH()

    def _Au(self):
        absQf = 2./3.
        s2t = self._sin2thetaueff()
        return self._Af(absQf, s2t)

    def _dAu_dDeltaalpha(self):
        return self._dAu_ds2t()*self._dsin2thetaueff_dDeltaalpha()

    def _dAu_dalphas(self):
        return self._dAu_ds2t()*self._dsin2thetaueff_dalphas()

    def _dAu_dmt(self):
        return self._dAu_ds2t()*self._dsin2thetaueff_dmt()

    def _dAu_dMZ(self):
        return self._dAu_ds2t()*self._dsin2thetaueff_dMZ()

    def _dAu_dMH(self):
        return self._dAu_ds2t()*self._dsin2thetaueff_dMH()

    def _Ad(self):
        absQf = 1./3.
        s2t = self._sin2thetadeff()
        return self._Af(absQf, s2t)

    def _dAd_dDeltaalpha(self):
        return self._dAd_ds2t()*self._dsin2thetadeff_dDeltaalpha()

    def _dAd_dalphas(self):
        return self._dAd_ds2t()*self._dsin2thetadeff_dalphas()

    def _dAd_dmt(self):
        return self._dAd_ds2t()*self._dsin2thetadeff_dmt()

    def _dAd_dMZ(self):
        return self._dAd_ds2t()*self._dsin2thetadeff_dMZ()

    def _dAd_dMH(self):
        return self._dAd_ds2t()*self._dsin2thetadeff_dMH()

    def _As(self):
        absQf = 1./3.
        s2t = self._sin2thetadeff()
        return self._Af(absQf, s2t)

    def _dAs_dDeltaalpha(self):
        return self._dAs_ds2t()*self._dsin2thetadeff_dDeltaalpha()

    def _dAs_dalphas(self):
        return self._dAs_ds2t()*self._dsin2thetadeff_dalphas()

    def _dAs_dmt(self):
        return self._dAs_ds2t()*self._dsin2thetadeff_dmt()

    def _dAs_dMZ(self):
        return self._dAs_ds2t()*self._dsin2thetadeff_dMZ()

    def _dAs_dMH(self):
        return self._dAs_ds2t()*self._dsin2thetadeff_dMH()

    def _Ac(self):
        absQf = 2./3.
        s2t = self._sin2thetaueff()
        return self._Af(absQf, s2t)

    def _dAc_dDeltaalpha(self):
        return self._dAc_ds2t()*self._dsin2thetaueff_dDeltaalpha()

    def _dAc_dalphas(self):
        return self._dAc_ds2t()*self._dsin2thetaueff_dalphas()

    def _dAc_dmt(self):
        return self._dAc_ds2t()*self._dsin2thetaueff_dmt()

    def _dAc_dMZ(self):
        return self._dAc_ds2t()*self._dsin2thetaueff_dMZ()

    def _dAc_dMH(self):
        return self._dAc_ds2t()*self._dsin2thetaueff_dMH()

    def _AFBl(self):
        return 3./4*self._Al()**2

    def _dAFBl_dMZ(self):
        return 6./4*self._Al()*self._dAl_dMZ()

    def _dAFBl_dMH(self):
        return 6./4*self._Al()*self._dAl_dMH()

    def _dAFBl_dDeltaalpha(self):
        return 6./4*self._Al()*self._dAl_dDeltaalpha()

    def _dAFBl_dalphas(self):
        return 6./4*self._Al()*self._dAl_dalphas()

    def _dAFBl_dmt(self):
        return 6./4*self._Al()*self._dAl_dmt()

    def _AFBb(self):
        return 3./4*self._Al()*self._Ab()

    def _dAFBb_dMZ(self):
        return 3./4*(self._Ab()*self._dAl_dMZ()+self._Al()*self._dAb_dMZ())

    def _dAFBb_dMH(self):
        return 3./4*(self._Ab()*self._dAl_dMH()+self._Al()*self._dAb_dMH())

    def _dAFBb_dDeltaalpha(self):
        return 3./4*(self._Ab()*self._dAl_dDeltaalpha()+self._Al()*self._dAb_dDeltaalpha())

    def _dAFBb_dalphas(self):
        return 3./4*(self._Ab()*self._dAl_dalphas()+self._Al()*self._dAb_dalphas())

    def _dAFBb_dmt(self):
        return 3./4*(self._Ab()*self._dAl_dmt()+self._Al()*self._dAb_dmt())

    def _AFBc(self):
        return 3./4*self._Al()*self._Ac()

    def _dAFBc_dMZ(self):
        return 3./4*(self._Ac()*self._dAl_dMZ()+self._Al()*self._dAc_dMZ())

    def _dAFBc_dMH(self):
        return 3./4*(self._Ac()*self._dAl_dMH()+self._Al()*self._dAc_dMH())

    def _dAFBc_dDeltaalpha(self):
        return 3./4*(self._Ac()*self._dAl_dDeltaalpha()+self._Al()*self._dAc_dDeltaalpha())

    def _dAFBc_dalphas(self):
        return 3./4*(self._Ac()*self._dAl_dalphas()+self._Al()*self._dAc_dalphas())

    def _dAFBc_dmt(self):
        return 3./4*(self._Ac()*self._dAl_dmt()+self._Al()*self._dAc_dmt())

    # Z boson partial widths and their input parameter dependence
    def _GammaObs(self, *a):
        return a[0]+a[1]*self._LH+a[2]*self._LH**2+a[3]*self._LH**4+a[4]*self._dH+a[5]*self._dt+a[6]*self._dt**2+a[7]*self._dt*self._LH+a[8]*self._dt*self._LH**2\
            + a[9]*self._dalphas+a[10]*self._dalphas**2+a[11]*self._dalphas*self._dH+a[12]*self._dalphas*self._dt+a[13]*self._dDeltaalpha+a[14]*self._dDeltaalpha*self._dH\
            + a[15]*self._dDeltaalpha*self._dt+a[16]*self._dZ

    def _dGammaObs_dalphas(self, *a):
        return (a[9]+2*a[10]*self._dalphas+a[11]*self._dH+a[12]*self._dt)/self._ALPHAS_REF

    def _dGammaObs_dDeltaalpha(self, *a):
        return (a[13]+a[14]*self._dH+a[15]*self._dt)/self._DELTAALPHA_REF

    def _dGammaObs_dMZ(self, *a):
        return (a[16])/self._MZ_REF

    def _dGammaObs_dMH(self, *a):
        return (a[1]+2*self._LH*a[2]+4*self._LH**3*a[3]+a[7]*self._dt+a[8]*self._dt*2*self._LH)/self._MH\
            + (a[4]+a[11]*self._dalphas+a[14]*self._dDeltaalpha)/self._MH_REF

    def _dGammaObs_dmt(self, *a):
        return (2*self._mt/self._MT_REF**2)*(a[5]+2*a[6]*self._dt+a[7]*self._LH+a[8]*self._LH**2+a[12]*self._dalphas+a[15]*self._dDeltaalpha)

    # Z boson partial widths, for right (left) handed fermions  with h=1 (h=-1)
    def _Gammae(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_E)/1e3
        else:
            gamma = self._GammaObs(*self._A_GAMMA_E)/1e3
            a = self._Al()
            return (1-h*a)*gamma/2.

    def _Gammamu(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_E)/1e3
        else:
            gamma = self._GammaObs(*self._A_GAMMA_E)/1e3
            a = self._Al()
            return (1-h*a)*gamma/2.

    def _Gammatau(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_TAU)/1e3
        else:
            gamma = self._GammaObs(*self._A_GAMMA_TAU)/1e3
            a = self._Al()
            return (1-h*a)*gamma/2.

    def _Gammanu(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_NU)/1e3
        elif h < 0:
            return 1.
        else:
            return 0.

    def _Gammau(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_U)/1e3
        else:
            gamma = self._GammaObs(*self._A_GAMMA_U)/1e3
            a = self._Au()
            return (1-h*a)*gamma/2.

    def _Gammad(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_D)/1e3
        else:
            gamma = self._GammaObs(*self._A_GAMMA_D)/1e3
            a = self._Ad()
            return (1-h*a)*gamma/2.

    def _Gammac(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_C)/1e3
        else:
            gamma = self._GammaObs(*self._A_GAMMA_C)/1e3
            a = self._Ac()
            return (1-h*a)*gamma/2.

    def _Gammas(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_D)/1e3
        else:
            gamma = self._GammaObs(*self._A_GAMMA_D)/1e3
            a = self._As()
            return (1-h*a)*gamma/2.

    def _Gammab(self, h=0):
        if h == 0:
            return self._GammaObs(*self._A_GAMMA_B)/1e3
        else:
            gamma = self._GammaObs(*self._A_GAMMA_B)/1e3
            a = self._Ab()
            return (1-h*a)*gamma/2.

    # W width and input parameter dependence

    def _GammaW(self):
        # https://arxiv.org/abs/1104.1769
        GF = 1.1663787e-5
        xs = (self._alphas-0.118)/0.003
        RW = 2.1940-0.16*(log(1+(23/self._MH)**2)-log(1+(23./100)**2))
        return 3.3904e-1*self.MW()**3*GF*(1+8.478*RW*1e-3+0.00065*xs)

    def _dGammaW_dalphas(self):
        # https://arxiv.org/abs/1104.1769
        GF = 1.1663787e-5
        xs = (self._alphas-0.118)/0.003
        RW = 2.1940-0.16*(log(1+(23/self._MH)**2)-log(1+(23./100)**2))
        return 3.3904e-1*self.MW()**3*GF*(0.00065/0.003) + self._dGammaW_dMW()*self._dMW_dalphas()

    def _dGammaW_dMH(self):
        # https://arxiv.org/abs/1104.1769
        GF = 1.1663787e-5
        xs = (self._alphas-0.118)/0.003
        RW = -0.16*log(1+(23/self._MH)**2)
        dRW_dMH = -0.16*1/(1+(23/self._MH)**2)*(-2)*23**2/self._MH**3
        return 3.3904e-1*self.MW()**3*GF*8.478e-3*dRW_dMH + self._dGammaW_dMW()*self._dMW_dMH()

    def _dGammaW_dMW(self):
        GF = 1.1663787e-5
        xs = (self._alphas-0.118)/0.003
        RW = 2.1940-0.16*(log(1+(23/self._MH)**2)-log(1+(23./100)**2))
        return 3.3904e-1*3*self.MW()**2*GF*(1+8.478*RW*1e-3+0.00065*xs)

    def _dGammaW_dDeltaalpha(self):
        return self._dGammaW_dMW()*self._dMW_dDeltaalpha()

    def _dGammaW_dmt(self):
        return self._dGammaW_dMW()*self._dMW_dmt()

    def _dGammaW_dMZ(self):
        return self._dGammaW_dMW()*self._dMW_dMZ()

    # W hadronic BR and input parameter dependence (TODO, not very relevant)
    def _BrWhad(self):
        return 0.67542

    def _dBrWhad_dalphas(self):
        return 0.

    def _dBrWhad_dMH(self):
        return 0.

    def _dBrWhad_dMW(self):
        return 0.

    def _dBrWhad_dDeltaalpha(self):
        return 0.

    def _dBrWhad_dmt(self):
        return 0.

    def _dBrWhad_dMZ(self):
        return 0.

    # W leptonic BR and input parameter dependence (TODO, not very relevant)
    def _BrWl(self):
        return 0.1081778

    def _dBrWl_dalphas(self):
        return 0.

    def _dBrWl_dMH(self):
        return 0.

    def _dBrWl_dMW(self):
        return 0.

    def _dBrWl_dDeltaalpha(self):
        return 0.

    def _dBrWl_dmt(self):
        return 0.

    def _dBrWl_dMZ(self):
        return 0.

    # theory uncertainties for specific observables (e.g. due to missing higher order corrections)
    def _theoerr_MW(self):
        if self._MWscheme:
            raise Exception("No theory error on MW in MW scheme")
        else:
            return self._THEOERR_M_W

    def _theoerr_Deltaalpha(self):
        if self._MWscheme:
            return self._dDeltaalpha_dMW()*self._THEOERR_M_W
        else:
            raise Exception("No theory error on alpha in alpha scheme")

    def _theoerr_GammaZ(self):
        return self._THEOERR_GAMMA_Z

    def _theoerr_Rl(self):
        return self._THEOERR_R_L

    def _theoerr_Rc(self):
        return self._THEOERR_R_C

    def _theoerr_Rb(self):
        return self._THEOERR_R_B

    def _theoerr_sigmahad(self):
        return self._THEOERR_SIGMA_HAD

    def _theoerr_Al(self):
        return 0.

    def _theoerr_Ab(self):
        return 0.

    def _theoerr_Ac(self):
        return 0.

    def _theoerr_AFBl(self):
        return 0.

    def _theoerr_AFBb(self):
        return 0.

    def _theoerr_AFBc(self):
        return 0.

    def _theoerr_GammaW(self):
        return 0.

    def _theoerr_BrWhad(self):
        return 0.

    def _theoerr_BrWl(self):
        return 0.

    def _theoerr_Rsin2thetamue(self):
        return 0.

    def _theoerr_sin2thetaleff(self):
        return 0.

    # common uncertainties due to sin2thetaleff theory uncertainty
    def _s2terr_sin2thetaleff(self):
        return self._THEOERR_SIN2THETA_LEFF

    def _s2terr_Al(self):
        return self._dAl_ds2t()*self._THEOERR_SIN2THETA_LEFF

    def _s2terr_Ab(self):
        return self._dAb_ds2t()*self._THEOERR_SIN2THETA_BEFF

    def _s2terr_Ac(self):
        return self._dAc_ds2t()*self._THEOERR_SIN2THETA_UDEFF

    def _s2terr_AFBl(self):
        return self._AFBl()*(self._s2terr_Al()/self._Al())*2

    def _s2terr_AFBb(self):
        return self._AFBb()*((self._s2terr_Al()/self._Al())**2 + (self._s2terr_Ab()/self._Ab())**2)**0.5

    def _s2terr_AFBc(self):
        return self._AFBc()*((self._s2terr_Al()/self._Al())**2 + (self._s2terr_Ac()/self._Ac())**2)**0.5

    # set all interpolation constants (DO NOT CHANGE (unless there are typos))
    def _set_interpol_constants(self):
        # reference values and interpolation constants from https://arxiv.org/abs/1906.08815
        self._ALPHAS_REF = 0.1184
        self._DELTAALPHA_REF = 0.059
        self._MZ_REF = 91.1876
        self._MH_REF = 125.7
        self._MT_REF = 173.2
        # corresponding to X0 and a1, a2, ... a16
        self._A_GAMMA_E = (83.983, -0.1202, -0.06919, 0.00383, 0.0597, 0.8037, -0.015, -
                           0.0195, 0.0032,  -0.0956, -0.0078, -0.0095, 0.25, -1.08, 0.056, -0.37, 286)
        self._A_GAMMA_TAU = (83.793, -0.1200, -0.06905, 0.00382, 0.0596, 0.8023, -0.015, -
                             0.0195, 0.0032,  -0.0954, -0.0078, -0.0094, 0.25, -1.08, 0.056, -0.37, 285)
        self._A_GAMMA_NU = (167.176, -0.1752, -0.1249, 0.00595, 0.1046, 1.253, -0.110, -
                            0.0232, 0.0064,  -0.187, -0.014, -0.014, 0.37, -0.085, 0.054, -0.30, 503)
        self._A_GAMMA_U = (299.994, -0.6152, -0.2771, 0.0174, 0.2341, 4.051, -
                           0.467, -0.0676, 0.017,  14.26, 1.6, -0.046, 1.82, -11.1, 0.16, -1.0, 1253)
        self._A_GAMMA_C = (299.918, -0.6152, -0.2771, 0.0174, 0.2340, 4.051, -
                           0.467, -0.0676, 0.017,  14.26, 1.6, -0.046, 1.82, -11.1, 0.16, -1.0, 1252)
        self._A_GAMMA_D = (382.829, -0.6685, -0.3322, 0.0193, 0.2792, 3.792, -
                           0.18, -0.0706, 0.020,  10.20, -2.4, -0.052, 0.71, -10.1, 0.16, -0.92, 1469)
        self._A_GAMMA_B = (375.890, -0.6017, -0.3158, 0.0190, 0.227, -2.174,
                           0.042, -0.027, 0.021,  10.53, -2.4, -0.056, 1.2, -10.1, 0.15, -0.95, 1458)
        self._A_GAMMA_Z = (2494.75, -4.055, -2.117, 0.122, 1.746, 19.68, -
                           1.63, -0.432, 0.12,  58.61, -4.0, -0.32, 8.1, -56.1, 1.1, -6.8, 9267)
        self._A_R_L = (20751.6, -8.112, -1.174, 0.155, 0.16, -37.59, -10.9,
                       1.27, 0.29,  732.30, -44, -0.61, 5.7, -358, -4.7, 37, 11649)
        self._A_R_C = (17222.2, -4.049, -0.749, 0.0832, 1.08, 98.956, -15.1, -
                       0.761, 0.080,  230.9, 125, 0.045, 36.9, -120, 1.2, -6.2, 3667)
        self._A_R_B = (21585.0, 4.904, 0.9149, -0.0535, -2.676, -292.21, 20.0,
                       1.97, -0.11,  -131.9, -84, -0.27, 4.4, 71.9, -0.77, -4.4, -1790)
        self._A_SIGMA_HAD = (41489.6, 0.408, -0.320, 0.0424, 1.32, 60.17,
                             16.3, -2.31, -0.19,  -579.58, 38, 0.010, 7.5, 85.2, 9.1, -68, -85957)
        self._THEOERR_R_L = 6e-3
        self._THEOERR_R_C = 5e-5
        self._THEOERR_R_B = 10e-5
        self._THEOERR_GAMMA_Z = 0.4e-3
        self._THEOERR_SIGMA_HAD = 6
        # corresponding to s0 and d1, d2, ... d10
        self._D_SIN2THETA_LEFF = (
            2314.64, 4.616, 0.539, -0.0737, 206, -25.71, 4.00, 0.288, 3.88, -6.49, -6560)
        self._D_SIN2THETA_BEFF = (
            2327.04, 4.638, 0.558, -0.07, 207., -9.554, 3.83, 0.179, 2.41, -8.24, -6630)
        self._THEOERR_SIN2THETA_BEFF = 5.3e-5  # not used
        self._THEOERR_SIN2THETA_LEFF = 4.3e-5

        # reference values and interpolation constants from https://arxiv.org/abs/hep-ph/0311148
        self._MH_REFMW = 100.
        self._MT_REFMW = 174.3
        self._MZ_REFMW = 91.1875
        self._DELTAALPHA_REFMW = 0.05907
        self._ALPHAS_REFMW = 0.119
        # corresponding to MW0 and c1, c2, ... c11
        self._C_MW = (80.3779, 0.05427, 0.008931, 0.0000882, 0.000161,
                      1.070, 0.5237, 0.0679, 0.00179, 0.0000664, 0.0795, 114.9)
        self._THEOERR_M_W = 4e-3

        # reference values and interpolation constants from https://arxiv.org/abs/hep-ph/0608099
        self._MH_REFOLD = 100.
        self._MT_REFOLD = 178.0
        self._MZ_REFOLD = 91.1876
        self._DELTAALPHA_REFOLD = 0.05907
        self._ALPHAS_REFOLD = 0.117
        self._D_SIN2THETA_DEFF = (0.2310286, 4.720e-4, 2.06e-5, 3.85e-6, -
                                  1.85e-6, 2.07e-2, -2.848e-3, 1.81e-4, -9.73e-6, 3.97e-4, -6.55e-1)
        self._D_SIN2THETA_UEFF = (0.2311395, 4.726e-4, 2.07e-5, 3.85e-6, -
                                  1.85e-6, 2.07e-2, -2.853e-3, 1.83e-4, -9.73e-6, 3.98e-4, -6.55e-1)
        self._THEOERR_SIN2THETA_UDEFF = 4.7e-5
