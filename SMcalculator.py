from enum import Enum
from math import log
import Utils
import numpy as np
import logging
logging.basicConfig(level=logging.INFO)

# Definition of input parameter schemes
class INPUTSCHEME(Enum):
    alpha = ['Deltaalpha', 'MZ', 'Gmu', 'MH', 'mt', 'alphas']
    MW = ['MW', 'MZ', 'Gmu', 'MH', 'mt', 'alphas']
    alphaMW = ['Deltaalpha', 'MW', 'MZ', 'MH', 'mt', 'alphas']
    sin2theta = ['sin2thetaleff', 'MZ', 'Gmu', 'MH', 'mt', 'alphas']
    alphasin2theta = ['Deltaalpha',
                      'sin2thetaleff', 'MZ', 'MH', 'mt', 'alphas']


# The EWPOs predicted or used as input by this class
EWPOS = ['MH', 'mt', 'alphas', 'MZ', 'Gmu', 'Deltaalpha', 'GammaZ', 'Rl', 'Re', 'Rmu', 'Rtau', 'Rc', 'Rb', 'sigmahad', 'Al', 'Ae', 'Amu', 'Atau',
         'AFBl', 'AFBe', 'AFBmu', 'AFBtau', 'AFBb', 'AFBc', 'Ab', 'Ac', 'MW', 'GammaW', 'BrWhad', 'BrWe', 'BrWmu', 'BrWtau', 'BrWl', 'sin2thetaleff']

# Default values of SM inputs
MH_DEFAULT = 125.2
MT_DEFAULT = 172.57
ALPHAS_DEFAULT = 0.118
MZ_DEFAULT = 91.1875
DELTAALPHA_DEFAULT = 0.05903
GMU_DEFAULT = 1.1663788e-5
MW_DEFAULT = 80.3692
SIN2THETALEFF_DEFAULT = 0.23153


class EWPOcalculator:
    """
    A class that calculates Standard Model predictions for electroweak precision observables (EWPOs)
    The class supports five "schemes" ('alpha', 'MW', 'sin2theta', 'alphaMW', and 'alphasin2theta') for parameter inputs.
    The full list of inputs in each scheme is defined in the enum INPUTSCHEME. Internally, all calculations are performed
    based on interpolation formulas corresponding to two-loop calculations in the  {alpha,MZ,Gmu} EW input parameter scheme.
    These predictions are inverted to obtain, e.g., predictions as a function of {alpha(MW),MZ,Gmu}, if {MW,MZ,Gmu} are given as input.
    The class provides functionality to compute various EWPOs, their derivatives, and theory and parametric covariance.
    """

    def __init__(self, MH=None, mt=None, alphas=None, MZ=None, Deltaalpha=None, Gmu=None, MW=None, sin2thetaleff=None, scheme=None, input_dict=None):
        """
        Initializes the EWPOcalculator with either a set of 6 input parameters  ({MH,mt,alphas} + {MZ,Gmu,Deeltaalpha}
        -- possibly replacing either Deeltaalpha or Gmu replaced by either MW or sin2theta) or a dictionary of the values
        + the scheme chosen. Units are GeV.

        Example usage:
        sm = EWPOcalculator(MH=125.,mt=173.,alphas=0.118,MZ=91.2,MW=80.4,Gmu=1663788e-5)
        or using default values and MW scheme:
        sm = EWPOcalculator(scheme=EWPOcalculator.INPUTSCHEME.MW)
        or using a dict (here everything but MW would be set to default)
        sm = EWPOcalculator(scheme=EWPOcalculator.INPUTSCHEME.MW,input_dict={'MW':80.4})

        Args:
            MH (float): Higgs mass.
            mt (float): Top quark mass.
            alphas (float): Strong coupling constant.
            MZ (float): Z boson mass.
            Deltaalpha (float): Electroweak alpha parameter.
            Gmu (float): Fermi coupling from muon decay Gmu.
            MW (float): W boson mass (optional).
            sin2thetaleff (float): Effective leptonic weak mixing angle sin^2(theta)_{l,eff}
            EWinputs (list): EW input parameters (optional)

        Raises:
            Exception: If the input values are inconsistent or missing required parameters.
        """

        # Set interpolation constants and initial values
        self._set_interpol_constants()

        # Set inputs (possibly to defaults) and possibly figure out scheme
        self._set_input_values(
            MH, mt, alphas, MZ, Deltaalpha, Gmu, MW, sin2thetaleff, scheme, input_dict)

    def reset(self):
        """
        Resets the EWPOcalculator to its original input values.
        """
        self.update(MH=self._MH_original, mt=self._mt_original, alphas=self._alphas_original,
                    MZ=self._MZ_original, Deltaalpha=self._Deltaalpha_original, Gmu=self._Gmu_original,
                    MW=self._MW_original, sin2thetaleff=self._sin2thetaleff_original)

    def update(self, MH=None, mt=None, alphas=None, MZ=None, Deltaalpha=None, Gmu=None, MW=None, sin2thetaleff=None, input_dict=None):
        """
        Updates the SM parameters and interpolates new values for observables based on input changes.

        Args:
            MH (float): Higgs mass [GeV].
            mt (float): Top quark mass [GeV].
            alphas (float): Strong coupling constant.
            MZ (float): Z boson mass [GeV].
            Deltaalpha (float): Electroweak alpha parameter.
            Gmu (float): Fermi coupling from muon decay Gmu [GeV^-2].
            MW (float): W boson mass [GeV].
            sin2thetaleff (float): Effective leptonic weak mixing angle sin^2(theta)_{l,eff}
            input_dict (dict): A dictinary, as an alternative to setting the other args.

        Raises:
            Exception: If the inputs are inconsistent (both dict and explicit args) or not compatible with EW input scheme.
        """
        if input_dict is not None:
            if not MH is None and mt is None and alphas is None and MZ is None and Deltaalpha is None and Gmu is None and MW is None:
                raise Exception("Update only with input_dict OR float values")
            if 'MH' in input_dict:
                MH = input_dict['MH']
            if 'mt' in input_dict:
                mt = input_dict['mt']
            if 'alphas' in input_dict:
                alphas = input_dict['alphas']
            if 'MZ' in input_dict:
                MZ = input_dict['MZ']
            if self._scheme == INPUTSCHEME.alpha:
                if 'Deltaalpha' in input_dict:
                    Deltaalpha = input_dict['Deltaalpha']
                if 'Gmu' in input_dict:
                    Gmu = input_dict['Gmu']
            if self._scheme == INPUTSCHEME.alphaMW:
                if 'Deltaalpha' in input_dict:
                    Deltaalpha = input_dict['Deltaalpha']
                if 'MW' in input_dict:
                    MW = input_dict['MW']
            if self._scheme == INPUTSCHEME.MW:
                if 'Gmu' in input_dict:
                    Gmu = input_dict['Gmu']
                if 'MW' in input_dict:
                    MW = input_dict['MW']
            if self._scheme == INPUTSCHEME.sin2theta:
                if 'Gmu' in input_dict:
                    Gmu = input_dict['Gmu']
                if 'sin2thetaleff' in input_dict:
                    sin2thetaleff = input_dict['sin2thetaleff']
            if self._scheme == INPUTSCHEME.alphasin2theta:
                if 'Deltaalpha' in input_dict:
                    Deltaalpha = input_dict['Deltaalpha']
                if 'sin2thetaleff' in input_dict:
                    sin2thetaleff = input_dict['sin2thetaleff']

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
            if self._scheme in (INPUTSCHEME.MW, INPUTSCHEME.sin2theta):
                raise Exception(
                    'Deltaalpha is not an input for this calculator')
            self._update_Deltaalpha(Deltaalpha)
        if Gmu is not None:
            if self._scheme in (INPUTSCHEME.alphaMW, INPUTSCHEME.alphasin2theta):
                raise Exception('Gmu is not an input for this calculator')
            self._update_Gmu(Gmu)
        if sin2thetaleff is not None:
            if self._scheme not in (INPUTSCHEME.sin2theta, INPUTSCHEME.alphasin2theta):
                raise Exception(
                    'sin2thetaleff is not an input for this calculator')
            self._update_sin2thetaleff(sin2thetaleff)
        if MW is not None:
            if self._scheme not in (INPUTSCHEME.MW, INPUTSCHEME.alphaMW):
                raise Exception('MW is not an input for this calculator')
            self._update_MW(MW)

        if self._scheme == INPUTSCHEME.MW:
            self._update_Deltaalpha_common(self._Deltaalpha_from_MW())
        if self._scheme == INPUTSCHEME.sin2theta:
            self._update_Deltaalpha_common(
                self._Deltaalpha_from_sin2thetaleff())
        if self._scheme == INPUTSCHEME.alphaMW:
            self._update_Gmu_common(self._Gmu_from_MW())
        if self._scheme == INPUTSCHEME.alphasin2theta:
            self._update_Gmu_common(self._Gmu_from_sin2thetaleff())

    def scheme(self):
        return self._scheme

    def MW(self):
        """
        Returns the value of the W boson mass.

        Returns:
            float: MW value if MW scheme is used, otherwise MW derived from Deltaalpha.
        """
        if self._scheme in (INPUTSCHEME.MW, INPUTSCHEME.alphaMW):
            return self._MW
        else:
            return self._MW_prediction()

    def Deltaalpha(self):
        """
        Returns the value of the running alpha parameter Deltaalpha.

        Returns:
            float: Deltaalpha value.
        """
        if self._scheme == INPUTSCHEME.MW:
            return self._Deltaalpha_from_MW()
        elif self._scheme == INPUTSCHEME.sin2theta:
            return self._Deltaalpha_from_sin2thetaleff()
        else:
            return self._Deltaalpha

    def Gmu(self):
        """
        Returns the value of the Fermi constant from muon decay Gmu.

        Returns:
            float: Gmu value.
        """
        if self._scheme == INPUTSCHEME.alphaMW:
            return self._Gmu_from_MW()
        elif self._scheme == INPUTSCHEME.alphasin2theta:
            return self._Gmu_from_sin2thetaleff()
        else:
            return self._Gmu

    def sin2thetaleff(self):
        """
        Returns the effective weak mixing angle (sin^2(theta)_eff).

        Returns:
            float: Effective weak mixing angle sin^2(theta)_eff.
        """
        if self._scheme in (INPUTSCHEME.sin2theta, INPUTSCHEME.alphasin2theta):
            return self._sin2thetaleff
        else:
            return self._sin2thetaleff_prediction()

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
        return self._GammaObs(*self._A_GAMMA_Z) / 1e3 + self._DGAMMAZ_DGMU*self._dGmu

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
        return self._GammaObs(*self._A_R_L) / 1e3 + self._DRL_DGMU*self._dGmu

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
        return self._GammaObs(*self._A_R_C) / 1e5 + self._DRC_DGMU*self._dGmu

    def Rb(self):
        """
        Returns the ratio of Z boson decays to bottom quarks over hadrons (Rb).

        Returns:
            float: Ratio of Z decays to bottom quarks.
        """
        return self._GammaObs(*self._A_R_B) / 1e5 + + self._DRB_DGMU*self._dGmu

    def sigmahad(self):
        """
        Returns the hadronic cross-section at the Z pole (sigma_had).

        Returns:
            float: Hadronic cross-section at the Z boson pole.
        """
        return self._GammaObs(*self._A_SIGMA_HAD) + self._DSIGMAHAD_DGMU*self._dGmu

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
        Retrieves the values for all observables in the EWPO calculator.

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

        # Update derived values that are input to interpolation formulas
        if self._scheme == INPUTSCHEME.MW:
            self._update_Deltaalpha_common(self._Deltaalpha_from_MW())
        elif self._scheme == INPUTSCHEME.sin2theta:
            self._update_Deltaalpha_common(
                self._Deltaalpha_from_sin2thetaleff())
        elif self._scheme == INPUTSCHEME.alphaMW:
            self._update_Gmu_common(self._Gmu_from_MW())
        elif self._scheme == INPUTSCHEME.alphasin2theta:
            self._update_Gmu_common(self._Gmu_from_sin2theta())

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

        if self._scheme == INPUTSCHEME.alpha:
            if hasattr(self, '_d{}_d{}'.format(obs, inp)):
                return getattr(self, '_d{}_d{}'.format(obs, inp))()
            else:
                logging.debug('No {} derivative of {}'.format(inp, obs))
                return 0
        else:
            # Chain rule for derivatives
            if self._scheme == INPUTSCHEME.MW:
                oldinp = 'Deltaalpha'
                newinp = 'MW'
            elif self._scheme == INPUTSCHEME.alphaMW:
                oldinp = 'Gmu'
                newinp = 'MW'
            elif self._scheme == INPUTSCHEME.sin2theta:
                oldinp = 'Deltaalpha'
                newinp = 'sin2thetaleff'
            elif self._scheme == INPUTSCHEME.alphasin2theta:
                oldinp = 'Gmu'
                newinp = 'sin2thetaleff'
            doldinp_dnewinp = 1/getattr(self, f'_d{newinp}_d{oldinp}')()
            if obs == oldinp:
                dobs_doldinp = 1.0
            elif hasattr(self, '_d{}_d{}'.format(obs, oldinp)):
                dobs_doldinp = getattr(
                    self, '_d{}_d{}'.format(obs, oldinp))()
            else:
                logging.debug('No {} derivative of {}'.format(inp, obs))
                dobs_doldinp = 0.0

            if inp == newinp:
                return dobs_doldinp * doldinp_dnewinp
            else:
                dobs_dinp = getattr(self, '_d{}_d{}'.format(obs, inp))(
                ) if hasattr(self, '_d{}_d{}'.format(obs, inp)) else 0.0
                if obs == oldinp:
                    return dobs_dinp
                else:
                    dnewinp_dinp = getattr(
                        self, '_d{}_d{}'.format(newinp, inp))()
                    return dobs_dinp - dobs_doldinp * doldinp_dnewinp * dnewinp_dinp

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
                    "Need to provide covariance for {}".format(inp))

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
        # Add theory error to covariance
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

    # set defaults
    def _set_input_values(self, MH, mt, alphas, MZ, Deltaalpha, Gmu, MW, sin2thetaleff, scheme, input_dict):

        # figure out scheme
        if scheme is None:
            if input_dict is not None:
                raise Exception(
                    "Need to set input 'scheme' along with 'input_dict'")
            nEWinputs = 0
            for o in Deltaalpha, Gmu, MW, sin2thetaleff:
                if o is not None:
                    nEWinputs += 1
            if nEWinputs < 1:
                raise Exception(
                    "Need to select an input scheme or set two of {Gmu, Deltaalpha, MW, sin2thetaleff}")
            elif nEWinputs < 2:
                if Gmu is None:
                    logging.info(f"Using default Gmu of {GMU_DEFAULT}")
                    Gmu = GMU_DEFAULT
                else:
                    raise Exception(
                        "Need to select an input scheme or set two of {Gmu, Deltaalpha, MW, sin2thetaleff}")
            elif nEWinputs > 2:
                raise Exception(
                    "Incompatible inputs, set only two of {Gmu, Deltaalpha, MW, sin2thetaleff}")

            if Deltaalpha is not None and Gmu is not None:
                scheme = INPUTSCHEME.alpha
            if MW is not None and Gmu is not None:
                scheme = INPUTSCHEME.MW
            if MW is not None and Deltaalpha is not None:
                scheme = INPUTSCHEME.alphaMW
            if sin2thetaleff is not None and Gmu is not None:
                scheme = INPUTSCHEME.sin2theta
            if sin2thetaleff is not None and alpha is not None:
                scheme = INPUTSCHEME.alphasin2theta

            logging.info(f"Calculation in {scheme.name} scheme")
        else:
            if not isinstance(scheme, INPUTSCHEME):
                if hasattr(INPUTSCHEME, scheme):
                    scheme = getattr(INPUTSCHEME, scheme)
                else:
                    raise Exception("Not a valid input scheme")
            if input_dict is not None:
                if not MH is None and mt is None and alphas is None and MZ is None and Deltaalpha is None and Gmu is None and MW is None:
                    raise Exception(
                        "Create calculator with either input_dict or float values but not both")
                MH = input_dict['MH'] if 'MH' in input_dict else None
                mt = input_dict['mt'] if 'mt' in input_dict else None
                MZ = input_dict['MZ'] if 'MZ' in input_dict else None
                alphas = input_dict['alphas'] if 'alphas' in input_dict else None

            if scheme == INPUTSCHEME.alpha:
                if input_dict is not None:
                    Deltaalpha = input_dict['Deltaalpha'] if 'Deltaalpha' in input_dict else None
                    Gmu = input_dict['Gmu'] if 'Gmu' in input_dict else None
                if Deltaalpha is None:
                    Deltaalpha = DELTAALPHA_DEFAULT
                if Gmu is None:
                    Gmu = GMU_DEFAULT
            if scheme == INPUTSCHEME.MW:
                if input_dict is not None:
                    MW = input_dict['MW'] if 'MW' in input_dict else None
                    Gmu = input_dict['Gmu'] if 'Gmu' in input_dict else None
                if MW is None:
                    MW = MW_DEFAULT
                if Gmu is None:
                    Gmu = GMU_DEFAULT
            if scheme == INPUTSCHEME.alphaMW:
                if input_dict is not None:
                    MW = input_dict['MW'] if 'MW' in input_dict else None
                    Deltaalpha = input_dict['Deltaalpha'] if 'Deltaalpha' in input_dict else None
                if Deltaalpha is None:
                    if MW is None:
                        MW = MW_DEFAULT
                    if Deltaalpha is None:
                        Deltaalpha = DELTAALPHA_DEFAULT
            if scheme == INPUTSCHEME.sin2theta:
                if input_dict is not None:
                    sin2thetaleff = input_dict['sin2thetaleff'] if 'sin2thetaleff' in input_dict else None
                    Gmu = input_dict['Gmu'] if 'Gmu' in input_dict else None
                if sin2thetaleff is None:
                    sin2thetaleff = SIN2THETALEFF_DEFAULT
                if Gmu is None:
                    Gmu = GMU_DEFAULT
            if scheme == INPUTSCHEME.alphasin2theta:
                if input_dict is not None:
                    sin2thetaleff = input_dict['sin2thetaleff'] if 'sin2thetaleff' in input_dict else None
                    Deltaalpha = input_dict['Deltaalpha'] if 'Deltaalpha' in input_dict else None
                if sin2thetaleff is None:
                    sin2thetaleff = SIN2THETALEFF_DEFAULT
                if Deltaalpha is None:
                    Deltaalpha = DELTAALPHA_DEFAULT
        if MH is None:
            MH = MH_DEFAULT
        if mt is None:
            mt = MT_DEFAULT
        if alphas is None:
            alphas = ALPHAS_DEFAULT
        if MZ is None:
            MZ = MZ_DEFAULT

        self._scheme = scheme

        self._SMinputs = scheme.value

        self._outputs = [x for x in EWPOS if x not in self._SMinputs]

        self._MH_original = MH
        self._mt_original = mt
        self._alphas_original = alphas
        self._MZ_original = MZ
        self._Deltaalpha_original = Deltaalpha
        self._Gmu_original = Gmu
        self._MW_original = MW
        self._sin2thetaleff_original = sin2thetaleff

        self.update(MH, mt, alphas, MZ, Deltaalpha, Gmu, MW, sin2thetaleff)

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
        self._LH_old = log(self._MH/self._MH_REFOLD)
        self._dH_old = (self._MH/self._MH_REFOLD)

    # update top mass and various related parameters used in interpolation formulas
    def _update_mt(self, mt):
        self._mt = mt
        self._dt = (self._mt/self._MT_REF)**2-1
        self._dt_MW = (self._mt/self._MT_REFMW)**2-1
        self._dt_old = (self._mt/self._MT_REFOLD)**2-1

    # update the strong coupling alpha_s and various related parameters used in interpolation formulas
    def _update_alphas(self, alphas):
        self._alphas = alphas
        self._dalphas = self._alphas/self._ALPHAS_REF-1
        self._dalphas_MW = self._alphas/self._ALPHAS_REFMW-1
        self._dalphas_old = self._alphas/self._ALPHAS_REFOLD-1

    # update Z mass and various related parameters used in interpolation formulas
    def _update_MZ(self, MZ):
        self._MZ = MZ
        self._dZ = self._MZ/self._MZ_REF-1
        self._dZ_MW = (self._MZ/self._MZ_REFMW)-1
        self._dZ_old = (self._MZ/self._MZ_REFOLD)-1

    # update Fermi coupling Gmu and various related parameters used in interpolation formulas
    # (only possible if Gmu taken as input)
    def _update_Gmu(self, Gmu):
        if self._scheme not in (INPUTSCHEME.alpha, INPUTSCHEME.sin2theta, INPUTSCHEME.MW):
            raise Exception("Gmu is not an input in this scheme")
        self._update_Gmu_common(Gmu)

    # update Fermi coupling Gmu and various related parameters used in interpolation formulas
    def _update_Gmu_common(self, Gmu):
        self._Gmu = Gmu
        self._dGmu = self._Gmu-self._GMU_REF

    # update Delta alpha, the paraemter determining running alpha(0)->alpha(MZ) (only possible if MW not taken as input)

    def _update_Deltaalpha(self, Deltaalpha):
        if self._scheme not in (INPUTSCHEME.alpha, INPUTSCHEME.alphaMW, INPUTSCHEME.alphasin2theta):
            raise Exception("alpha is not an input in this scheme")
        else:
            self._Deltaalpha = Deltaalpha
            self._update_Deltaalpha_common(Deltaalpha)

    # update various parameters related to Deltaalpha and used in interpolation formulas
    def _update_Deltaalpha_common(self, Deltaalpha):
        self._dDeltaalpha = Deltaalpha/self._DELTAALPHA_REF-1
        self._dDeltaalpha_MW = (Deltaalpha/self._DELTAALPHA_REFMW)-1
        self._dDeltaalpha_old = (Deltaalpha/self._DELTAALPHA_REFOLD)-1

    # update W mass (only possible if MW taken as input)
    def _update_MW(self, MW):
        if self._scheme not in (INPUTSCHEME.MW, INPUTSCHEME.alphaMW):
            raise Exception("MW is not an input in this scheme")
        else:
            self._MW = MW

    # update sin2thetaleff mass (only possible if sin2thetaleff taken as input)
    def _update_sin2thetaleff(self, sin2thetaleff):
        if self._scheme not in (INPUTSCHEME.sin2theta, INPUTSCHEME.alphasin2theta):
            raise Exception("sin2thetaleff is not an input in this scheme")
        else:
            self._sin2thetaleff = sin2thetaleff

    # predict MW from Deltaalpha and other SM inputs
    def _MW_prediction(self):
        c = self._C_MW
        return c[0] - c[1]*self._dH_MW - c[2]*self._dH_MW**2 + c[3]*self._dH_MW**4 + c[4]*(self._dh_MW - 1) - c[5]*self._dDeltaalpha_MW + c[6]*self._dt_MW - c[7]*self._dt_MW**2 - c[8]*self._dH_MW*self._dt_MW + c[9]*self._dh_MW*self._dt_MW - c[10]*self._dalphas_MW + c[11]*self._dZ_MW + self._DMW_DGMU*self._dGmu

    # predict sin2thetaleff from Deltaalpha and other SM inputs
    def _sin2thetaleff_prediction(self):
        return self._sin2theta_lbeff(*self._D_SIN2THETA_LEFF)

    # predict Deltaalpha from MW and other SM inputs
    def _Deltaalpha_from_MW(self):
        c = self._C_MW
        dDeltaalpha = (c[0] - self._MW - c[1]*self._dH_MW - c[2]*self._dH_MW**2 + c[3]*self._dH_MW**4 + c[4]*(self._dh_MW - 1) + c[6]*self._dt_MW -
                       c[7]*self._dt_MW**2 - c[8]*self._dH_MW*self._dt_MW + c[9]*self._dh_MW*self._dt_MW - c[10]*self._dalphas_MW + c[11]*self._dZ_MW + self._DMW_DGMU*self._dGmu)/c[5]
        return (1+dDeltaalpha)*self._DELTAALPHA_REFMW

    def _Deltaalpha_from_sin2thetaleff(self):
        s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10 = self._D_SIN2THETA_LEFF
        dDeltaalpha = -(s0+d1*self._LH+d2*self._LH**2+d3*self._LH**4+d5*self._dt+d6 * self._dt**2+d7*self._dt*self._LH+d8 *
                        self._dalphas + d9*self._dalphas*self._dt+d10*self._dZ+self._DSIN2THETA_DGMU*self._dGmu*1e4 - self._sin2thetaleff*1e4)/d4
        return (1+dDeltaalpha)*self._DELTAALPHA_REF

    def _Gmu_from_MW(self):
        c = self._C_MW
        dGmu = -(c[0] - c[1]*self._dH_MW - c[2]*self._dH_MW**2 + c[3]*self._dH_MW**4 + c[4]*(self._dh_MW - 1) - c[5]*self._dDeltaalpha_MW + c[6]*self._dt_MW -
                 c[7]*self._dt_MW**2 - c[8]*self._dH_MW*self._dt_MW + c[9]*self._dh_MW*self._dt_MW - c[10]*self._dalphas_MW + c[11]*self._dZ_MW - self._MW)/self._DMW_DGMU
        return (self._GMU_REF+dGmu)

    def _Gmu_from_sin2thetaleff(self):
        s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10 = self._D_SIN2THETA_LEFF
        dGmu = (self._sin2thetaleff - (s0+d1*self._LH+d2*self._LH**2+d3*self._LH**4+d4*self._dDeltaalpha+d5*self._dt+d6 *
                self._dt**2+d7*self._dt*self._LH+d8*self._dalphas + d9*self._dalphas*self._dt+d10*self._dZ)*1e-4)/self._DSIN2THETA_DGMU
        return (self._GMU_REF+dGmu)

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
        return (-c[1] - 2*c[2]*self._dH_MW + 4*c[3]*self._dH_MW**3 - c[8]*self._dt_MW)/self._MH + (c[4] + c[9]*self._dt_MW)*(2*self._MH/self._MH_REFMW**2)

    def _dMW_dGmu(self):
        return self._DMW_DGMU

    # Linearized SM input dependence of sin2theta
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

    def _dsin2thetaleff_dGmu(self):
        return self._DSIN2THETA_DGMU

    # Linearized SM input depedence of Deltaalpha
    # for the two schemes in which Deltaalpha is not an input

    def _dDeltaalpha_dMZ(self):
        if self._scheme == INPUTSCHEME.MW:
            return -self._dDeltaalpha_dMW()*self._dMW_dMZ()
        if self._scheme == INPUTSCHEME.sin2theta:
            return -self._dDeltaalpha_dsin2thetaleff()*self._dsin2thetaleff_dMZ()

    def _dDeltaalpha_dalphas(self):
        if self._scheme == INPUTSCHEME.MW:
            return -self._dDeltaalpha_dMW()*self._dMW_dalphas()
        if self._scheme == INPUTSCHEME.sin2theta:
            return -self._dDeltaalpha_dsin2thetaleff()*self._dsin2thetaleff_dalphas()

    def _dDeltaalpha_dmt(self):
        if self._scheme == INPUTSCHEME.MW:
            return -self._dDeltaalpha_dMW()*self._dMW_dmt()
        if self._scheme == INPUTSCHEME.sin2theta:
            return -self._dDeltaalpha_dsin2thetaleff()*self._dsin2thetaleff_dmt()

    def _dDeltaalpha_dMH(self):
        if self._scheme == INPUTSCHEME.MW:
            return -self._dDeltaalpha_dMW()*self._dMW_dMH()
        if self._scheme == INPUTSCHEME.sin2theta:
            return -self._dDeltaalpha_dsin2thetaleff()*self._dsin2thetaleff_dMH()

    def _dDeltaalpha_dGmu(self):
        if self._scheme == INPUTSCHEME.MW:
            return -self._dDeltaalpha_dMW()*self._dMW_dGmu()
        if self._scheme == INPUTSCHEME.sin2theta:
            return -self._dDeltaalpha_dsin2thetaleff()*self._dsin2thetaleff_dGmu()

    def _dDeltaalpha_dsin2thetaleff(self):
        if self._scheme == INPUTSCHEME.sin2theta:
            return 1/self._dsin2thetaleff_dDeltaalpha()

    def _dDeltaalpha_dMW(self):
        if self._scheme == INPUTSCHEME.MW:
            return 1/self._dMW_dDeltaalpha()

    # Linearized SM input depedence of Gmu
    # for the two schemes in which Gmu is not an input
    def _dGmu_dMZ(self):
        if self._scheme == INPUTSCHEME.alphasin2theta:
            return -self._dGmu_dsin2thetaleff()*self._dsin2thetaleff_dMZ()
        if self._scheme == INPUTSCHEME.alphaMW:
            return -self._dGmu_dMW()*self._dMW_dMZ()

    def _dGmu_dalphas(self):
        if self._scheme == INPUTSCHEME.alphasin2theta:
            return -self._dGmu_dsin2thetaleff()*self._dsin2thetaleff_dalphas()
        if self._scheme == INPUTSCHEME.alphaMW:
            return -self._dGmu_dMW()*self._dMW_dalphas()

    def _dGmu_dmt(self):
        if self._scheme == INPUTSCHEME.alphasin2theta:
            return -self._dGmu_dsin2thetaleff()*self._dsin2thetaleff_dmt()
        if self._scheme == INPUTSCHEME.alphaMW:
            return -self._dGmu_dMW()*self._dMW_dmt()

    def _dGmu_dMH(self):
        if self._scheme == INPUTSCHEME.alphasin2theta:
            return -self._dGmu_dsin2thetaleff()*self._dsin2thetaleff_dMH()
        if self._scheme == INPUTSCHEME.alphaMW:
            return -self._dGmu_dMW()*self._dMW_dMH()

    def _dGmu_dMW(self):
        if self._scheme == INPUTSCHEME.alphasin2theta:
            return -self._dGmu_dsin2thetaleff()*self._dsin2thetaleff_dMW()
        if self._scheme == INPUTSCHEME.alphaMW:
            return 1/self._dMW_dGmu()

    def _dGmu_dsin2thetaleff(self):
        if self._scheme == INPUTSCHEME.alphasin2theta:
            return 1/self._dsin2thetaleff_dGmu()

    def _dGmu_dDeltaalpha(self):
        if self._scheme == INPUTSCHEME.alphaMW:
            return -self._dGmu_dMW()*self._dMW_dDeltaalpha()
        if self._scheme == INPUTSCHEME.alphasin2theta:
            return -self._dGmu_dsin2thetaleff()*self._dsin2thetaleff_dDeltaalpha()

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

    # Linearized Gmu dependence of EWPOs
    def _dGammaZ_dGmu(self):
        return self._DGAMMAZ_DGMU

    def _dRl_dGmu(self):
        return self._DRL_DGMU

    def _dRc_dGmu(self):
        return self._DRC_DGMU

    def _dRb_dGmu(self):
        return self._DRB_DGMU

    def _dsigmahad_dGmu(self):
        return self._DSIGMAHAD_DGMU

    # various interpolation formulas -- see below for interpolation constants and references
    def _sin2theta_lbeff(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        res = s0+d1*self._LH+d2*self._LH**2+d3*self._LH**4+d4*self._dDeltaalpha+d5*self._dt+d6 * \
            self._dt**2+d7*self._dt*self._LH+d8*self._dalphas + \
            d9*self._dalphas*self._dt+d10*self._dZ
        return res*1e-4+self._DSIN2THETA_DGMU*self._dGmu

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

    def _dsin2theta_lbeff_dGmu(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return self._DSIN2THETA_DGMU

    def _sin2theta_udeff(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return s0+d1*self._LH_old+d2*self._LH_old**2+d3*self._LH_old**4+d4*(self._dH_old**2-1)+d5*self._dDeltaalpha_old+d6*self._dt_old+d7*self._dt_old**2+d8*self._dt_old*(self._dH_old-1)+d9*self._dalphas_old+d10*self._dZ_old + self._DSIN2THETA_DGMU*self._dGmu

    def _dsin2theta_udeff_dDeltaalpha(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return d5/self._DELTAALPHA_REFOLD

    def _dsin2theta_udeff_dalphas(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return d9/self._ALPHAS_REFOLD

    def _dsin2theta_udeff_dmt(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return (2*self._mt/self._MT_REFOLD**2)*(d6+2*d7*self._dt_old+d8*(self._dH_old-1))

    def _dsin2theta_udeff_dMZ(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return d10/self._MZ_REFOLD

    def _dsin2theta_udeff_dMH(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return (1./self._MH)*(d1+2*d2*self._LH_old+4*d3*self._LH_old**3)+(2*d4*self._dH_old+d8*self._dt_old)/self._MH_REFOLD

    def _dsin2theta_udeff_dGmu(self, s0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10):
        return self._DSIN2THETA_DGMU

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

    def _dsin2thetabeff_dGmu(self):
        return self._dsin2theta_lbeff_dGmu(*self._D_SIN2THETA_BEFF)

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

    def _dsin2thetaueff_dGmu(self):
        return self._dsin2theta_udeff_dGmu(*self._D_SIN2THETA_UEFF)

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

    def _dsin2thetadeff_dGmu(self):
        return self._dsin2theta_udeff_dGmu(*self._D_SIN2THETA_DEFF)

    # Asymmetries and their input parameter dependence
    # Calculated from effective sin2theta values

    def _Af(self, absQf, s2t):
        return (1-4*absQf*s2t)/(1-4*absQf*s2t+8*(absQf*s2t)**2)

    def _dAf_ds2t(self, absQf, s2t):
        return 16*absQf**2*s2t*(2*absQf*s2t-1)/(8*absQf**2*s2t**2-4*absQf*s2t+1)**2

    def _dAl_ds2t(self):
        absQf = 1
        s2t = self.sin2thetaleff()
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
        s2t = self.sin2thetaleff()
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

    def _dAl_dGmu(self):
        return self._dAl_ds2t()*self._dsin2thetaleff_dGmu()

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

    def _dAb_dGmu(self):
        return self._dAb_ds2t()*self._dsin2thetabeff_dGmu()

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

    def _dAu_dGmu(self):
        return self._dAu_ds2t()*self._dsin2thetaueff_dGmu()

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

    def _dAd_dGmu(self):
        return self._dAd_ds2t()*self._dsin2thetadeff_dGmu()

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

    def _dAs_dGmu(self):
        return self._dAs_ds2t()*self._dsin2thetadeff_dGmu()

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

    def _dAc_dGmu(self):
        return self._dAc_ds2t()*self._dsin2thetaueff_dGmu()

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

    def _dAFBl_dGmu(self):
        return 6./4*self._Al()*self._dAl_dGmu()

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

    def _dAFBb_dGmu(self):
        return 3./4*(self._Ab()*self._dAl_dGmu()+self._Al()*self._dAb_dGmu())

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

    def _dAFBc_dGmu(self):
        return 3./4*(self._Ac()*self._dAl_dGmu()+self._Al()*self._dAc_dGmu())

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
        GF = self._Gmu
        xs = (self._alphas-0.118)/0.003
        RW = 2.1940-0.16*(log(1+(23/self._MH)**2)-log(1+(23./100)**2))
        return 3.3904e-1*self.MW()**3*GF*(1+8.478*RW*1e-3+0.00065*xs)

    def _dGammaW_dalphas(self):
        # https://arxiv.org/abs/1104.1769
        GF = self._Gmu
        xs = (self._alphas-0.118)/0.003
        RW = 2.1940-0.16*(log(1+(23/self._MH)**2)-log(1+(23./100)**2))
        return 3.3904e-1*self.MW()**3*GF*(0.00065/0.003) + self._dGammaW_dMW()*self._dMW_dalphas()

    def _dGammaW_dMH(self):
        # https://arxiv.org/abs/1104.1769
        GF = self._Gmu
        xs = (self._alphas-0.118)/0.003
        RW = -0.16*log(1+(23/self._MH)**2)
        dRW_dMH = -0.16*1/(1+(23/self._MH)**2)*(-2)*23**2/self._MH**3
        return 3.3904e-1*self.MW()**3*GF*8.478e-3*dRW_dMH + self._dGammaW_dMW()*self._dMW_dMH()

    def _dGammaW_dMW(self):
        GF = self._Gmu
        xs = (self._alphas-0.118)/0.003
        RW = 2.1940-0.16*(log(1+(23/self._MH)**2)-log(1+(23./100)**2))
        return 3.3904e-1*3*self.MW()**2*GF*(1+8.478*RW*1e-3+0.00065*xs)

    def _dGammaW_dGmu(self):
        # https://arxiv.org/abs/1104.1769
        GF = self._Gmu
        xs = (self._alphas-0.118)/0.003
        RW = 2.1940-0.16*(log(1+(23/self._MH)**2)-log(1+(23./100)**2))
        return 3.3904e-1*self.MW()**3*(1+8.478*RW*1e-3+0.00065*xs) + self._dGammaW_dMW()*self._dMW_dGmu()

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

    def _dBrWhad_dGmu(self):
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

    def _dBrWl_dGmu(self):
        return 0.

    # theory uncertainties on parameters that can also serve as input
    def _theoerr_MW(self):
        if self._scheme in (INPUTSCHEME.MW, INPUTSCHEME.alphaMW):
            return 0.
        return self._THEOERR_M_W

    def _s2terr_sin2thetaleff(self):
        if self._scheme in (INPUTSCHEME.sin2theta, INPUTSCHEME.alphasin2theta):
            return 0.
        return self._THEOERR_SIN2THETA_LEFF

    def _theoerr_Deltaalpha(self):
        if self._scheme in (INPUTSCHEME.alpha, INPUTSCHEME.alphaMW, INPUTSCHEME.alphasin2theta):
            return 0.
        if self._scheme == INPUTSCHEME.MW:
            return self._dDeltaalpha_dMW()*self._THEOERR_M_W
        if self._scheme == INPUTSCHEME.sin2theta:
            return self._dDeltaalpha_dsin2thetaleff()*self._THEOERR_SIN2THETA_LEFF

    def _theoerr_Gmu(self):
        if self._scheme in (INPUTSCHEME.alpha, INPUTSCHEME.MW, INPUTSCHEME.sin2theta):
            return 0.
        if self._scheme == INPUTSCHEME.alphasin2theta:
            return self._dGmu_dsin2thetaleff()*self._THEOERR_SIN2THETA_LEFF
        if self._scheme == INPUTSCHEME.alphaMW:
            return self._dGmu_dMW()*self._THEOERR_M_W

    # common uncertainties due to sin2thetaleff theory uncertainty
    def _s2terr_Al(self):
        return self._dAl_ds2t()*self._s2terr_sin2thetaleff()

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

    # theory uncertainties for specific observables (e.g. due to missing higher order corrections)

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

        # Gmu dependence at one-loop level, to augment two-loop parametrizations
        self._GMU_REF = 1.166379e-5
        self._DGAMMAW_DGMU = 298374
        self._DGAMMAZ_DGMU = 289696
        self._DMW_DGMU = 1547204
        self._DRB_DGMU = -1268
        self._DRC_DGMU = 1704
        self._DRL_DGMU = 530555
        self._DSIGMAHAD_DGMU = -170001740
        self._DSIN2THETA_DGMU = -28708
