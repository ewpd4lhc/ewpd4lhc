import shutil
import os
import cmath
import sys
import yaml
import importlib
from subprocess import call
from math import pi, log10

DEFAULT_NEVENTS = 10000
OUTPUT_SUBDIR = 'Events/run_01/'
BANNER_FILE = 'run_01_tag_1_banner.txt'
INTEGRATED_WEIGHT_KEY = 'Integrated weight (pb)'
RESTRICTION = '-massless'
INITIAL_MG_COMMANDS = """set automatic_html_opening False
set notification_center False
"""
CMD_TEMPLATE = INITIAL_MG_COMMANDS+"""set auto_convert_model T
import model {model}{restriction}
generate {process} {order}
output {output} -f
    quit"""

PARAM_CARD_NAME = 'Cards/param_card.dat'
PARAM_CARD_DEFAULT = 'param_card_massless.dat'
MG5_EXECUTABLE = 'bin/mg5_aMC'

ZDECAYS = [
    ('e-', 'e+'), ('mu-', 'mu+'), ('ta-', 'ta+'),
    ('u', 'u~'), ('d', 'd~'), ('c', 'c~'),
    ('s', 's~'), ('b', 'b~'), ('ve', 've~'),
    ('vm', 'vm~'), ('vt', 'vt~')
]
WDECAYS = [
    ('e+', 've'), ('mu+', 'vm'), ('ta+', 'vt'),
    ('u', 'd~'), ('c', 's~')
]


D6LINKEY = 'd6lin'
D6QUADKEY = 'd6quad'
SMKEY = 'sm'
D6LINORDER = 'NP^2==1'
D6QUADORDER = 'NP==1'
SMORDER = 'NP==0'

ALPHA0 = 1/137.035999084

MODELS = [('SMEFTsim_top_MwScheme_UFO', 'cll1221 cHDD cHd cHl111 cHl122 cHl133 cHl311 cHl322 cHl333 cHe11 cHe22 cHe33 cHj1 cHj3 cHu cHWB cHQ3 cHQ1 cHbq'.split(), 'MW', 'noLFU'),
          ('SMEFTsim_top_alphaScheme_UFO',
           'cll1221 cHDD cHd cHl111 cHl122 cHl133 cHl311 cHl322 cHl333 cHe11 cHe22 cHe33 cHj1 cHj3 cHu cHWB cHQ3 cHQ1 cHbq'.split(), 'alpha', 'noLFU'),
          ('SMEFTsim_topU3l_MwScheme_UFO',
           'cll1 cHDD cHd cHl1 cHl3 cHe cHj1 cHj3 cHu cHWB cHQ3 cHQ1 cHbq'.split(), 'MW', 'LFU'),
          ('SMEFTsim_topU3l_alphaScheme_UFO',
           'cll1 cHDD cHd cHl1 cHl3 cHe cHj1 cHj3 cHu cHWB cHQ3 cHQ1 cHbq'.split(), 'alpha', 'LFU'),
          ('SMEFTsim_U35_MwScheme_UFO',
           'cll1 cHDD cHd cHl1 cHl3 cHe cHq1 cHq3 cHu cHWB'.split(), 'MW', 'LFU'),
          ('SMEFTsim_U35_alphaScheme_UFO',
           'cll1 cHDD cHd cHl1 cHl3 cHe cHq1 cHq3 cHu cHWB'.split(), 'alpha', 'LFU'),
          ]


def rel_para(para):
    return multiply_para(1/para[SMKEY], para)


def multiply_para(factor, para):
    res = {}
    res[D6LINKEY] = {}
    res[D6QUADKEY] = {}
    for x in para[D6LINKEY]:
        res[D6LINKEY][x] = para[D6LINKEY][x]*factor
    for x in para[D6QUADKEY]:
        res[D6QUADKEY][x] = para[D6QUADKEY][x]*factor
    res[SMKEY] = para[SMKEY]*factor
    return res


def add_paras(para1, para2):
    res = {}
    res[D6LINKEY] = {}
    res[D6QUADKEY] = {}
    for x in para1[D6LINKEY]:
        res[D6LINKEY][x] = para1[D6LINKEY][x]+para2[D6LINKEY][x]
    for x in para1[D6QUADKEY]:
        res[D6QUADKEY][x] = para1[D6QUADKEY][x]+para2[D6QUADKEY][x]
    res[SMKEY] = para1[SMKEY]+para2[SMKEY]
    return res


def add_dicts(d1, d2):
    res = {}
    for x in d1:
        res[x] = d1[x]+d2[x]
    return res


def multiply_dicts(d1, d2):
    res = {}
    for c1 in d1:
        for c2 in d2:
            c = c1+'*'+c2 if c1 < c2 else c2+'*'+c1
            if c in res:
                res[c] += d1[c1]*d2[c2]
            else:
                res[c] = d1[c1]*d2[c2]
    return res


def multiply_paras(para1, para2):
    sm1 = para1[SMKEY]
    sm2 = para2[SMKEY]
    rel_para1 = multiply_para(1/sm1, para1)
    rel_para2 = multiply_para(1/sm2, para2)
    res = add_paras(rel_para1, rel_para2)
    newd62 = multiply_dicts(rel_para1[D6LINKEY], rel_para2[D6LINKEY])
    res[D6QUADKEY] = add_dicts(res[D6QUADKEY], newd62)
    res = multiply_para(sm1*sm2, res)
    res[SMKEY] = sm1*sm2
    return res


def invert_para(para):
    sm = para[SMKEY]
    rel_para = multiply_para(1/sm, para)
    # 1/(1+x) ~ 1 -x + x^2
    res = multiply_para(-1, rel_para)
    newd62 = multiply_dicts(rel_para[D6LINKEY], rel_para[D6LINKEY])
    res[D6QUADKEY] = add_dicts(res[D6QUADKEY], newd62)
    res = multiply_para(1./sm, res)
    res[SMKEY] = 1./sm
    return res


def sqrt_para(para):
    sm = para[SMKEY]
    rel_para = multiply_para(1/sm, para)
    # sqrt(1+x) ~ 1 +x/2 - x^2/8
    res = multiply_para(0.5, rel_para)
    newd62 = multiply_dicts(rel_para[D6LINKEY], rel_para[D6LINKEY])
    newd62 = {k: -1./8*v for k, v in newd62.items()}
    res[D6QUADKEY] = add_dicts(res[D6QUADKEY], newd62)
    res = multiply_para(sm**0.5, res)
    res[SMKEY] = sm**0.5
    return res


def rnd(r, n=6):
    if D6LINKEY in r:
        cmax = 0.
        for c in r[D6LINKEY]:
            cmax = max(cmax, abs(r[D6LINKEY][c]))
        n -= int(round(log10(cmax)))

    for r1 in r:
        if isinstance(r[r1], dict):
            r[r1] = rnd(r[r1], n)
        else:
            r[r1] = round(r[r1], n)
    return r


class MGrunner:

    def __init__(self, mg5dir, tmpdir, outdir, SMparams={}):
        self.MG5_DIR = mg5dir
        self.TMP_DIR = tmpdir
        self.SMparams = SMparams

        # Ensure TMP_DIR exists
        os.makedirs(self.TMP_DIR, exist_ok=True)

        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)

        # Define fermion pairs and doublets
        self.ffbars = ZDECAYS
        self.fdoublets = WDECAYS

    def run_process(self, process, name, nevents, model, restriction, scanpars, doQuad=False):
        """
        Runs the process and manages the output for SM, linear, and quadratic results.
        """
        output = os.path.join(self.TMP_DIR, name)
        sm = self.run(process, output + '_sm', order=SMORDER, coefficients={},
                      nevents=nevents, model=model, restriction=restriction)
        lin = {}
        quad = {}

        for c in scanpars:
            lin[c] = self.run(
                process, output + f'_{c}_'+D6LINKEY, order=D6LINORDER,
                coefficients={c: '1'}, nevents=nevents, model=model, restriction=restriction
            )

        if doQuad:
            for c1 in scanpars:
                for c2 in scanpars:
                    if c2 < c1:
                        continue
                    quad[c1 + '*' + c2] = self.run(
                        process, output + f'_{c1}_{c2}_'+D6QUADKEY, order=D6QUADORDER,
                        coefficients={c1: '1', c2: '1'}, nevents=nevents, model=model, restriction=restriction
                    )
            for c1 in scanpars:
                for c2 in scanpars:
                    if c1 < c2:
                        quad[c1 + '*' + c2] = quad[c1 + '*' + c2] - \
                            quad[c1 + '*' + c1] - quad[c2 + '*' + c2]

        return sm, lin, quad

    def run_all(self, doQuad=False, doZ=True, doW=True, models=None):
        """
        Runs all processes for given models and generates output data for W and Z processes.
        """

        if models is None:
            all_models = [m[0] for m in MODELS]
        else:
            all_models = models

        scan_pars = None
        for model in all_models:
            for m in MODELS:
                if m[0] == model:
                    scanpars = m[1]
                    break
        assert scanpars is not None

        # Run W processes
        if doW:
            for f in self.fdoublets:
                path = os.path.join(self.outdir, model)
                os.makedirs(path, exist_ok=True)
                n = f'w{f[0]}{f[1]}'
                sm, lin, quad = self.run_process(
                    process=f'w+ > {f[0]} {f[1]}', name=n, nevents=DEFAULT_NEVENTS, model=model,
                    restriction=RESTRICTION, scanpars=scanpars, doQuad=doQuad
                )
                self.save_output(path, n, sm, lin, quad)

        # Run Z processes
        if doZ:
            for f in self.ffbars:
                for h in (1, -1):
                    h0 = 'R' if h > 0 else 'L'
                    h1 = 'R' if h < 0 else 'L'
                    p = f'z > {f[1]}{{{h1}}} {f[0]}{{{h0}}}'
                    n = f'Z{f[0]}{h0}'
                    path = os.path.join(self.outdir, model)
                    os.makedirs(path, exist_ok=True)
                    sm, lin, quad = self.run_process(
                        process=p, name=n, nevents=DEFAULT_NEVENTS, model=model,
                        restriction=RESTRICTION, scanpars=scanpars, doQuad=doQuad
                    )
                    self.save_output(path, n, sm, lin, quad)

    def save_output(self, path, name, sm, lin, quad):
        """
        Saves the simulation output to a YAML file.
        """
        data = {SMKEY: sm, D6LINKEY: lin}
        if len(quad) > 0:
            data[D6QUADKEY] = quad
        with open(f'{path}/{name}.yml', 'w') as file:
            yaml.dump(data, file)

    def run(self, process, output, order, coefficients, nevents, model, restriction):
        """
        Generates, prepares, launches, and parses the output of a process.
        """
        if os.path.exists(output):
            shutil.rmtree(output)
        self.generate(process, order, output, model, restriction)
        if not os.path.exists(output):
            return 0.
        self.prepare(output, model)
        self.launch(output, nevents, coefficients)
        return self.parse(output)

    def parse(self, output):
        """
        Parses the cross-section from the output banner file.
        """
        xs = 0.
        resultfile = os.path.join(output, OUTPUT_SUBDIR, BANNER_FILE)
        if not os.path.exists(resultfile):
            if os.path.exists(output):
                shutil.rmtree(output)
            return 0.
        with open(resultfile) as f:
            for line in f:
                if INTEGRATED_WEIGHT_KEY in line:
                    xs = float(line.split(':')[-1].strip())
        shutil.rmtree(output)
        return xs

    def generate(self, process, order, output, model, restriction):
        """
        Generates the process using MG5 with the given settings.
        """
        mg5_cmd = CMD_TEMPLATE.format(
            model=model, restriction=restriction, process=process, order=order, output=output
        )
        self.call_mg5(output, mg5_cmd)

    def prepare(self, output, model):
        """
        Prepares the parameter card by copying the default card for the model.
        """
        the_card = os.path.join(output, PARAM_CARD_NAME)
        default_card = os.path.join(
            self.MG5_DIR, 'models', model, PARAM_CARD_DEFAULT)
        shutil.copy(default_card, the_card)

    def launch(self, output, nevents, coefficients={}):
        """
        Launches the simulation with the specified settings.
        """
        mg5_cmd = INITIAL_MG_COMMANDS + \
            f"\nlaunch {output}\ndone\nset nevents {nevents}\nset use_syst False\n"
        mg5_cmd += '\n'.join([f"set {c} {coefficients[c]}" for c in coefficients])
        mg5_cmd += '\n' + \
            '\n'.join([f"set {p} {self.SMparams[p]}" for p in self.SMparams])
        mg5_cmd += '\ndone\nquit'
        self.call_mg5(output, mg5_cmd)

    def call_mg5(self, name, cmd):
        """
        Calls MG5 with a given command string.
        """
        cmdcard = os.path.join(self.TMP_DIR, name.replace('/', '_') + '.cmd')
        with open(cmdcard, 'w') as f:
            f.write(cmd)
        mg5_exec = os.path.join(self.MG5_DIR, MG5_EXECUTABLE)
        call([mg5_exec, cmdcard])
        os.remove(cmdcard)


class CalcPara:
    def __init__(self, paradir, smdata=None, rescale=False, MWscheme=False):
        self.gamma = {}
        self.gammaW = {}
        self.rescale = rescale
        for f in [x[0] for x in ZDECAYS]:
            for h in 'L', 'R':
                if os.path.exists(paradir+'/Z'+f+h+'.yml'):
                    with open(paradir+'/Z'+f+h+'.yml', 'r') as inf:
                        self.gamma[f+h] = yaml.safe_load(inf)
                        if not D6QUADKEY in self.gamma[f+h]:
                            self.gamma[f+h][D6QUADKEY] = {}
                else:
                    print('path does not exist', paradir+'/Z'+f+h+'.yml')
        for ff in [''.join(x) for x in WDECAYS]:
            if os.path.exists(paradir+'/w'+ff+'.yml'):
                with open(paradir+'/w'+ff+'.yml', 'r') as inf:
                    self.gammaW[ff] = yaml.safe_load(inf)
                    if not D6QUADKEY in self.gammaW[ff]:
                        self.gammaW[ff][D6QUADKEY] = {}
            else:
                print('path does not exist', paradir+'/Z'+f+h+'.yml')
        if os.path.exists(paradir+'/MW.yml'):
            with open(paradir+'/MW.yml', 'r') as inf:
                self.MWpara = yaml.safe_load(inf)
        else:
            self.MWpara = None
        if os.path.exists(paradir+'/Deltaalpha.yml'):
            with open(paradir+'/Deltaalpha.yml', 'r') as inf:
                self.alpha = yaml.safe_load(inf)
        else:
            self.alpha = None

        if smdata is not None:

            thisdir = os.path.dirname(__file__)
            sys.path.insert(0, '../'+thisdir)
            import SMcalculator

            if not os.path.exists(smdata):
                raise Exception(str(smdata)+' not found')
            with open(smdata, 'r') as f:
                data = yaml.safe_load(f)
            if MWscheme:
                self.sm = SMcalculator.EWPOcalculator(input_dict={
                                                    o: data['measured'][o]['central'] for o in SMcalculator.INPUTS_MWSCHEME}, scheme='MW')
            else:
                self.sm = SMcalculator.EWPOcalculator(input_dict={
                                                    o: data['measured'][o]['central'] for o in SMcalculator.INPUTS_ALPHASCHEME}, scheme='alpha')
        else:
            self.sm = None

    def GammaPara(self, f, h):
        if h == 0:
            res = {}
            res[SMKEY] = self.gamma[f+'L'][SMKEY]+self.gamma[f+'R'][SMKEY]
            for o in [D6LINKEY, D6QUADKEY]:
                res[o] = {}
                for x in self.gamma[f+'L'][o]:
                    res[o][x] = self.gamma[f+'L'][o][x]+self.gamma[f+'R'][o][x]
            return res
        elif h == -1:
            return self.gamma[f+'L']
        elif h == 1:
            return self.gamma[f+'R']

    def GammaZpol(self, f, h=0):
        res = {}
        para = self.GammaPara(f, h)
        if not self.rescale:
            return para
        k = self.GammaSM(f, h)/para[SMKEY] if para[SMKEY] != 0 else 0
        for x in para:
            if isinstance(para[x], dict):
                res[x] = {}
                for y in para[x]:
                    res[x][y] = para[x][y]*k
            else:
                res[x] = k*para[x]
        return res

    def GammaW(self, into=[''.join(x) for x in WDECAYS]):
        res = {}
        for ff in into:
            para = self.gammaW[ff]
            for x in para:
                if isinstance(para[x], dict):
                    if not x in res:
                        res[x] = {}
                    for y in para[x]:
                        if y in res[x]:
                            res[x][y] += para[x][y]
                        else:
                            res[x][y] = para[x][y]
                else:
                    if x in res:
                        res[x] += para[x]
                    else:
                        res[x] = para[x]
        return res

    def GammaWhad(self):
        return self.GammaW(['ud~', 'cs~'])

    def GammaWe(self):
        return self.GammaW(['e+ve'])

    def GammaWmu(self):
        return self.GammaW(['mu+vm'])

    def GammaWtau(self):
        return self.GammaW(['ta+vt'])

    def GammaZe(self):
        return self.GammaZpol('e-')

    def GammaZmu(self):
        return self.GammaZpol('mu-')

    def GammaZtau(self):
        return self.GammaZpol('ta-')

    def BrWhad(self):
        return multiply_paras(self.GammaWhad(), invert_para(self.GammaW()))

    def BrWe(self):
        return multiply_paras(self.GammaWe(), invert_para(self.GammaW()))

    def BrWmu(self):
        return multiply_paras(self.GammaWmu(), invert_para(self.GammaW()))

    def BrWtau(self):
        return multiply_paras(self.GammaWtau(), invert_para(self.GammaW()))

    def RWmue(self):
        return multiply_paras(self.GammaWmu(), invert_para(self.GammaWe()))

    def RWtaue(self):
        return multiply_paras(self.GammaWtau(), invert_para(self.GammaWe()))

    def RWtaumu(self):
        return multiply_paras(self.GammaWtau(), invert_para(self.GammaWmu()))

    def RZmue(self):
        return multiply_paras(self.GammaZmu(), invert_para(self.GammaZe()))

    def RWZmue(self):
        return multiply_paras(self.RWmue(), invert_para(sqrt_para(self.RZmue())))

    def MW(self):
        return self.MWpara

    def Deltaalpha(self):
        return self.alpha

    def sin2thetaleff(self, f='e-'):
        # calculate from Al:
        # ((1-Al**2)**0.5+Al-1)/(4*Al)
        Alpara = self.A(f)
        Al2para = multiply_paras(Alpara, Alpara)
        res = multiply_para(-1, Al2para)
        res[SMKEY] = res[SMKEY]+1
        res = sqrt_para(res)
        res = add_paras(res, Alpara)
        res[SMKEY] = res[SMKEY]-1
        res = multiply_para(0.25, res)
        res = multiply_paras(res, invert_para(Alpara))
        return res

    def Rsin2thetamue(self):
        return multiply_paras(self.sin2thetaleff('mu-'), invert_para(self.sin2thetaleff('e-')))

    def GammaSM(self, f, h):
        f = f.replace('-', '').replace('ta', 'tau').replace('vt',
                                                            'nu').replace('vm', 'nu').replace('ve', 'nu')
        return getattr(self.sm, '_Gamma'+f)(h) if self.sm is not None else None

    def A(self, f):
        # L-R / L + R
        dinv = invert_para(self.GammaZpol(f))
        n = add_paras(self.GammaZpol(f, -1),
                      multiply_para(-1, self.GammaZpol(f, +1)))
        res = multiply_paras(dinv, n)
        return res

    def ASM(self, f):
        return getattr(self.sm, 'A'+f.replace('e-', 'l').replace('mu-', 'l').replace('ta-', 'l'))()

    def AFB(self, f):
        res = multiply_para(3/4., multiply_paras(self.A('e-'), self.A(f)))
        return res

    def Gammahad(self):
        return self.Gammasum(['u', 'd', 'c', 's', 'b'])

    def GammaZ(self):
        return self.Gammasum(['u', 'd', 'c', 's', 'b', 'e-', 'mu-', 'ta-', 've', 'vm', 'vt'])

    def Gammasum(self, fs):
        res = {}
        for f in fs:
            para = self.GammaZpol(f)
            for x in para:
                if isinstance(para[x], dict):
                    if not x in res:
                        res[x] = {}
                    for y in para[x]:
                        if y in res[x]:
                            res[x][y] += para[x][y]
                        else:
                            res[x][y] = para[x][y]
                else:
                    if x in res:
                        res[x] += para[x]
                    else:
                        res[x] = para[x]
        return res

    def R(self, f):
        if f in ['e-', 'mu-', 'ta-', 'l-', 'l', 'ell']:
            flip = True
        else:
            flip = False
        # Gammahad/Gammaf
        res = multiply_paras(self.GammaZpol(f), invert_para(self.Gammahad()))
        if flip:
            res = invert_para(res)
        return res

    def sigmahad(self):
        mZ = 91.1875
        if self.sm is not None:
            mZ = self.sm.MZ()
        pb_per_GeV2 = 0.389379e9
        prefact = 12*pi/mZ**2*pb_per_GeV2
        res = multiply_paras(multiply_paras(invert_para(multiply_paras(
            self.GammaZ(), self.GammaZ())), self.Gammahad()), self.GammaZpol('e-'))
        res = multiply_para(prefact, res)
        return res


class MWalphaShifts:

    def __init__(self, modelpath, outdir):
        self.modelpath = modelpath
        self.outdir = outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    def MW(self, coeffs):
        import parameters
        importlib.reload(parameters)

        for c in coeffs:
            globals()[c] = coeffs[c]
        LambdaSMEFT = 1e3
        Gf = parameters.Gf.value
        aEW = parameters.aEW.value
        MZ = parameters.MZ.value

        sth2 = eval(parameters.sth2.value)
        sth = eval(parameters.sth.value)
        cth = eval(parameters.cth.value)
        MWsm = eval(parameters.MWsm.value)
        vevhat = eval(parameters.vevhat.value)
        dMZ2 = eval(parameters.dMZ2.value)
        dGf = eval(parameters.dGf.value)
        dgw = eval(parameters.dgw.value)
        dMW = eval(parameters.dMW.value)
        ee = eval(parameters.ee.value)
        g1 = eval(parameters.g1.value)
        gw = eval(parameters.gw.value)
        cth2 = cth**2
        DeltaAlpha = -2*g1*gw/(g1**2+gw**2)*cHWB*vevhat**2/LambdaSMEFT**2
        c2th = cth**2-sth2
        dMW_ = MWsm*(1/2.*(-cth2/c2th*dMZ2+sth2/c2th*(DeltaAlpha-dGf)))
        for c in coeffs:
            globals().pop(c)

        return dMW.real

    def MWSM(self):
        import parameters
        importlib.reload(parameters)

        LambdaSMEFT = 1e3
        Gf = parameters.Gf.value
        aEW = parameters.aEW.value
        MZ = parameters.MZ.value

        sth2 = eval(parameters.sth2.value)
        sth = eval(parameters.sth.value)
        cth = eval(parameters.cth.value)
        return eval(parameters.MWsm.value).real

    def alphaSM(self):
        import parameters
        importlib.reload(parameters)
        LambdaSMEFT = 1e3
        Gf = parameters.Gf.value
        MZ = parameters.MZ.value
        MW = parameters.MW.value

        return eval(parameters.aEW.value).real

    def alpha(self, coeffs):
        import parameters
        importlib.reload(parameters)

        for c in coeffs:
            globals()[c] = coeffs[c]
        LambdaSMEFT = 1e3
        Gf = parameters.Gf.value
        MZ = parameters.MZ.value
        MW = parameters.MW.value
        aEW = eval(parameters.aEW.value)
        sth2 = eval(parameters.sth2.value)
        sth = eval(parameters.sth.value)
        cth = eval(parameters.cth.value)
        cth2 = cth**2
        MWsm = eval(parameters.MWsm.value)
        vevhat = eval(parameters.vevhat.value)
        ee = eval(parameters.ee.value)
        g1 = eval(parameters.g1.value)
        gw = eval(parameters.gw.value)

        dMZ2 = eval(parameters.dMZ2.value)
        dGf = eval(parameters.dGf.value)
        dgw = eval(parameters.dgw.value)
        dg1 = eval(parameters.dg1.value)
        DeltaAlpha = -2*g1*gw/(g1**2+gw**2)*cHWB*vevhat**2/LambdaSMEFT**2
        alpha = 1/(4*pi)*g1**2*gw**2/(g1**2+gw**2) * \
            (2*cth2*dg1+2*sth2*dgw+DeltaAlpha)
        for c in coeffs:
            globals().pop(c)
        return alpha.real

    def run(self, models=None):
        if models is None:
            models = [m[0] for m in MODELS]
        alphamodels = []
        alphapars = []
        MWmodels = []
        MWpars = []
        for m in MODELS:
            if not m[0] in models:
                continue
            if m[2] == 'alpha':
                alphamodels.append(m[0])
                alphapars.append(m[1])
            elif m[2] == 'MW':
                MWmodels.append(m[0])
                MWpars.append(m[1])
        self.runMW(alphamodels, alphapars)
        self.runalpha(MWmodels, MWpars)

    def runMW(self, models, pars):
        for model, par in zip(models, pars):
            if not os.path.exists(self.modelpath+'/'+model):
                raise Exception('{} does not exist'.format(
                    self.modelpath+'/'+model))
            sys.path.insert(1, self.modelpath+model)

            result = {}
            result[SMKEY] = self.MWSM()
            result[D6LINKEY] = {}
            result[D6QUADKEY] = {}
            for p in par:
                d = {}
                for p1 in par:
                    d[p1] = 0
                d[p] = 1
                result[D6LINKEY][p] = self.MW(d)
            sys.path.remove(self.modelpath+model)
            path = self.outdir+'/'+model
            if not os.path.exists(path):
                os.makedirs(path)
            n = 'MW'
            with open('{}/{}.yml'.format(path, n), 'w') as file:
                yaml.dump(result, file)

    def runalpha(self, models, pars):
        alpha0 = ALPHA0
        for model, par in zip(models, pars):
            sys.path.insert(1, self.modelpath+model)

            result = {}
            result[SMKEY] = self.alphaSM()/alpha0-1
            result[D6LINKEY] = {}
            result[D6QUADKEY] = {}
            for p in split():
                d = {}
                for p1 in split():
                    d[p1] = 0
                d[p] = 1
                result[D6LINKEY][p] = self.alpha(d)/alpha0
            sys.path.remove(self.modelpath+model)
            path = self.outdir+'/'+model
            if not os.path.exists(path):
                os.makedirs(path)
            n = 'Deltaalpha'
            with open('{}/{}.yml'.format(path, n), 'w') as file:
                yaml.dump(result, file)


class CalcParaWrapper():
    def __init__(self, paradir, outdir, smdata=None, rescale=False):
        self.paradir = paradir
        self.outdir = outdir
        self.smdata = smdata
        self.rescale = rescale
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    def run(self, models=None):
        if models is None:
            models = [m[0] for m in MODELS]
        for m in MODELS:
            if not m[0] in models:
                continue
            MWscheme = m[2] == 'MW'
            modelname = m[0]
            lfu = m[3] == 'LFU'
            c = CalcPara(self.paradir+'/'+modelname, self.smdata,
                         rescale=self.rescale, MWscheme=MWscheme)

            with open(self.outdir+'/'+modelname+'.yml', 'w') as outf:
                yaml.dump({'GammaZ': rnd(c.GammaZ())}, outf)
                if lfu:
                    yaml.dump({'Rl': rnd(c.R('e-'))}, outf)
                else:
                    yaml.dump({'Re': rnd(c.R('e-'))}, outf)
                    yaml.dump({'Rmu': rnd(c.R('mu-'))}, outf)
                    yaml.dump({'Rtau': rnd(c.R('ta-'))}, outf)
                yaml.dump({'Rc': rnd(c.R('c'))}, outf)
                yaml.dump({'Rb': rnd(c.R('b'))}, outf)
                yaml.dump({'sigmahad': rnd(c.sigmahad())}, outf)
                if lfu:
                    yaml.dump({'Al': rnd(c.A('e-'))}, outf)
                else:
                    yaml.dump({'Ae': rnd(c.A('e-'))}, outf)
                    yaml.dump({'Amu': rnd(c.A('mu-'))}, outf)
                    yaml.dump({'Atau': rnd(c.A('ta-'))}, outf)

                yaml.dump({'Ab': rnd(c.A('b'))}, outf)
                yaml.dump({'Ac': rnd(c.A('c'))}, outf)
                if lfu:
                    yaml.dump({'AFBl': rnd(c.AFB('e-'))}, outf)
                else:
                    yaml.dump({'AFBe': rnd(c.AFB('e-'))}, outf)
                    yaml.dump({'AFBmu': rnd(c.AFB('mu-'))}, outf)
                    yaml.dump({'AFBtau': rnd(c.AFB('ta-'))}, outf)

                yaml.dump({'AFBb': rnd(c.AFB('b'))}, outf)
                yaml.dump({'AFBc': rnd(c.AFB('c'))}, outf)
                yaml.dump({'GammaW': rnd(c.GammaW())}, outf)
                if lfu:
                    yaml.dump({'BrWhad': rnd(c.BrWhad())}, outf)
                else:
                    yaml.dump({'BrWe': rnd(c.BrWe())}, outf)
                    yaml.dump({'BrWmu': rnd(c.BrWmu())}, outf)
                    yaml.dump({'BrWtau': rnd(c.BrWtau())}, outf)
                if MWscheme:
                    yaml.dump({'Deltaalpha': rnd(c.Deltaalpha())}, outf)
                else:
                    yaml.dump({'MW': rnd(c.MW())}, outf)
                if not lfu:
                    yaml.dump({'RWZmue': rnd(c.RWZmue())}, outf)
                    yaml.dump({'RZmue': rnd(c.RZmue())}, outf)
                    yaml.dump({'RWmue': rnd(c.RWmue())}, outf)
                    yaml.dump({'RWtaumu': rnd(c.RWtaumu())}, outf)
                    yaml.dump({'RWtaue': rnd(c.RWtaue())}, outf)
                if lfu:
                    yaml.dump({'sin2thetaleff': rnd(c.sin2thetaleff())}, outf)
                else:
                    yaml.dump({'sin2thetaeeff': rnd(
                        c.sin2thetaleff('e-'))}, outf)
                    yaml.dump({'sin2thetamueff': rnd(
                        c.sin2thetaleff('mu-'))}, outf)
                    yaml.dump({'sin2thetataueff': rnd(
                        c.sin2thetaleff('ta-'))}, outf)
                    yaml.dump({'Rsin2thetamue': rnd(c.Rsin2thetamue())}, outf)
