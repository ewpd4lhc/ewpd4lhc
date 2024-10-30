# EWPD for LHC

This code creates a SMEFT likelihood for EWPD (Z pole and W pole observables), taylored to the needs of the LHC collaborations.

* Supports not only the {alpha,MZ,GF} input scheme but also {MW,MZ,GF} (preferred at LHC), {alpha,MZ,MW},  {sin2theta,MZ,GF}, {alpha,MZ,sin2theta} 
* Dynamically calculates state-of-the art SM predictions from input parameters
* Includes of parametric and theory uncertainties as well as their correlated effect
* Loads flexible EFT parametrizations from text files
* Parametrizations with dimension-eight as well as NLO perturbative contributions are included (up to 5 input parameter schemes)
* Baseline parametrization and notation is based on SMEFTsim, with other parametrizations validated against it
* SMEFT symmetry assumptions used in LHC interpretations are available (e.g. "top" symmetries with U(2)^3 symmetry in the quark sector)
* Output as yaml file and Roofit workspace

## Usage

Clone with
```
git clone https://github.com/ewpd4lhc/ewpd4lhc.git
cd ewpd4lhc
```

Requires `python3` with the  `numpy` and `yaml` modules. `ROOT` for Roofit output.
Runs without further preparation on CERN lxplus.

### Configuration
The main file for configuration is [config/ewpd4lhc.cfg](config/ewpd4lhc.cfg).
Here one can choose:

* The list of variables (many precision observables are implemented, both with and without lepton flavour universality)
* The input scheme (possible to choose between {MW,MZ,Gmu} and {alpha,MZ,Gmu})
* Treatment of theoretical and parametric uncertainties (either as part of the covariance of the Gaussian or as nuisance parameters)
* The SMEFT symmetry
* Possibly a subset of coefficients
* The source for dimension-six linear, dimension-six squared, and dimension-eight contributions 

Measurement data is stored as yaml files in [data/](data).
In the same folder there are several SMEFT parametrizations to choose from.

### Running the code
A yaml file describing the likelihood can be created like so:

```
./ewpd4lhc.py
```

where by default [config/ewpd4lhc.cfg](config/ewpd4lhc.cfg) is taken as input and the SMEFT likelihood is stored in a textfile name `ewpd_out.yml`.
Different configs and output names can be specified with command line arguments.
All arguments can be listed with:

```
./ewpd4lhc.py --help
```


For example, a ROOT workspace can be created either directly:

```
./ewpd4lhc.py --root_output ROOTFILENAME.root
```

Or in a seperate step from the output textfile (allowing the user to modify or build more complex likelihoods) with:

```
./ewpd4lhc.py --output YAMLFILE.yml
./yaml2root.py --input YAMLFILE.yml --output ROOTFILENAME.root
```

A simple example for fitting the ROOT file is implemented.
Workspace contents are printed with:

```
./ROOTfit/fit.py --input ROOTFILENAME.root
```

A multi-parameter fit is performed, e.g., with:

```
./ROOTfit/fit.py --input ROOTFILENAME.root --pois=cHu,cHd,cHj3,cHj1
```

It is possible to specify `--poi=all` but this will usually not converge without extra constraints as the EWPD likelihood is degenerate.
1D and 2D scans of the likelihood can be performed with the following command, with scan points being stored in a textfile

```
./ROOTfit/fit.py --input ROOTFILENAME.root --pois=cHu --scan=-0.1:0.1 --outfolder=1Dscan
./ROOTfit/fit.py --input ROOTFILENAME.root --pois=cHu,cHd --scan=-0.1:0.1,0.2:0.2 --outfolder=2Dscan
```

One can additionally profile parameters in the scan and plot:

```
./ROOTfit/fit.py --input ROOTFILENAME.root --pois=cHu --scan=-0.5:0.5 --float=cHd,cHj1,cHj3 --outfolder=1DscanProfiled --plot
./ROOTfit/fit.py --input ROOTFILENAME.root --pois=cHu,cHd --scan=-0.5:0.5,-1:1 --float=cHj1,cHj3 --outfolder=2DscanProfiled --plot
```

It is also possible to use the low level classes `EWPOcalculator` and `SMEFTlikelihood`. For example, SM predictions can be calculated like so:

```import SMcalculator
sm=SMcalculator.EWPOcalculator(MH=125.25,mt=172.69,alphas=0.118,MZ=91.1875,MW=80.377)
print('AFBb:',sm.AFBb())
sm.update(MW=80.3)
print('AFBb(MW=80.3):',sm.AFBb())
sm.reset()
print('Also AFBb:',sm.get('AFBb'))
print('dAFBb/dMW:',sm.derivative('AFBb','MW'))
print('All observables:', sm.getall())
```

## More detailed description

### Experimental data

Z pole observables from LEP and SLD [^1] are included, both assuming LFU ([data/pole_observables_lfu.yml](data/pole_observables_lfu.yml)) and not assuming LFU ([data/pole_observables_nolfu.yml](data/pole_observables_nolfu.yml)). The correct set corresponding to the SMEFT symmetry assumption is used automatically.

For these observables, the values from Ref [^2], extrated using an updated luminosity calculation are used.

PDG [^3] values for Higgs mass, top mass, as well as W mass and width are taken.

For the top mass a theoretical uncertainty of 0.5 GeV is added in quadrature to the experimental uncertainty.

The value of Deltaalpha (where Deltaalpha is defined as Deltaalpha = alpha(MZ)/alpha(0)-1, the shift in the electromagnetic coupling alpha due to the running from q2=0 to q2=MZ) is calculated as the sum of hadronic contributions [^4] (extracted using experimental data) and the theoretical calculation for the leptonic contribution [^7].
The experimental correlation of Refs [^1], [^2] is taken into account and noted down in the yaml files.

For alphaS the value is taken from lattic calculation [^5], according to recommendation in Ref. [^6]

The full list of observables currently implemented can be understood from the `data` files.

### SM prediction

SM predictions are calculated based on the SM input variables MZ, MH, mt, alphas, and either MW (MW scheme) or Deltaalpha (alpha scheme), which are extracted from the tml files in [data/](data/).
The interpolation formulas from Refs. [^8]-[^11] are used. Theoretical uncertainties arising in the calculation of MW, GammaZ, Rb, Rc, Rl, sigmahad are taken from the relevant publications and treated as uncorrelated. A correlated uncertainity is implemented for $sin^{2}\theta_{l}^{eff}$, using the value of Ref. [^8].

In the MW scheme, SM predictions are also derived using Refs [^8][^9][^11], but the formula of Ref [^10] is inverted to calculated Deltaalpha(MW), following the methodology of Ref [^12].

The formulas are also used to derive the linear dependence of all observables on the input parameters (where the derivative is evaluated using the central values of the experimental data).
Given the experimentally precision of the five input quantities, the input parameter dependence beyond linear order is negligible.

### Model for uncertainties

The likelihood L is modelled as multivariate Gaussian, with
$$log L = \Delta x_{i} C^{-1}_{ij} \Delta x_{i}$$,
where $\Delta x_i$ is the difference between measured and predicted value for the observable $i$.

The covariance $C$ is based on the experimental uncertainty in the measured value and potentially additional contributions.

For the inclusion of theoretical uncertainties (e.g. due to missing higher orders) and parametric uncertainties (due to uncertainties in the five input parameters) two options are provided.

The first method is their inclusion as nuisance parameters. In that case the covariance also encompasses the five observables serving as SM inputs.
The SM inputs are parameters of the model that modify predictions (and thus $\Delta x_i$) and are constrained by the multivariate Gaussian. Only the linear dependence on input parameters is included.
Theoretical uncertainties are also modelled with nuisance parameters that shift predictions, each constrained by a Gaussian distributions.

The second method is to add them to the experimental covariance.
In that case the model is simpler, the covariance does not encompas SM inputs and the no nuisance parameters are introduced.
The parametric (or similar the theory) covariance $ C^{param}_{ij} $ of two observables $O_i$ and $O_j$ is calcuated as
$C^{param}_{ij} = \sum_{k} C^{param,k}_{ij}$
where
$C^{param,k}_{ij} = (dO_i/dO_k)\times\sigma_k\times(dO_j/dO_jk)\times\sigma_k$, 
$\sigma_k$ is the uncertainty in the input parameter $k$, and $dO_i/dO_k$ the analytically calculated partial derivate.
The total covariance is the sum of experimental, parametric, and theory uncertainty covariance.

The calculation of the above covariance has been validated as follows. Based on the experimental covariance 10 M pseudo data sets have been created. For each the difference between SM prediction (based on the five inputs) and pseudo observation is calculated. The covariance of the ensemble and the analytical calculation described above agree almost perfectly, with permille-level difference that can be attributed to non-linear effects.

In particular the relatively large uncertainty on MW scheme introduces, in the MW scheme, a large correlated effect in the difference between SM prediction and observation (which is ultimately the relevant quantity in the SMEFT interpretation).

#### SM fit 

To validate the inputs and the implementation of SM prediction formulas, an analytic SM fit is performed automatically when running the tool, using the configured variable list (+SM inputs).
The SM parametrizations are used in the linear approximation, which is sufficient in most cases, except for the indirect determination of the Higgs mass, where non-linear effects are important due to the low precision.


### SMEFT parametrization

It is possible to choose from three parametrizations. There is also the option to mix parametrizations, e.g. taking linear dimension-six effects from SMEFTsim and only higher-order effects from EWPD2dim8.
The parametrizations generally agree within 10% for linear dimension-six predictions (the only level at which they are comparable).

#### SMEFTsim
Parametrizations extracted from the SMEFTsim 3.0 model [^14].

Features: Both MW and alpha schemes and all symmetries, all variables. Dimension-six and dimension-six squared terms.

The parametrizations for W and Z bosons have been derived from SMEFTsim 3.0 by generating calculating polarized decay rates in MadGraph5_aMCatNLO [^15]. Shifts in MW and alpha have been extracted from the UFO model.
The advantage of this parametrization is small numerical uncertainties and the option to include "quadratic terms", which are partial 1/Lambda^4 contributions that are typically included LHC interpretations.

#### EWPD2dim8
Parametrization extracted from the tables in "EWPD to dimension eight" [^13].

Features: Both schemes and all symmetries, Z pole variables only. Dimension-six, dimension-six squared terms, and dimension-eight terms.

Parametrizations have been calculated from the tables in Ref [^13] and converted to various symmetries.
There is also the option to include quadratic terms. Compared to SMEFTsim additional terms quadratic in dimenson six Wilson coefficients are taken into account.
The only parametrizations that includes dimension eight contributions. Relatively small numerical precision due to the presentation of numbers.

#### EWPDatNLO
Parametrization taken from the supplementary material to "Electroweak and QCD corrections to Z and W pole observables in the SMEFT" [^16]

Features: Only alpha scheme and U35 symmetry, all variables. Dimension-six terms only.

The only parametrization with NLO perturbative corrections. NLO corrections modify the values of coefficients appearing at LO and introduce additional coefficients.


## References

[^1]: S. Schael et al. [ALEPH, DELPHI, L3, OPAL, SLD, LEP Electroweak Working Group, SLD Electroweak Group and SLD Heavy Flavour Group], "Precision electroweak measurements on the $Z$ resonance," Phys. Rept. 427 (2006), 257-454 [https://arxiv.org/abs/hep-ex/0509008]

[^2]: P. Janot and S. Jadach, "Improved Bhabha cross section at LEP and the number of light neutrino species," Phys. Lett. B 803 (2020), 135319 [https://arxiv.org/abs/1912.02067]

[^3]: R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01 (2022) [https://pdg.lbl.gov/]

[^4]: M. Davier, A. Hoecker, B. Malaescu and Z. Zhang, "Reevaluation of the hadronic vacuum polarisation contributions to the Standard Model predictions of the muon $g-2$ and ${\alpha (m_Z^2)}$ using newest hadronic cross-section data," Eur. Phys. J. C 77 (2017) no.12, 827 [https://arxiv.org/abs/1706.09436]

[^5]: Y. Aoki et al. [Flavour Lattice Averaging Group (FLAG)], "FLAG Review 2021", Eur. Phys. J. C 82 (2022) no.10, 869 [https://arxiv.org/abs/2111.09849]

[^6]: M. Trott, "$\alpha_s$ as an input parameter in the SMEFT", [https://arxiv.org/abs/2306.14784].

[^7]: M. Steinhauser, "Leptonic contribution to the effective electromagnetic coupling constant up to three loops," Phys. Lett. B 429 (1998), 158-161 [https://arxiv.org/abs/hep-ph/9803313]

[^8]: I. Dubovyk, A. Freitas, J. Gluza, T. Riemann and J. Usovitsch, "Electroweak pseudo-observables and Z-boson form factors at two-loop accuracy," JHEP 08 (2019), 113 [https://arxiv.org/abs/1906.08815]

[^9]: M. Awramik, M. Czakon and A. Freitas, "Electroweak two-loop corrections to the effective weak mixing angle", JHEP 11 (2006), 048 [https://arxiv.org/abs/hep-ph/0608099].

[^10]: M. Awramik, M. Czakon, A. Freitas and G. Weiglein, "Precise prediction for the W boson mass in the standard model," Phys. Rev. D 69 (2004), [https://arxiv.org/abs/hep-ph/0311148]

[^11]: G. C. Cho, K. Hagiwara, Y. Matsumoto and D. Nomura, "The MSSM confronts the precision electroweak data and the muon g-2," JHEP 11 (2011), 068 [https://arxiv.org/abs/1104.1769]

[^12]:  I. Brivio and M. Trott, "Scheming in the SMEFT... and a reparameterization invariance!," JHEP 07 (2017), 148 [https://arxiv.org/abs/1701.06424]

[^13]: T. Corbett, A. Helset, A. Martin and M. Trott, "EWPD in the SMEFT to dimension eight," JHEP 06 (2021), 076 [https://arxiv.org/abs/2102.02819]

[^14]: I. Brivio, "SMEFTsim 3.0 - a practical guide," JHEP 04 (2021), 073 [https://arxiv.org/abs/2012.11343]

[^15]: J. Alwall, R. Frederix, S. Frixione, V. Hirschi, F. Maltoni, O. Mattelaer, H. S. Shao, T. Stelzer, P. Torrielli and M. Zaro, "The automated computation of tree-level and next-to-leading order differential cross sections, and their matching to parton shower simulations," JHEP 07 (2014), 079 [https://arxiv.org/abs/1405.0301]

[^16]: S. Dawson and P. P. Giardino, "Electroweak and QCD corrections to $Z$ and $W$ pole observables in the standard model EFT," Phys. Rev. D 101 (2020) no.1, 013001 [https://arxiv.org/abs/1909.02000]

[^17]: W. Verkerke and D. P. Kirkby, "The RooFit toolkit for data modeling," eConf C0303241 (2003), MOLT007 [https://arxiv.org/abs/physics/0306116]

[^18]: L. Moneta, K. Belasco, K. S. Cranmer, S. Kreiss, A. Lazzaro, D. Piparo, G. Schott, W. Verkerke and M. Wolf, "The RooStats Project", PoS ACAT2010 (2010), 057  [https://arxiv.org/abs/1009.1003]
