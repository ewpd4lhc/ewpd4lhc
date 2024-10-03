import para_factory
models = ['SMEFTsim_U35_alphaScheme_UFO']
mg5dir = 'MG5_aMC_v3_5_6/'
runner = para_factory.MGrunner(mg5dir=mg5dir, tmpdir='/tmp/MGtmp/',
                               outdir='testrawparas', SMparams={'aEW': 0.0072973525205055605, 'aS': 0.1181})
runner.run_all(doQuad=False, models=models)
shifts = para_factory.MWalphaShifts(
    modelpath=mg5dir+'/models/', outdir='testrawparas')
shifts.run(models=models)
calc = para_factory.CalcParaWrapper(paradir='testrawparas', outdir='testparas')
calc.run(models=models)
