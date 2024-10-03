import sys,os
OUTPATH='tmpout/testSMEFT'
path=os.path.join(os.path.dirname(__file__), '../../')

sys.path.append(path)
import SMEFTlikelihood

def printAndWrite(f,*text):
    print(*text)
    f.write(" ".join([str(x) for x in text])+'\n')

with open(OUTPATH,'w') as f:
    smeftEWPD=SMEFTlikelihood.EWPDlikelihood('MW', 'top',datapath=path+'/data/pole_observables_nolfu.yml',parapath=path+'/data/pole_observables_nolfu_SMEFTpara')
    printAndWrite(f,'measured',smeftEWPD.measured())
    printAndWrite(f,'predicted',smeftEWPD.predicted())
    printAndWrite(f,'error',smeftEWPD.error())
    printAndWrite(f,'correlation',smeftEWPD.correlation())

with open(OUTPATH+'_olddata','w') as f:
    smeftEWPD=SMEFTlikelihood.EWPDlikelihood('MW', 'top',datapath='data/pole_observables_nolfu.yml',parapath=path+'/data/pole_observables_nolfu_SMEFTpara')
    printAndWrite(f,'measured',smeftEWPD.measured())
    printAndWrite(f,'predicted',smeftEWPD.predicted())
    printAndWrite(f,'error',smeftEWPD.error())
    printAndWrite(f,'correlation',smeftEWPD.correlation())

with open(OUTPATH+'_oldpara','w') as f:
    smeftEWPD=SMEFTlikelihood.EWPDlikelihood('MW', 'top',datapath=path+'/data/pole_observables_nolfu.yml',parapath='data/pole_observables_nolfu_SMEFTpara')
    printAndWrite(f,'measured',smeftEWPD.measured())
    printAndWrite(f,'predicted',smeftEWPD.predicted())
    printAndWrite(f,'error',smeftEWPD.error())
    printAndWrite(f,'correlation',smeftEWPD.correlation())

with open(OUTPATH+'_oldall','w') as f:
    smeftEWPD=SMEFTlikelihood.EWPDlikelihood('MW', 'top',datapath='data/pole_observables_nolfu.yml',parapath='data/pole_observables_nolfu_SMEFTpara')
    printAndWrite(f,'measured',smeftEWPD.measured())
    printAndWrite(f,'predicted',smeftEWPD.predicted())
    printAndWrite(f,'error',smeftEWPD.error())
    printAndWrite(f,'correlation',smeftEWPD.correlation())
