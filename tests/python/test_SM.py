import sys,os
OUTDIR='tmpout/testSM'

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
import SMcalculator

def printAndWrite(f,*text):
    print(*text)
    f.write(" ".join([str(x) for x in text])+'\n')

with open(OUTDIR,'w') as f:
    sm_MWscheme=SMcalculator.EWPOcalculator(MH=125.,mt=173.,alphas=0.118,MZ=91.1875,MW=80.37)
    printAndWrite(f,'AFBb:',sm_MWscheme.AFBb())
    sm_MWscheme.update(MW=80.3)
    printAndWrite(f,'AFBb(MW=80.3):',sm_MWscheme.AFBb())
    sm_MWscheme.set('MW',80.4) # same as update
    printAndWrite(f,'AFBb(MW=80.4):',sm_MWscheme.AFBb())
    sm_MWscheme.reset()
    printAndWrite(f,'AFBb:',sm_MWscheme.get('AFBb'))
    printAndWrite(f,'dAFBb/dMW:',sm_MWscheme.derivative('AFBb','MW'))
    printAndWrite(f,'All observables:', sm_MWscheme.getall())
    sm_alphascheme=SMcalculator.EWPOcalculator(MH=125.25,mt=172.69,alphas=0.118,MZ=91.1875,Deltaalpha=0.0590)
    printAndWrite(f,'AFBb(alpha scheme):',sm_alphascheme.get('AFBb'))
    sm_alphascheme2=SMcalculator.EWPOcalculator(scheme='alpha',input_dict={'MH':125.25,'mt':172.69,'alphas':0.118,'MZ':91.1875,'Deltaalpha':sm_MWscheme.Deltaalpha()}) # alternative constructor using dict
    printAndWrite(f,'AFBb(alpha(MW) scheme):',sm_alphascheme2.get('AFBb'))
    printAndWrite(f,'(dAFBb/dDeltaalpha)/(dMW/ddDeltaalpha):',sm_alphascheme2.derivative('AFBb','Deltaalpha')/sm_alphascheme2.derivative('MW','Deltaalpha'))
