import sys,os,yaml
OUTDIR='tmpout/testSchemes'

mainpath=os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(mainpath)

import SMcalculator

sm_alpha=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.alpha)
sm_MW=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.MW)
sm_alphaMW=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.alphaMW)
sm_sin2theta=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.sin2theta)
sm_alphasin2theta=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.alphasin2theta)

sm_MW0=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.MW,input_dict=sm_alpha.getall())
sm_alphaMW0=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.alphaMW,input_dict=sm_alpha.getall())
sm_sin2theta0=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.sin2theta,input_dict=sm_alpha.getall())
sm_alphasin2theta0=SMcalculator.EWPOcalculator(scheme=SMcalculator.INPUTSCHEME.alphasin2theta,input_dict=sm_alpha.getall())


print('Prediction with alpha scheme values as input')
for scheme in '', 'alpha', 'MW', 'alphaMW', 'sin2theta', 'alphasin2theta':
    print(f'{scheme: <20}',end='')
print()
for o in SMcalculator.EWPOS:
    print(f'{o: <20}',end='')
    for sm in sm_alpha, sm_MW0, sm_alphaMW0, sm_sin2theta0, sm_alphasin2theta0:
        print('{0: <20}'.format(round(sm.get(o),15)),end='')
    print()

print('Prediction with PDG values for inputs')
for scheme in '', 'alpha', 'MW', 'alphaMW', 'sin2theta', 'alphasin2theta':
    print(f'{scheme: <20}',end='')
print()
for o in SMcalculator.EWPOS:
    print(f'{o: <20}',end='')
    for sm in sm_alpha, sm_MW, sm_alphaMW, sm_sin2theta, sm_alphasin2theta:
        print('{0: <20}'.format(round(sm.get(o),15)),end='')
    print()

print('Now some numerical checks')
for sm1 in sm_alpha, sm_MW, sm_alphaMW, sm_sin2theta, sm_alphasin2theta:
    for sm2 in sm_alpha, sm_MW, sm_alphaMW, sm_sin2theta, sm_alphasin2theta:
        sm1.reset()
        sm2.reset()
        if sm1==sm2:
            continue
        sm2.update(input_dict=sm1.getall())
        for o in SMcalculator.EWPOS:
            #print(o,sm1.scheme(),sm2.scheme(),sm1.get(o),sm2.get(o))
            assert abs(sm1.get(o)/sm2.get(o)-1)<1e-12

datapath=os.path.join(mainpath,'data/pole_observables_lfu.yml')
data=yaml.safe_load(open(datapath))

for sm in sm_alpha, sm_MW, sm_alphaMW, sm_sin2theta, sm_alphasin2theta:
    for inp in sm.inputs():
        outputs = [o for o in sm.outputs() if o in data['measured']]
        if inp=='sin2thetaleff':
            x=0.23150
            err=0.00020
        else:
            x=data['measured'][inp]['central']
            err=data['measured'][inp]['error']
        sm.reset()
        pred_default={}
        pred_up={}
        pred_down={}
        derivative={}
        sigmas=3
        for o in outputs:
            pred_default[o]=sm.get(o)
            derivative[o]=sm.derivative(o,inp)
        sm.update(input_dict={inp:x+err*sigmas})
        for o in outputs:
            pred_up[o]=sm.get(o)
        sm.update(input_dict={inp:x-err*sigmas})
        for o in outputs:
            pred_down[o]=sm.get(o)
        for o in outputs:
            linear = sigmas*derivative[o]*err/pred_default[o]
            up = (pred_up[o]-pred_default[o])/pred_default[o]
            down = (pred_down[o]-pred_default[o])/pred_default[o]
            devup=up/linear-1 if up!=0 else 0
            devdown = -down/linear-1 if down!=0 else 0
            dev=(up-down)/2/linear-1 if linear!=0 else 0
            if linear==0:
                assert abs(up)<1e-12 and abs(down)<1e-12
            # check, if within expected 1% if symmetrized
            assert abs(dev)<0.01 or abs(up)<1e-12
            cutoff=0.02
            if not (abs(devup)<cutoff or abs(up)<1e-12):
                print('More than {} deviation for {} sigma up variation'.format(cutoff,sigmas),sm.scheme(),'d'+o+'_d'+inp,devup)
            if not (abs(devdown)<cutoff or abs(down)<1e-12):
                print('More than {} deviation for {} sigma down variation'.format(cutoff,sigmas),sm.scheme(),'d'+o+'_d'+inp,devdown)


                

            
            
            
        
