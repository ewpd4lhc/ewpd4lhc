import sys,os,yaml
OUTDIR='tmpout/'

scheme='MW'

path='../data/'

D6QUADKEY='d6quad'
D8LINKEY='d8lin'
SMKEY='sm'
D6LINKEY='d6lin'

paras=[
    ('SMEFTsim     ',f'{path}pole_observables_lfu_SMEFTpara/SMEFTsim_U35_{scheme}.yml'),
    ('EWPD2dim8    ',f'{path}pole_observables_lfu_SMEFTpara/EWPD2dim8_U35_{scheme}.yml'),
    ('EWPDatNLO LO',f'{path}pole_observables_lfu_SMEFTpara/EWPDatLO_U35_{scheme}.yml'),
    ('EWPDatNLO NLO',f'{path}pole_observables_lfu_SMEFTpara/EWPDatNLO_U35_{scheme}.yml'),]



paras = {name:yaml.safe_load(open(fn)) for name,fn in paras}


coeffsd6_set=set()
coeffsd62_set=set()
coeffsd8_set=set()
obs = set()
for p in paras:
    for o in paras[p]:
        obs.add(o)
        for c in paras[p][o][D6LINKEY]:
            coeffsd6_set.add(c)
        if D6QUADKEY in paras[p][o]:
            for c in paras[p][o][D6QUADKEY]:
                coeffsd62_set.add(c)
        if D8LINKEY in paras[p][o]:
            for c in paras[p][o][D8LINKEY]:
                coeffsd8_set.add(c)
coeffs_order='cHWB cHDD cHl1 cHl3 cHe cHq1 cHq3 cHu cHd cll1 cHbox cW cHB cHW cuB cuW cle cee clq1 clq3 clu cld cqe ceu ced cqq1 cqq3 cqu1 cqd1 cud1 cuu cdd'.split()
coeffs_order_d8='c8HWB c8HDD c8HDD2 c8Hl1 c8Hl2 c8Hl3 c8He c8Hq1 c8Hq2 c8Hq3 c8Hd dGF8'.split()


coeffsd6=[]
coeffsd62=[]
coeffsd8=[]
for c in coeffs_order:
    if c in coeffsd6_set:
        coeffsd6.append(c)
    for c2 in coeffs_order:
        cc = c+'*'+c2
        if cc in coeffsd62_set:
            coeffsd62.append(cc)
for c in coeffs_order_d8:
    if c in coeffsd8_set:
        coeffsd8.append(c)



for o in obs:
    print('='*80)
    print(o)
    print('-'*80)
    print(' '*5,end='\t')
    for p in paras:
        print(p,end='\t')
    print()    
    print('-'*80)
    for c in coeffsd6:
        print(c,end='\t')
        for p in paras:
            x = str(round(paras[p][o][D6LINKEY][c],3) if o in paras[p] and c in paras[p][o][D6LINKEY] else ' - ')+'\t'
            print(x,end='\t')
        print()

    print('-'*80)
    for c in coeffsd62:
        print(c,end='\t')
        for p in paras:
            x = str(round(paras[p][o][D6QUADKEY][c],3) if o in paras[p] and D6QUADKEY in paras[p][o] and c in paras[p][o][D6QUADKEY] else ' - ')+'\t'
            print(x,end='\t')
        print()

    print('-'*80)
    for c in coeffsd8:
        print(c,end='\t')
        for p in paras:
            x = str(round(paras[p][o][D8LINKEY][c],3) if o in paras[p] and D8LINKEY in paras[p][o] and c in paras[p][o][D8LINKEY] else ' - ')+'\t'
            print(x,end='\t')
        print()
    
    print('-'*80)
