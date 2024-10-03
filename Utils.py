import os
import yaml


def correlation_dict(names, correlation):
    # symmetrize
    for n in names:
        for m in names:
            if n in correlation and m in correlation[n]:
                if not m in correlation:
                    correlation[m] = {}
                else:
                    if n in correlation[m] and abs(correlation[m][n]-correlation[n][m]) > 1e-9:
                        raise Exception('Two incompatible values for correlation of {} and {} ({} and {})'.format(
                            n, m, correlation[m][n], correlation[n][m]))
                correlation[m][n] = correlation[n][m]
    # add 1 for diagonal and 0 if no correlation
    for n in names:
        for m in names:
            if m == n:
                if not n in correlation:
                    correlation[n] = {}
                correlation[n][n] = 1.
            elif not n in correlation or not m in correlation[n]:
                if not n in correlation:
                    correlation[n] = {}
                correlation[n][m] = 0.
    return correlation


def load_data(fname):
    if not os.path.exists(fname):
        raise Exception(str(fname)+' not found')
    with open(fname, 'r') as f:
        data = yaml.safe_load(f)
    if not 'measured' in data:
        raise Exception('No "measured" in '+str(fname))
    else:
        measured = data['measured']
    if not 'correlation' in data:
        logging.warning('No correlation information in '+str(name))
        correlation = {}
    else:
        correlation = data['correlation']
    names = measured.keys()
    dcorr = correlation_dict(names, correlation)
    covariance = {}
    central = {}
    for n in names:
        covariance[n] = {}
        central[n] = float(measured[n]['central'])
        for m in names:
            covariance[n][m] = float(
                dcorr[n][m]*measured[n]['error']*measured[m]['error'])
    observable = {}
    for n in names:
        if 'observable' in measured[n]:
            observable[n] = measured[n]['observable']
        else:
            observable[n] = n
    prediction = {}
    for n in names:
        if 'prediction' in measured[n]:
            prediction[n] = measured[n]['prediction']
        else:
            prediction[n] = None
    return list(names), central, covariance, observable, prediction


def load_para(parapath, key):
    if not os.path.exists(parapath):
        raise Exception("Did not find parametrization at {}".format(parapath))
    with open(parapath) as inf:
        para = yaml.safe_load(inf)
        out = {}
        for o in para:
            if key in para[o]:
                out[o] = para[o][key]
            else:
                out[o] = {}
        return out


def toList(names, d):
    l = []
    for n in names:
        if n in d:
            l.append(d[n])
        else:
            l.append(0.)
    return l


def toMatrix(names1, names2, d):
    l = []
    for n in names1:
        l.append(toList(names2, d[n]))
    return l


def printTable(namesx, namesy, matrix, nround, width=12, pad=False):
    logging.info(
        '\t'.join([' '*width]+[(' ' if pad else '')+n[:width] for n in namesx]))
    for m in namesy:
        logging.info('\t'.join([m[:width]+' '*max(width-len(m), 0)]+[(' ' if matrix[n]
                     [m] >= 0 and pad else '')+str(round(matrix[n][m], nround)) for n in namesx]))


def printMatrix(names, matrix, nround, width=12, pad=False):
    printTable(names, names, matrix, nround, width, pad)
