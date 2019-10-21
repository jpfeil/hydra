#!/usr/bin/env python2.7
import bnpy
import logging
import pandas as pd
import numpy as np
import os
from library.fit import get_assignments

def predict(expression, args):
    """
    Takes a prefit model and returns the assignments

    :param expression:
    :param model_path:
    :return:
    """
    logger = logging.getLogger('root')
    logger.info('Loading model:\n%s' % args.model_path)

    name = os.path.basename(args.model_path)
    model = bnpy.ioutil.ModelReader.load_model_at_prefix(args.model_path,
                                                         prefix=name)

    pth = os.path.join(args.model_path, 'training-data.tsv')
    train = pd.read_csv(pth, sep='\t', index_col=0)

    data = expression.reindex(train.index).dropna()
    logger.info('Centering data.')
    data = data.sub(train.mean(axis=1), axis=0)
    if data.ndim == 1:
        data = data.values.reshape(1, data.shape[0])
    else:
        data = data.T.values

    xdata = bnpy.data.XData(data)

    logger.info('Applying model to data:\n%s' % args.expression)
    unclass = 1 - np.sum(model.allocModel.get_active_comp_probs())
    LP = model.calc_local_params(xdata)
    asnmts = pd.DataFrame(index=expression.columns,
                          columns=[name])
    for sample, row in zip(asnmts.index.values, range(LP['resp'].shape[0])):
        _max = np.max(LP['resp'][row, :])
        if _max < unclass:
            logging.info("WARNING: At least one sample was not classified!!!")
            asnmts.loc[sample, name] = -1
        else:
            _arg = np.argmax(LP['resp'][row, :])
            asnmts.loc[sample, name] = _arg

    pth = os.path.join(args.output_dir, 'PREDICT', name, 'assignments.tsv')
    asnmts.to_csv(pth, sep='\t')
    return asnmts
