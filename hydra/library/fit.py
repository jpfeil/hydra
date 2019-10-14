#!/usr/bin/env python2.7
import os
import re
import bnpy
import tempfile
import shutil
import shlex
import sys
import numpy as np
import logging

from subprocess import Popen, PIPE
from threading import Timer

from library.notebook import create_notebook


# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def run(cmd, timeout_sec=900):
    """
    https://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
    :param cmd:
    :param timeout_sec:
    :return:
    """
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
    timer = Timer(timeout_sec, proc.kill)

    try:
        timer.start()
        stdout, stderr = proc.communicate()

    finally:
        timer.cancel()

    return stdout, stderr


converged_regex = re.compile("... done. converged.")


def subprocess_fit(name,
                   dataset,
                   gamma=1.0,
                   sF=1.0,
                   K=3,
                   nLap=1000,
                   bstart=0,
                   mstart=2,
                   dstart=2,
                   save_output=False,
                   timeout_sec=900):
    """

    :param name:
    :param dataset: bnpy.data.Xdata object
    :param gamma:
    :param sF:
    :param K:
    :param save_output:
    :param timeout_sec:
    :return:
    """
    logger = logging.getLogger('root')

    name = re.sub(r'[^\w\d-]', '_', name)

    workdir = tempfile.mkdtemp(prefix="%s_" % name)
    output_dir = 'K={K}-gamma={G}-ECovMat={Cov}-moves=birth,merge,delete,shuffle/'.format(K=K,
                                                                                          G=gamma,
                                                                                          Cov=sF)
    output_path = os.path.join(workdir, output_dir)

    csv_file_path = os.path.join(workdir, "%s.csv" % name)
    assert dataset.dim > 0, 'Dataset has not dimension'
    assert dataset.get_size() > 2, 'Dataset needs more than two observations!'
    logger.debug("Gene %s\nSize %d\n" % (name, dataset.get_size()))
    dataset.to_csv(csv_file_path)

    cmd = """
          python -m bnpy.Run {data} 
                    DPMixtureModel 
                    Gauss 
                    memoVB 
                    --K {K}
                    --nLap {nLap}
                    --nTask 1
                    --nBatch 1
                    --gamma0 {G}
                    --sF {sF}
                    --b_startLap {bstart}
                    --m_startLap {mstart}
                    --d_startLap {dstart}
                    --ECovMat eye
                    --moves birth,merge,delete,shuffle
                    --output_path {output}
          """.format(data=csv_file_path,
                     K=K,
                     nLap=nLap,
                     G=gamma,
                     sF=sF,
                     bstart=bstart,
                     mstart=mstart,
                     dstart=dstart,
                     output=output_path)

    stdout, stderr = run(shlex.split(cmd), timeout_sec=timeout_sec)

    logger.debug(cmd)
    logger.debug(stdout)
    logger.debug(stderr)

    converged = False
    m = converged_regex.search(stdout)
    if m:
        converged = True

    else:
        logger.debug(stdout)
        logger.debug(stderr)

    try:
        hmodel = bnpy.ioutil.ModelReader.load_model_at_prefix(os.path.join(output_path, '1'),
                                                              prefix='Best')
    except IOError:
        print(stdout)
        raise ValueError("%s model not found!" % name)

    params = "Gamma: {G}\nK: {K}\nsF: {sF}\nbStart: {b}\n" \
             "mStart: {m}\ndStart: {d}\nnLap: {n}".format(G=gamma,
                                                          K=K,
                                                          sF=sF,
                                                          b=bstart,
                                                          m=mstart,
                                                          d=dstart,
                                                          n=nLap)
    if not save_output:
        shutil.rmtree(workdir)
    else:
        print("Output:\n%s" % workdir)
    return hmodel, converged, params, stdout, stderr


def get_assignments(model, data):
    """
    Takes model and data and classifies samples

    Will label samples with NaN if they do not
    fit in any of the model components

    :param model:
    :param data:
    :return:
    """
    unclass = 1 - np.sum(model.allocModel.get_active_comp_probs())
    # Get the sample assignments
    LP = model.calc_local_params(data)
    asnmts = []
    for row in range(LP['resp'].shape[0]):
        _max = np.max(LP['resp'][row, :])
        if _max < unclass:
            asnmts.append(np.nan)

        else:
            _arg = np.argmax(LP['resp'][row, :])
            asnmts.append(_arg)

    return asnmts


def apply_multivariate_model(input, args, output, name='MultivariateModel'):
    """

    :param input:
    :param args:
    :param output:
    :param name:
    :return:
    """
    logger = logging.getLogger('root')
    logger.info("Centering input to multivariate clustering.")

    name = re.sub(r'[^\w\d-]', '_', name)

    center = input.apply(lambda x: x - x.mean(), axis=1)

    # Need to take the transpose
    # Samples x Genes
    data = center.T.values

    # Create dataset object for inference
    dataset = bnpy.data.XData(data)

    # Set the prior for creating a new cluster
    gamma = args.gamma

    # Start with a standard identity matrix
    sF = args.sF

    # Starting with 5 cluster because starting with
    # 1 cluster biases the fit towards not finding clusters.
    K = args.K

    nLap = args.num_laps

    logger.info("Multivariate Model Params:\n"
                 "gamma: %.2f\n"
                 "sF: %.2f\n"
                 "K: %d\n"
                 "nLaps: %d" % (gamma, sF, K, nLap))

    # Fit multivariate model
    hmodel, converged, params, stdout, stderr = subprocess_fit(name,
                                                               dataset,
                                                               gamma,
                                                               sF,
                                                               K,
                                                               nLap=nLap,
                                                               timeout_sec=args.max_fit_time)

    if converged is False:
        logging.info("WARNING: Multivariate model did not converge!")
        pth = os.path.join(output, 'NOT_CONVERGED')
        with open(pth, 'w') as f:
            f.write("WARNING: Multivariate model did not converge!")

    asnmts = get_assignments(hmodel, dataset)

    pth = os.path.join(output, 'assignments.tsv')
    with open(pth, 'w') as f:
        for sample, asnmt in zip(center.columns, asnmts):
            f.write('{sample}\t{assignment}\n'.format(sample=sample,
                                                      assignment=asnmt))

    # Save model
    bnpy.ioutil.ModelWriter.save_model(hmodel,
                                       output,
                                       prefix=name)

    try:
        bnpy.ioutil.ModelReader.load_model_at_prefix(output,
                                                     prefix=name)

    except IOError:
        logging.error("Error in writing model!")
        raise

    create_notebook(name, output)

    with open(os.path.join(output, 'PARAMS'), 'w') as f:
        f.write(params)

    with open(os.path.join(output, 'STDOUT'), 'w') as f:
        f.write(stdout)
