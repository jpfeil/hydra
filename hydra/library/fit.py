import os
import re
import bnpy
import tempfile
import shutil
import shlex
import sys

from subprocess import Popen, PIPE
from threading import Timer


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
    :param dataset:
    :param gamma:
    :param sF:
    :param K:
    :param save_output:
    :param timeout_sec:
    :return:
    """

    workdir = tempfile.mkdtemp(prefix="%s_" % name)
    output_dir = 'K={K}-gamma={G}-ECovMat={Cov}-moves=birth,merge,delete,shuffle/'.format(K=K,
                                                                                          G=gamma,
                                                                                          Cov=sF)
    output_path = os.path.join(workdir, output_dir)

    csv_file_path = os.path.join(workdir, "%s.csv" % name)
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

    converged = False
    m = converged_regex.search(stdout)
    if m:
        converged = True

    hmodel = bnpy.ioutil.ModelReader.load_model_at_prefix(os.path.join(output_path, '1'),
                                                          prefix='Best')

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

    return hmodel, converged, params, stdout
