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


def run(cmd, timeout_sec):
    """
    https://stackoverflow.com/questions/1191374/using-module-subprocess-with-timeout
    :param cmd:
    :param timeout_sec:
    :return:
    """
    proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE)
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
                    --b_startLap 10
                    --m_startLap 10
                    --d_startLap 20
                    --ECovMat eye
                    --moves birth,merge,delete,shuffle
                    --output_path {output}
          """.format(data=csv_file_path,
                     K=K,
                     nLap=nLap,
                     G=gamma,
                     sF=sF,
                     output=output_path)

    stdout, _ = run(cmd, timeout_sec=timeout_sec)

    converged = False
    m = converged_regex.search(stdout)
    if m:
        converged = True

    hmodel = bnpy.ioutil.ModelReader.load_model_at_prefix(os.path.join(output_path, '1'),
                                                          prefix='Best')

    if not save_output:
        shutil.rmtree(workdir)

    else:
        print("Output:\n%s" % workdir)

    return hmodel, converged
