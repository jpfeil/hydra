import logging
import errno
import os

from collections import defaultdict



def mkdir_p(path):
    """
    https://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python

    :param path:
    :return:
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def get_genesets(dirs, src):
    """
    Formats the paths to the gene set files

    :param list dirs: Gene set directories
    :return: Path to gene set files
    :rtype: list
    """
    pths = []
    gs_map = {}
    for d in dirs:
        logging.info("Pulling %s gene sets:" % d)
        gs_dir = os.path.join(src, 'gene-sets', d)
        gss = os.listdir(gs_dir)
        for s in gss:
            logging.info("\t%s" % s)
            gs_pth = os.path.join(gs_dir, s)
            pths.append(gs_pth)
            gs_map[s] =  d

    return pths, gs_map


def get_test_genesets(src):
    d = os.path.join(src, 'test', 'gene-sets')
    sets = os.listdir(d)
    logging.info("Available Gene Sets:")

    for s in sets:
        logging.info(s)

    pths = []
    gs_map  = {}
    for s in sets:
        gs_map[s] = 'test'
        pth =  os.path.join(d, s)
        pths.append(pth)

    return pths, gs_map


def read_genesets(sets):
    gs = defaultdict(set)

    for s in sets:
        name = os.path.basename(s)
        with open(s) as f:
            for line in f:
                gene = line.strip()
                gs[name].add(gene)
    return gs
