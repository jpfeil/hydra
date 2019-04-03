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


def get_genesets(gmt):
    """
    Formats the paths to the gene set files

    :param gmt: Path to gene set database in GMT format
    :return: dictionary of gene sets
    :rtype: defaultdict(set)
    """
    genesets = {}

    if not os.path.exists(gmt):
        raise ValueError('Cannot locate GMT file: %s' % gmt)

    with open(gmt) as f:
        for line in f:
            fields = line.strip().split('\t')
            name = fields[0]
            logging.info("Loading gene set: %s" % name)
            genes = fields[2:]
            genesets[name] = set(genes)
    return genesets


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
