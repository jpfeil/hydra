import logging
import errno
import os
import re

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


def get_genesets(gmt, gmt_regex=None):
    """
    Formats the paths to the gene set files

    :param gmt: Path to gene set database in GMT format
    :return: dictionary of gene sets
    :rtype: defaultdict(set)
    """
    if not os.path.exists(gmt):
        raise ValueError('Cannot locate GMT file: %s' % gmt)

    genesets = {}

    logging.info('GMT Regex: %s' % gmt_regex)

    if gmt_regex:
        regex = re.compile(gmt_regex)

    with open(gmt) as f:
        for line in f:
            fields = line.strip().split('\t')
            name = fields[0]
            # Skip gene sets that do not match regex
            if gmt_regex and not regex.search(name):
                continue
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


def find_aliases(gss, mapper, index):
    """
    Removes genes that do not overlap with matrix. Will
    try to rescue genes by mapping them to an alias that
    overlaps.

    :param gss (dict): Gene set database; {names: genes}
    :param mapper (pandas.DataFrame): DataFrame mapping genes to aliases
    :param index (pandas.DataFrame.index): Expression DataFrame gene index
    :return: Dictionary of gene sets including overlapping gene aliases
    :rtype: dict
    """
    logger = logging.getLogger('root')
    for gs, genes in gss.items():
        filtered = []
        for gene in genes:
            if gene in index:
                filtered.append(gene)

            else:
                aliases = mapper.loc[mapper['gene'] == gene, 'aliases']
                if len(aliases) > 0:
                    for alias in aliases:
                        matches = [x for x in alias.split('|') if x in index]
                        if len(matches) > 0:
                            for match in matches:
                                logger.info("Replacing %s in %s with %s" % (gene,
                                                                             gs,
                                                                             match))
                                filtered.append(match)
        # Remove duplicates
        gss[gs] = list(set(filtered))
    return gss