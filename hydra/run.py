#!/usr/bin/env python2.7
import argparse
import bnpy
import os
import logging
import multiprocessing
import numpy as np
import pandas as pd

from collections import defaultdict

from library.utils import mkdir_p, parallel_fit

src = os.path.dirname(os.path.abspath(__file__))


def get_genesets():
    d = os.path.join(src, 'gene-sets')
    sets = os.listdir(d)
    logging.info("Available Gene Sets:")
    for s in sets:
        logging.info(s)
    sets = [os.path.join(d, s) for s in sets]
    return sets


def read_genesets(sets):
    gs = defaultdict(set)

    for s in sets:
        name = os.path.basename(s)
        with open(s) as f:
            for line in f:
                gene = line.strip()
                gs[name].add(gene)
    return gs


def is_multimodal(name, expression):
    model = parallel_fit(name, expression)
    probs = model.allocModel.get_active_comp_probs()

    if len(probs) > 1:
        return name, True

    else:
        return name, False


def filter_geneset(lst, matrx, CPU=1):
    """

    :param lst: List of genes
    :param matrx: Expression dataframe
    :return:
    """
    pool = multiprocessing.Pool(processes=CPU)

    results = []
    for gene in lst:
        data = matrx.loc[gene, :].values
        data = data - np.mean(data)
        data = data.reshape(len(data), 1)
        res = pool.apply_async(is_multimodal, args=(gene, data,))
        results.append(res)

    output = [x.get() for x in results]
    return [x[0] for x in output if x[1] is True]



def main():
    """
    Fits one, two, and three component mixture models and returns fit information in output dir.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)

    parser.add_argument('-e', '--expression',
                        help='Gene symbol by sample matrix')

    parser.add_argument('--CPU',
                        dest='CPU',
                        type=int,
                        default=1)

    parser.add_argument('--output-dir',
                        dest='output_dir',
                        default='hydra-out')

    args = parser.parse_args()

    mkdir_p(args.output_dir)

    logging.basicConfig(filename=os.path.join(args.output_dir, 'hydra.log'),
                        level=logging.INFO)

    logging.info("Started Hydra...")

    sets = get_genesets()
    genesets = read_genesets(sets)

    matrx = pd.read_csv(args.expression,
                        sep='\t',
                        index_col=0)

    # Filter genes that overlap
    for gs, genes in genesets.items():
        genesets[gs] = [x for x in genes if x in matrx.index]

    filtered_genesets = defaultdict(set)
    for gs, genes in genesets.items():

        start = len(genes)
        filtered_genesets[gs] = filter_geneset(list(genes), matrx, CPU=args.CPU)
        end = len(genesets[gs])

        logging.info("Filtering: {gs} went from {x} to {y} genes".format(gs=gs,
                                                                         x=start,
                                                                         y=end))

        center = matrx.loc[filtered_genesets[gs]].apply(lambda x: x - x.mean())

        hmodel = parallel_fit(gs, center.T.values)

        pth = os.path.join(args.output_dir, gs)
        mkdir_p(pth)

        bnpy.ioutil.ModelWriter.save_model(hmodel, pth, prefix=gs)

if __name__ == '__main__':
    main()