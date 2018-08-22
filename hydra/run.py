#!/usr/bin/env python2.7
import argparse
import bnpy
import os
import logging
import multiprocessing
import numpy as np
import pandas as pd

from collections import defaultdict
from scipy.stats import kruskal

from library.utils import mkdir_p, parallel_fit
from library.notebook import create_notebook

src = os.path.dirname(os.path.abspath(__file__))


def get_genesets(dirs):
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


def get_test_genesets():
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


def distinct_covariates(data, model, alpha=0.01):
    """
    Determine if clustering separates data by components

    :param data:
    :param model:
    :return:
    """
    LP = model.calc_local_params(data)
    asnmts = LP['resp'].argmax(axis=1)

    exp_groups = [[] for _ in range(max(asnmts) + 1)]
    cov_groups = [[] for _ in range(max(asnmts) + 1)]

    exp = data.X[:, 0]
    cov = data.X[:, 1]

    for e, c, a in zip(exp, cov, asnmts):
        exp_groups[a].append(e)
        cov_groups[a].append(c)

    found_exp = False
    found_cov = False

    _, exp_p = kruskal(*exp_groups, nan_policy='raise')

    if exp_p < alpha / (max(asnmts) + 1.0):
        found_exp = True

    _, cov_p = kruskal(*cov_groups, nan_policy='raise')
    if cov_p < alpha / (max(asnmts) + 1.0):
        found_cov = True

    return found_exp and found_cov


analyzed = {}
def is_multimodal(name, data, covariate=None, min_prob_filter=None, output_dir=False):
    """
    This function determines if there is a multimodal pattern in the data

    :param name:
    :param data:
    :param covariate:
    :param min_prob_filter:
    :return:
    """
    if name in analyzed:
        return analyzed[name]

    X = bnpy.data.XData(data)
    model = parallel_fit(name, X)
    probs = model.allocModel.get_active_comp_probs()
    min_prob = np.min(probs)

    # Remove genes that have a low component frequency
    if min_prob < min_prob_filter:
        analyzed[name] = (name, False)
        return name, False

    elif len(probs) > 1:
        result = True

        # Make sure gene is multimodal with respect to covariate
        if covariate is not None:
            assert X.dim == 2, 'Covariate data is not two dimensional'
            result = distinct_covariates(X, model)

        # Save genes for future analysis
        if result is True and output_dir:
            _dir = os.path.join(output_dir, 'MultiModalGenes', name)
            mkdir_p(_dir)
            bnpy.ioutil.ModelWriter.save_model(model,
                                               _dir,
                                               prefix=name)
        analyzed[name] = (name, result)
        return name, result

    else:
        analyzed[name] = (name, False)
        return name, False


def filter_geneset(lst, matrx, covariate=None, CPU=1, gene_mean_filter=None, min_prob_filter=None, output_dir=False):
    """

    :param lst: List of genes
    :param matrx: Expression dataframe
    :return:
    """
    pool = multiprocessing.Pool(processes=CPU)

    results = []
    for gene in lst:
        data = matrx.loc[gene, :].values
        mean = np.mean(data)

        # Skip Genes that have a mean below the mean filter
        if mean < gene_mean_filter:
            continue

        data = data - mean
        data = data.reshape(len(data), 1)

        if covariate is not None:
            data = np.hstack((data, covariate))
            data = data.reshape(len(data), 2)

        #res = is_multimodal(gene, data, covariate, min_prob_filter, output_dir)

        #print res

        res = pool.apply_async(is_multimodal, args=(gene, data, covariate, min_prob_filter, output_dir,))
        results.append(res)

    output = [x.get() for x in results]

    print output

    # Remove duplicate genes
    return list(set([x[0] for x in output if x[1] is True]))


def find_aliases(gss, mapper, index):
    """
    Removes genes that do not overlap with matrix. Will
    try to rescue genes by mapping them to an alias that
    overlaps.

    :param gss:
    :param mapper:
    :param index:
    :return:
    """
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
                                logging.info("Replacing %s in %s with %s" % (gene,
                                                                             gs,
                                                                             match))
                                filtered.append(match)
        # Remove duplicates
        gss[gs] = list(set(filtered))
    return gss


def main():
    """
    Fits one, two, and three component mixture models and returns fit information in output dir.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)

    parser.add_argument('-e', '--expression',
                        help='Gene symbol by sample matrix.\nDo not center the data.',
                        required=True)

    parser.add_argument('-c', '--covariate',
                        help='sample X covariate matrix.',
                        required=False)

    parser.add_argument('--CPU',
                        dest='CPU',
                        type=int,
                        default=1)

    parser.add_argument('--hallmark',
                        help='Uses hallmark gene sets',
                        dest='hallmark',
                        default=False,
                        action='store_true')

    parser.add_argument('--msigdb',
                        help='Uses MSigDB gene sets',
                        dest='msigdb',
                        default=False,
                        action='store_true')

    parser.add_argument('--reactome',
                        help='Uses Reactome gene sets',
                        dest='reactome',
                        default=False,
                        action='store_true')

    parser.add_argument('--treehouse',
                        help='Uses Treehouse gene sets',
                        dest='treehouse',
                        default=False,
                        action='store_true')

    parser.add_argument('--immune',
                        help='Uses curated immune gene sets',
                        dest='immune',
                        default=False,
                        action='store_true')

    parser.add_argument('--all-genes',
                        help='Uses all genes in expression matrix',
                        dest='all_genes',
                        default=False,
                        action='store_true')

    parser.add_argument('--save-genes',
                        help='Saves multimodal gene fits',
                        dest='save_genes',
                        default=False,
                        action='store_true')

    parser.add_argument('--min-mean-filter',
                        dest='min_mean_filter',
                        help='Removes genes with an average expression below value.',
                        type=float,
                        default=None)

    parser.add_argument('--min-prob-filter',
                        dest='min_prob_filter',
                        help='Removes genes with a minor component less than value.',
                        type=float,
                        default=None)

    parser.add_argument('--min-num-filter',
                        dest='min_num_filter',
                        help='Skips gene sets with fewer than X multimodal genes.',
                        type=int,
                        default=5)

    parser.add_argument('--debug',
                        action='store_true')

    parser.add_argument('--output-dir',
                        dest='output_dir',
                        default='hydra-out')

    args = parser.parse_args()

    mkdir_p(args.output_dir)

    logging.basicConfig(filename=os.path.join(args.output_dir, 'hydra.log'),
                        level=logging.INFO)

    logging.getLogger().addHandler(logging.StreamHandler())

    logging.info("Started Hydra...")

    logging.info("Parameters:")
    for key, value in vars(args).items():
        logging.info('\t%s: %s' % (key, value))

    matrx = pd.read_csv(args.expression,
                        sep='\t',
                        index_col=0)

    # Remove duplicates in index
    matrx = matrx[~matrx.index.duplicated(keep='first')]

    # Determine which gene sets are included.
    if args.debug:
        sets, gs_map = get_test_genesets()
    	genesets = read_genesets(sets)

    elif args.all_genes:
        gs_map = {'ALL_GENES': 'ALL_GENES' }
        genesets = {'ALL_GENES': set(matrx.index.values)}

    else:
        dirs = ['misc']

        if args.msigdb:
            dirs = ['msigdb'] + dirs

        if args.reactome:
            dirs = ['reactome'] + dirs

        if args.hallmark:
            dirs = ['hallmark'] + dirs

        if args.treehouse:
            dirs = ['treehouse'] + dirs

        if args.immune:
            dirs = ['immune'] + dirs

        sets, gs_map = get_genesets(dirs)

        if len(sets) == 0:
            raise ValueError("Need to specify gene sets for analysis.")

    	genesets = read_genesets(sets)


    # Find overlap in alias space
    pth = os.path.join(src, 'data/alias-mapper.gz')
    alias_mapper = pd.read_csv(pth, sep='\t')
    genesets = find_aliases(genesets, alias_mapper, matrx.index)

    if args.min_mean_filter is not None:
        logging.info("Minimum gene expression mean: %0.2f" % args.min_mean_filter)

    if args.min_prob_filter is not None:
        logging.info("Minimum component probability: %0.2f" % args.min_prob_filter)

    if args.covariate is not None:
        logging.info("Reading in covariate:\n%s" % args.covariate)
        covariate = pd.read_csv(args.covariate,
                                sep='\t',
                                header=None,
                                index_col=0)

        covariate = covariate.reindex(matrx.columns).dropna()
        matrx = matrx.reindex(covariate.index, axis='columns').dropna()

        # Center the covariate data
        covariate = covariate.values
        covariate = covariate - np.mean(covariate)
        covariate.shape = (len(covariate), 1)

    filtered_genesets = defaultdict(set)
    for gs, genes in genesets.items():

        start = len(genes)
        filtered_genesets[gs] = filter_geneset(list(genes),
                                               matrx,
                                               covariate=covariate,
                                               CPU=args.CPU,
                                               gene_mean_filter=args.min_mean_filter,
                                               min_prob_filter=args.min_prob_filter,
                                               output_dir=args.output_dir)
        end = len(filtered_genesets[gs])

        logging.info("Filtering: {gs} went from {x} to {y} genes".format(gs=gs,
                                                                         x=start,
                                                                         y=end))
        if end < args.min_num_filter:
            logging.info("Skipping {gs} because there are not enough genes".format(gs=gs))
            continue

        gs_root = gs_map[gs]

        # Make directory for output
        gsdir = os.path.join(args.output_dir, gs_root, gs)
        mkdir_p(gsdir)

        # Center data to make inference easier
        center = matrx.loc[filtered_genesets[gs], :].apply(lambda x: x - x.mean(), axis=1)

        # Trying to figure out where the duplicate index is coming from
        center = center[~center.index.duplicated(keep='first')]

        # Save the expression data for future analysis
        pth = os.path.join(gsdir, 'expression.tsv')
        matrx.loc[filtered_genesets[gs]].to_csv(pth, sep='\t')

        # Need to take the transpose
        # Samples x Genes
        data = center.T.values

        # Create dataset object for inference
        dataset = bnpy.data.XData(data)

        gamma = 1.0
        sF = 0.5
        K = 5

        logging.info("Multivariate Model Params:\ngamma: %.2f\nsF: %.2f\nK: %d" % (gamma, sF, K))

        # Fit multivariate model
        hmodel = parallel_fit(gs, dataset, gamma, sF, K, save_output=args.debug)

        # Get the sample assignments
        LP = hmodel.calc_local_params(dataset)
        asnmts = LP['resp'].argmax(axis=1)

        pth = os.path.join(gsdir, 'assignments.tsv')
        with open(pth, 'w') as f:
            for sample, asnmt in zip(center.columns, asnmts):
                f.write('%s\t%s\n' % (sample, asnmt))

        # Save model
        bnpy.ioutil.ModelWriter.save_model(hmodel,
                                           gsdir,
                                           prefix=gs)

        create_notebook(gs, gsdir)

        if args.debug:
            break

if __name__ == '__main__':
    main()
