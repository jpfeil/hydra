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

from library.fit import subprocess_fit
from library.utils import mkdir_p, get_genesets, get_test_genesets, read_genesets
from library.notebook import create_notebook


src = os.path.dirname(os.path.abspath(__file__))


def distinct_covariates(gene, data, model, covariate, keep=None, alpha=0.01, debug=False):
    """
    Determine if clustering separates data into statistically different components

    :param data: XData object
    :param model: bnpy model
    :param covariate: numpy array in the same order as the data object
    """

    if keep is None:
        raise ValueError("Need to specify whether or not to keep gene based on covariate!")

    logger = logging.getLogger('root')

    # Fit the data again to determine the assignments
    LP = model.calc_local_params(data)
    asnmts = LP['resp'].argmax(axis=1)

    keeps = []
    for k, cov in zip(keep, covariate.columns):

        # Create a separate group for each assignment
        cov_groups = [[] for _ in range(max(asnmts) + 1)]

        # Add data to the appropriate elements in the list
        for c, a in zip(covariate[cov], asnmts):
            if pd.isnull(c):
                continue

            cov_groups[a].append( float(c) )

        # Remove groups with a small sample size
        cov_groups = [x for x in cov_groups if len(x) > 5]

        # Determine if the covariate data is statistically
        # different by a non-parametric Kruskal-Wallis test
        found_cov = False
        try:
            _, cov_p = kruskal(*cov_groups, nan_policy='raise')
            logger.debug("Gene: %s, Covariate: %s, P-value: %.8f, Alpha: %.8f" % (gene,
                                                                                  cov,
                                                                                  cov_p,
                                                                                  alpha))

            # Kruskal-Wallis is an omnibus test which apparently controls for the
            # false positive rate:
            # https://stats.stackexchange.com/questions/133444/bonferroni-correction-on-multiple-kruskal-wallis-tests

            if cov_p < alpha:
                found_cov = True

            elif pd.isnull(cov_p):
                raise ValueError("Kruskal-Wallis p-value is NaN")

        except ValueError:
            if len(cov_groups) == 1:
                logger.info("%s covariate data clustered into one group." % cov)

        if k is True:
            keeps.append(found_cov)

        elif k is False:
            keeps.append(not found_cov)

    return keeps


# Global variable for keeping track of multimodally expressed genes
analyzed = {}


def is_multimodal(gene,
                  matrx,
                  covariate=None,
                  gene_mean_filter=None,
                  min_prob_filter=None,
                  output_dir=None,
                  save_genes=False,
                  alpha=0.01):
    """
    This function determines if there is a multimodal pattern in the data. Also has a bunch of other
    functions that should be factored out. For example, this function also skips genes that do not
    pass the minimum probability filter.

    :param gene: gene name
    :param data: expression data
    :param covariate: covariate data
    :param min_prob_filter: minimium probability filter
    :param output_dir: path to output directory
    :param save_genes: whether to save the gene fit
    :param alpha: significance threshold
    :return:
    """
    # If we have analyzed this sample before,
    # then just take that value and save some time
    if gene in analyzed:
        return analyzed[gene]

    logger = logging.getLogger('root')

    data = matrx.loc[gene, :].values
    mean = np.mean(data)

    # Skip Genes that have a mean below the mean filter
    if mean < gene_mean_filter:
        logger.debug("Gene %s was removed by mean filter." % gene)
        return gene, False

    # Center the expression data. This data should not be
    # centered before this step
    data = data - mean

    # Reshape the data so that it is compatible with the
    # mixture model
    data = data.reshape(len(data), 1)

    keeps = None
    if covariate is not None:
        keeps = [bool(x) for x in covariate.loc['keep', :].values]
        covariate = covariate.reindex(matrx.columns)

    # Convert data to a bnpy XData object
    X = bnpy.data.XData(data)

    # Run the parallel fit model
    model, converged = subprocess_fit(gene, X, gamma=5.0, K=3)
    probs = model.allocModel.get_active_comp_probs()
    min_prob = np.min(probs)

    # Do not consider genes where the model did not converge
    if converged is False:
        logger.info("Model did not converge for %s" % gene)
        return gene, False

    # Remove genes that have a low component frequency
    if min_prob < min_prob_filter:
        logger.info("Gene %s was removed by min prob filter." % gene)
        analyzed[gene] = (gene, False)
        return gene, False

    elif len(probs) > 1:
        result = True

        # Make sure gene is multimodal with respect to covariate
        if covariate is not None:
            results = distinct_covariates(gene,
                                          X,
                                          model,
                                          covariate,
                                          keep=keeps,
                                          alpha=alpha)

            for r, k, name in zip(results, keeps, covariate.columns):
                if r is False:
                    logger.debug("%s was removed by correlation filter %s" % (gene, name))
                    result = False

                elif r is True:
                    logger.debug("%s passed correlation filter %s" % (gene, name))

                else:
                    raise ValueError()

        # Save genes for future analysis
        if result is True and output_dir and save_genes:
            _dir = os.path.join(output_dir, 'MultiModalGenes', gene)
            mkdir_p(_dir)
            bnpy.ioutil.ModelWriter.save_model(model,
                                               _dir,
                                               prefix=gene)

        analyzed[gene] = (gene, result)
        return gene, result

    else:
        logger.debug("Gene %s was removed because it is unimodal" % gene)
        analyzed[gene] = (gene, False)
        return gene, False


def filter_geneset(lst,
                   matrx,
                   covariate=None,
                   CPU=1,
                   gene_mean_filter=None,
                   min_prob_filter=None,
                   output_dir=None,
                   save_genes=False):
    """
    Applies non-parametric mixture model to expression data. Can optionally add
    a covariate (e.g. survival, IC50) to select genes that vary with a variable of
    interest.

    Loops over a list of genes and selects rows in the expression matrix. Creates
    a bivariate data set for analyzing genes that covary with a variable of interest.

    :param lst (list): list of genes
    :param matrx (pd.DataFrame): expression dataframe
    :param covariate (np.array): covariate vector
    :param CPU (int): number of CPUs available
    :param gene_mean_filter (float): mean threshold for filtering genes
    :param min_prob_filter (float): min probability threshold for filtering genes
    :param output_dir (str): path to output directory for saving intermediate files
    :param save_genes (bool): flag for deciding whether or not to save the gene models
    :return:
    """
    pool = multiprocessing.Pool(processes=CPU)

    results = []
    for gene in lst:

        # Determine if gene and covariate is multimodal
        res = pool.apply_async(is_multimodal, args=(gene,
                                                    matrx,
                                                    covariate,
                                                    gene_mean_filter,
                                                    min_prob_filter,
                                                    output_dir,
                                                    save_genes))
        results.append(res)

    output = [x.get() for x in results]

    # Select multimodal genes
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
    Fits non-parametric mixture models to expression data. Can add a covariate to filter
    genes that are differentially expressed with respect to another continuous variable.
    Finally, there is support for variants, but the variants must be in a format that
    is similar to expression data. One can sample from a normal distribution to introduce
    noise to the variant status.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)

    parser.add_argument('-e', '--expression',
                        help='Gene symbol by sample matrix.\nDo not center the data.',
                        required=True)

    parser.add_argument('-v', '--variants',
                        help='Variant identifier by sample matrix.\nName must be distinct from expression.',
                        required=False)

    parser.add_argument('-c', '--covariate',
                        help='sample X covariate matrix. There must be a row "keep" that specifies whether to keep '
                             'or remove genes that correlate with variable using a boolean. 1: keep. 0: remove.',
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

    parser.add_argument('--num-laps',
                        dest='num_laps',
                        help='Number of laps for VB algorithm.',
                        type=int,
                        default=500)

    parser.add_argument('--test',
                        action='store_true')

    parser.add_argument('--debug',
                        action='store_true')

    parser.add_argument('--output-dir',
                        dest='output_dir',
                        default='hydra-out')

    args = parser.parse_args()

    # Make the output directory if it doesn't already exist
    mkdir_p(args.output_dir)

    # Start logging
    level = logging.INFO
    if args.debug:
        level = logging.DEBUG

    logging.basicConfig(filename=os.path.join(args.output_dir, 'hydra.log'),
                        level=level)

    logging.getLogger().addHandler(logging.StreamHandler())

    logging.info("Started Hydra...")

    logging.info("Parameters:")
    for key, value in vars(args).items():
        logging.info('\t%s: %s' % (key, value))

    # Read in expression data
    logging.info("Reading in expression data:\n%s" % args.expression)
    matrx = pd.read_csv(args.expression,
                        sep='\t',
                        index_col=0)

    # Remove duplicates in index
    logging.info("Removing duplicate genes:")
    matrx = matrx[~matrx.index.duplicated(keep='first')]

    for gene in matrx.index:
        if '/' in gene:
            raise ValueError("Gene names cannot contain forward slashes!")

    # Determine which gene sets are included.
    if args.test:
        logging.info("Loading debug gene sets...")
        sets, gs_map = get_test_genesets(src)
        genesets = read_genesets(sets)

    elif args.all_genes:
        logging.info("Creating ALL_GENES gene set...")
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

        sets, gs_map = get_genesets(dirs, src)

        if len(sets) == 0:
            raise ValueError("Need to specify gene sets for analysis.")

        genesets = read_genesets(sets)

    # Find overlap in alias space
    pth = os.path.join(src, 'data/alias-mapper.gz')
    alias_mapper = pd.read_csv(pth, sep='\t')
    logging.info("Looking for gene aliases...")
    genesets = find_aliases(genesets, alias_mapper, matrx.index)

    if args.min_mean_filter is not None:
        logging.info("Minimum gene expression mean: %0.2f" % args.min_mean_filter)

    if args.min_prob_filter is not None:
        logging.info("Minimum component probability: %0.2f" % args.min_prob_filter)

    covariate = None
    if args.covariate is not None:
        logging.info("Reading in covariate:\n%s" % args.covariate)
        covariate = pd.read_csv(args.covariate,
                                sep='\t',
                                index_col=0)

        if 'keep' not in covariate.index:
            raise ValueError("Need a keep row in the covariate matrix!")

        #
        covariate.loc['keep', :] = covariate.loc['keep', :].astype('bool')

    # Iterate over the gene sets and select for multimodally expressed genes
    filtered_genesets = defaultdict(set)
    for gs, genes in genesets.items():

        start = len(genes)
        filtered_genesets[gs] = filter_geneset(list(genes),
                                               matrx,
                                               covariate=covariate,
                                               CPU=args.CPU,
                                               gene_mean_filter=args.min_mean_filter,
                                               min_prob_filter=args.min_prob_filter,
                                               output_dir=args.output_dir,
                                               save_genes=args.save_genes)
        end = len(filtered_genesets[gs])

        logging.info("Filtering: {gs} went from {x} to {y} genes".format(gs=gs,
                                                                         x=start,
                                                                         y=end))
        if end < args.min_num_filter:
            logging.info("Skipping {gs} because there are not enough genes to cluster".format(gs=gs))
            continue

        gs_root = gs_map[gs]

        # Make directory for output
        gsdir = os.path.join(args.output_dir, gs_root, gs)
        mkdir_p(gsdir)

        # Create training data for the model fit
        training = matrx.loc[filtered_genesets[gs], :]

        # Save the expression data for future analysis
        pth = os.path.join(gsdir, 'training-data.tsv')
        training.to_csv(pth, sep='\t')

        # Center data to make inference easier
        center = training.apply(lambda x: x - x.mean(), axis=1)

        # Need to take the transpose
        # Samples x Genes
        data = center.T.values

        # Create dataset object for inference
        dataset = bnpy.data.XData(data)

        # Set the prior for creating a new cluster
        gamma = 5.0

        # Start with a standard identity matrix
        sF = 1.0

        # Starting with one cluster because most
        # distributions will likely have only one cluster.
        K = 3

        logging.info("Multivariate Model Params:\ngamma: %.2f\nsF: %.2f\nK: %d" % (gamma, sF, K))

        # Fit multivariate model
        hmodel, converged = subprocess_fit(gs,
                                           dataset,
                                           gamma,
                                           sF,
                                           K,
                                           save_output=args.save_genes)

        if converged is False:
            logging.info("WARNING: Multivariate model did not converge!")

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
