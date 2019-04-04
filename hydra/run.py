#!/usr/bin/env python2.7
import argparse
import bnpy
import datetime
import logging
import multiprocessing
import numpy as np
import os
import pandas as pd

from collections import defaultdict
from scipy.stats import kruskal

from library.analysis import EnrichmentAnalysis, MultivariateMixtureModel
from library.fit import subprocess_fit
from library.utils import mkdir_p, get_genesets, get_test_genesets, read_genesets
from library.notebook import create_notebook

date = str(datetime.date.today())

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
            # Currently I'm keeping genes that fall into this category
            # but maybe I should consider skipping them...
            if len(cov_groups) == 1:
                logger.info("%s %s covariate data clustered into one group." % (gene, cov))

        if k is True:
            keeps.append(found_cov)

        elif k is False:
            keeps.append(not found_cov)

    return keeps


# Global variable for keeping track of multimodally expressed genes
manager = multiprocessing.Manager()
analyzed = manager.dict()


def is_multimodal(gene,
                  matrx,
                  covariate=None,
                  gene_mean_filter=None,
                  min_prob_filter=None,
                  output_dir=None,
                  save_genes=False,
                  alpha=0.01,
                  sensitive=False):
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
    :param sensitive
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

    bstart = 10
    mstart = 10
    dstart = 20
    K = 1
    sF = 2.0

    # Run with sensitive parameterization
    if sensitive is True:
        bstart = 0
        mstart = 2
        dstart = 2
        K = 5
        sF = 1.0

    # Run the parallel fit model
    # This is parameterized to be sensitive about
    # identifying multimodally expressed distributions
    model, converged, params, stdout = subprocess_fit(gene,
                                                      X,
                                                      gamma=5.0,
                                                      K=K,
                                                      sF=sF,
                                                      bstart=bstart,
                                                      mstart=mstart,
                                                      dstart=dstart)

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

            pth = os.path.join(_dir, "PARAMS")
            with open(pth, "w") as f:
                f.write(params)

            pth = os.path.join(_dir, "STDOUT")
            with open(pth, "w") as f:
                f.write(stdout)

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
                   save_genes=False,
                   sensitive=False):
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


    logger = logging.getLogger('root')

    if sensitive:
        logger.debug('SENSITIVE MODE...')

    results = []
    for gene in lst:

        # Determine if gene and covariate is multimodal
        res = pool.apply_async(is_multimodal, args=(gene,
                                                    matrx,
                                                    covariate,
                                                    gene_mean_filter,
                                                    min_prob_filter,
                                                    output_dir,
                                                    save_genes,
                                                    0.01,
                                                    sensitive))
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


def get_assignments(model, data):
    """
    Takes model and data and classifies samples

    Will label samples with -1 cluster if they do not
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

def run_notebook():
    import subprocess
    p = subprocess.check_call(['jupyter',
                               'notebook',
                               '--allow-root',
                               '--ip', '0.0.0.0',
                               '--no-browser'])

def apply_multivariate_model(input, args, output, name='MultivariateModel'):
    # Center data to make inference easier

    logging.info("Centering input to multivariate clustering.")
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

    logging.info("Multivariate Model Params:\n"
                 "gamma: %.2f\n"
                 "sF: %.2f\n"
                 "K: %d\n"
                 "nLaps: %d" % (gamma, sF, K, nLap))

    # Fit multivariate model
    hmodel, converged, params, stdout = subprocess_fit(name,
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

    create_notebook(src, name, output)

    with open(os.path.join(output, 'PARAMS'), 'w') as f:
        f.write(params)

    with open(os.path.join(output, 'STDOUT'), 'w') as f:
        f.write(stdout)


def run_filter_gene_set(genesets, matrx, args, covariate=None):
    filtered_genesets = defaultdict(set)
    for gs, genes in genesets.items():

        start = len(genes)
        res = filter_geneset(list(genes),
                             matrx,
                             covariate=covariate,
                             CPU=args.CPU,
                             gene_mean_filter=args.min_mean_filter,
                             min_prob_filter=args.min_prob_filter,
                             output_dir=args.output_dir,
                             save_genes=args.save_genes,
                             sensitive=args.sensitive)
        end = len(res)

        logging.info("Filtering: {gs} went from {x} to {y} genes".format(gs=gs,
                                                                         x=start,
                                                                         y=end))
        if end < args.min_num_filter:
            logging.info("Skipping {gs} because there "
                         "are not enough genes to cluster".format(gs=gs))
            continue

        filtered_genesets[gs] = res

    return filtered_genesets


def gene_set_clustering(genesets, matrx, args, covariate=None):
    """
    Iterate over the gene sets and select for multimodally expressed genes

    :param genesets:
    :param matrx:
    :param args:
    :param covariate:
    :return:
    """

    filtered_genesets = run_filter_gene_set(genesets, matrx, args, covariate)
    for gs, genes in filtered_genesets.items():
        # Make directory for output
        gsdir = os.path.join(args.output_dir, 'OUTPUT', gs)
        mkdir_p(gsdir)

        # Create training data for the model fit
        training = matrx.reindex(genes)

        # Save the expression data for future analysis
        pth = os.path.join(gsdir, 'training-data.tsv')
        training.to_csv(pth, sep='\t')

        apply_multivariate_model(training, args, gsdir, gs)


def enrichment_analysis(matrx, args, covariate=None, reanalysis=False):

    assert args.save_genes is True, 'Need to save genes for enrichment analysis'

    genes = matrx.index.values
    genesets = {'AllGenes': genes}

    if not reanalysis:
        _ = run_filter_gene_set(genesets, matrx, args, covariate)
        mm_genes = os.path.join(args.output_dir, 'MultiModalGenes')

    else:
        mm_genes = reanalysis

    post = EnrichmentAnalysis(mm_genes,
                              args.expression,
                              args.min_prob_filter,
                              gmt=args.gmt)

    output = os.path.join(args.output_dir, 'EnrichmentAnalysis', date)
    mkdir_p(output)

    terms = post.get_enriched_terms()

    if terms.shape[0] == 0:
        logging.info("No enriched terms were found!\n"
                     "Try decreasing the min-prob-filter parameter.\n"
                     "Repeat enrichment analysis using jupyter notebook")
        return

    else:
        pth = os.path.join(output, 'EnrichedTerms')
        terms.to_csv(pth, sep='\t')

    training = matrx.reindex(post.get_enriched_term_genes())

    pth = os.path.join(output, 'training-data.tsv')
    training.to_csv(pth, sep='\t')

    apply_multivariate_model(training, args, output)


def main():
    """
    Fits non-parametric mixture models to expression data. Can add a covariate to filter
    genes that are differentially expressed with respect to another continuous variable.
    Finally, there is support for variants, but the variants must be in a format that
    is similar to expression data. One can sample from a normal distribution to introduce
    noise to the variant status.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)

    parser.add_argument('mode',
                        help='run - starts clustering pipeline\n'
                             'notebook - spins up jupyter notebook\n'
                             'reanalyze - repeats multivariate clustering')

    parser.add_argument('-e', '--expression',
                        help='Gene symbol by sample matrix.\nDo not center the data.',
                        required=True)

    parser.add_argument('--go-enrichment',
                        help='Performs GO enrichment analysis using gene set database',
                        dest='go_enrichment',
                        action='store_true',
                        required=False)

    parser.add_argument('--enrichment',
                        help='Performs enrichment analysis using gene set database',
                        action='store_true',
                        required=False)

    parser.add_argument('--gmt',
                        help='Gene set database in GMT format.',
                        required=False)

    parser.add_argument('--all-genes',
                        help='Uses all multimodal genes in expression matrix',
                        dest='all_genes',
                        default=False,
                        action='store_true')

    parser.add_argument('--save-genes',
                        help='Saves multimodal gene fits',
                        dest='save_genes',
                        default=True,
                        action='store_true')

    parser.add_argument('--min-mean-filter',
                        dest='min_mean_filter',
                        help='Removes genes with an average expression below value.'
                             'Used in gene-wise clustering step.',
                        type=float,
                        default=None)

    parser.add_argument('--min-prob-filter',
                        dest='min_prob_filter',
                        help='Removes genes with a minor component less than value. '
                             'Used in gene-wise clustering step.',
                        type=float,
                        default=None)

    parser.add_argument('--min-num-filter',
                        dest='min_num_filter',
                        help='Skips gene sets with fewer than X multimodal genes.',
                        type=int,
                        default=5)

    parser.add_argument('-K',
                        help='Number of clusters to start with for multivariate clustering',
                        type=int,
                        default=5)

    parser.add_argument('--gamma',
                        help='Prior for dirichlet dispersion parameter gamma.',
                        type=float,
                        default=5.0)

    parser.add_argument('--sF',
                        help='Prior for diagonal of covariance matrix.',
                        type=float,
                        default=2.0)

    parser.add_argument('--num-laps',
                        dest='num_laps',
                        help='Number of laps for VB algorithm.',
                        type=int,
                        default=1000)

    parser.add_argument('--max-fit-time',
                        dest='max_fit_time',
                        help='Maximum number of seconds for multivariate fit before timeout',
                        type=int,
                        default=900)

    parser.add_argument('--sensitive',
                        help='Runs univariate filter in sensitive mode.',
                        action='store_true')

    parser.add_argument('--output-dir',
                        dest='output_dir',
                        default='hydra-out')

    parser.add_argument('--CPU',
                        dest='CPU',
                        type=int,
                        default=1)

    parser.add_argument('--test',
                        action='store_true')

    parser.add_argument('--debug',
                        action='store_true')

    parser.add_argument('-v', '--variants',
                        help='Variant identifier by sample matrix.\nName must be distinct from expression.',
                        required=False)

    parser.add_argument('-c', '--covariate',
                        help='sample X covariate matrix. There must be a row "keep" that specifies whether to keep '
                             'or remove genes that correlate with variable using a boolean. 1: keep. 0: remove.',
                        required=False)

    parser.add_argument('--gmt-regex',
                        dest='gmt_regex',
                        help='Regex for subsetting gmt file',
                        required=False)

    parser.add_argument('--mm-path',
                        dest='mm_path',
                        help='Path to MultiModal directory for reanalysis',
                        required=False)

    args = parser.parse_args()

    if args.mode == 'notebook':
        run_notebook()

    elif args.mode == 'reanalyze':
        raise NotImplementedError


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
    logging.info("Number of genes: %d" % matrx.shape[0])
    matrx = matrx[~matrx.index.duplicated(keep='first')]
    logging.info("Number of genes after removing duplicates: %d" % matrx.shape[0])

    drop_genes = set()
    for gene in matrx.index:
        if '/' in gene:
            logging.info("Gene names cannot contain forward slashes!")
            drop_genes.add(gene)

        if "'" in gene or "\"" in gene:
            logging.info("Gene names cannot contain quotation marks!")
            drop_genes.add(gene)

    matrx = matrx.drop(list(drop_genes), axis=0)
    logging.info("Number of genes after "
                 "removing misformatted genes: %d" % matrx.shape[0])

    # Determine which gene sets are included.
    if args.test:
        logging.info("Loading debug gene sets...")
        sets, _ = get_test_genesets(src)
        genesets = read_genesets(sets)
        args.gmt = 'TEST'

    elif args.gmt:
        genesets = get_genesets(args.gmt, args.gmt_regex)

        if len(genesets) == 0:
            raise ValueError("Need to specify gene sets for analysis.")

    else:
        logging.warn("No gene set database given.")
        genesets = {'ALL_GENES': set(matrx.index.values)}

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

    if args.enrichment:
        assert args.gmt is not None, "Need GMT file for enrichment analysis"
        logging.info("Starting unsupervised enrichment analysis:\n%s" % args.gmt)
        enrichment_analysis(matrx, args, covariate=covariate, reanalysis=args.mm_path)

    elif args.go_enrichment:
        logging.info("Starting unsupervised enrichment analysis using clusterProfiler GO analysis")
        args.gmt = "GO"
        enrichment_analysis(matrx, args, covariate=covariate, reanalysis=args.mm_path)

    else:
        logging.info("Starting gene set clustering analysis:\n%s" % args.gmt)
        gene_set_clustering(genesets, matrx, args, covariate)



if __name__ == '__main__':
    main()
