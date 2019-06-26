#!/usr/bin/env python2.7
import argparse
import bnpy
import datetime
import logging
import multiprocessing
import numpy as np
import os
import pandas as pd
import subprocess
import sys
import textwrap
import re

from collections import defaultdict

from library.analysis import EnrichmentAnalysis
from library.fit import subprocess_fit, apply_multivariate_model
from library.notebook import enrichment_notebook
from library.utils import mkdir_p, \
                          get_genesets, \
                          get_test_genesets, \
                          read_genesets, \
                          find_aliases

date = str(datetime.date.today())

src = os.path.dirname(os.path.abspath(__file__))

# Global variable for keeping track of multimodally expressed genes
manager = multiprocessing.Manager()
analyzed = manager.dict()


def filter_gene(gene,
                matrx,
                gene_mean_filter=None,
                min_prob_filter=None,
                output_dir=None,
                sensitive=False):
    """
    Filters genes using multiomodal and simple expression filters.

    :param gene: gene name
    :param data: expression data
    :param covariate: covariate data
    :param min_prob_filter: minimium probability filter
    :param output_dir: path to output directory
    :param sensitive
    :return:
    """

    logger = logging.getLogger('root')

    # If we have analyzed this sample before,
    # then just take that value and save some time
    if gene in analyzed:
        logger.debug("Gene %s was previously analyzed" % gene)
        return analyzed[gene]

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
    elif min_prob < min_prob_filter:
        logger.info("Gene %s was removed by min prob filter." % gene)
        analyzed[gene] = (gene, False)
        return gene, False

    elif len(probs) > 1:
        result = True

        if result is True and output_dir:
            _dir = os.path.join(output_dir, 'MultiModalGenes', gene)
            assert not os.path.exists(_dir), 'Tried to overwrite previous gene-level fit!'
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

    elif len(probs) == 1:
        logger.debug("Gene %s was removed because it is unimodal" % gene)
        analyzed[gene] = (gene, False)
        return gene, False

    else:
        raise ValueError()


def filter_gene_set(lst,
                    matrx,
                    CPU=1,
                    gene_mean_filter=None,
                    min_prob_filter=None,
                    output_dir=None,
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
    :param sensitive (bool): Boolean to determine if mixture model is run in sensitive mode
    :return:
    """
    logger = logging.getLogger('root')

    pool = multiprocessing.Pool(processes=CPU)

    if sensitive:
        logger.debug('SENSITIVE MODE...')

    results = []
    for gene in lst:

        # Determine if gene and covariate is multimodal
        res = pool.apply_async(filter_gene, args=(gene,
                                                  matrx,
                                                  gene_mean_filter,
                                                  min_prob_filter,
                                                  output_dir,
                                                  sensitive))
        results.append(res)

    output = [x.get() for x in results]

    # Select multimodal genes
    return list(set([x[0] for x in output if x[1] is True]))


def filter_gene_sets(genesets, matrx, args):
    """
    Takes gene sets and filters down to the multimodally expressed genes

    :param genesets (dict): Dictionary containing gene sets
    :param matrx (pandas.DataFrame): Expression dataframe
    :param args (argparse.Namespace): Input parameters
    :return: Filtered gene sets
    :rtype: collections.defaultdict
    """
    logger = logging.getLogger('root')
    filtered_genesets = defaultdict(set)
    for gs, genes in genesets.items():

        start = len(genes)
        res = filter_gene_set(list(genes),
                              matrx,
                              CPU=args.CPU,
                              gene_mean_filter=args.min_mean_filter,
                              min_prob_filter=args.min_prob_filter,
                              output_dir=args.output_dir,
                              sensitive=args.sensitive)
        end = len(res)

        logger.info("Filtering: {gs} went from {x} to {y} genes".format(gs=gs,
                                                                         x=start,
                                                                         y=end))
        if end < args.min_num_filter:
            logger.info("Skipping {gs} because there "
                         "are not enough genes to cluster".format(gs=gs))
            continue

        filtered_genesets[gs] = res

    return filtered_genesets


def gene_set_clustering(genesets, matrx, args):
    """
    Iterate over the gene sets and select for multimodally expressed genes

    :param genesets:
    :param matrx:
    :param args:
    :param covariate:
    :return:
    """

    filtered_genesets = filter_gene_sets(genesets, matrx, args)
    for gs, genes in filtered_genesets.items():
        # Make directory for output
        gsdir = os.path.join(args.output_dir, 'MultivariateAnalysis', gs)
        mkdir_p(gsdir)

        # Create training data for the model fit
        training = matrx.reindex(genes)

        # Save the expression data for future analysis
        pth = os.path.join(gsdir, 'training-data.tsv')
        training.to_csv(pth, sep='\t')

        apply_multivariate_model(training, args, gsdir, src, name=gs)


def enrichment_analysis(matrx, args):
    """
    Runs unsuperivsed enrichment analysis. Will run the multimodal filter
    unless a path to the multimodal output directory is provided.

    :param matrx:
    :param args:
    :return:
    """
    logger = logging.getLogger('root')
    genes = matrx.index.values
    if args.mm_path is None:
        start = len(genes)
        res = filter_gene_set(list(genes),
                              matrx,
                              CPU=args.CPU,
                              gene_mean_filter=args.min_mean_filter,
                              min_prob_filter=args.min_prob_filter,
                              output_dir=args.output_dir,
                              sensitive=args.sensitive)
        end = len(res)

        logger.info("Filtering: {gs} went from {x} to {y} genes".format(gs='Enrichment Analysis',
                                                                        x=start,
                                                                        y=end))
        _ = filter_gene_set(list(genes),
                             matrx,
                             CPU=args.CPU,
                             gene_mean_filter=args.min_mean_filter,
                             min_prob_filter=args.min_prob_filter,
                             output_dir=args.output_dir,
                             sensitive=args.sensitive)

        mm_genes = os.path.join(args.output_dir, 'MultiModalGenes')

    elif os.path.exists(args.mm_path):
        mm_genes = args.mm_path

    else:
        raise ValueError('Could not identify multimodally expressed genes for enrichment analysis')

    post = EnrichmentAnalysis(mm_genes,
                              args.expression,
                              min_comp_filter=args.min_prob_filter,
                              min_effect_filter=args.min_effect_filter,
                              gmt_path=args.gmt)

    output = os.path.join(args.output_dir, 'EnrichmentAnalysis', date)
    mkdir_p(output)

    terms = post.get_enriched_terms()

    if terms.shape[0] == 0:
        logger.info("No enriched terms were found!\n"
                    "Try decreasing the min-prob-filter parameter.\n"
                    "Repeat enrichment analysis using jupyter notebook")
        return

    else:
        pth = os.path.join(output, 'EnrichedTerms')
        terms.to_csv(pth, sep='\t')

    training = matrx.reindex(post.get_enriched_term_genes())

    pth = os.path.join(output, 'training-data.tsv')
    training.to_csv(pth, sep='\t')

    apply_multivariate_model(training, args, output, src)




def filter(matrx, args):
    """
    Runs the multimodal filter and stops

    :param matrx pandas.DataFrame: Expression matrix
    :param args argparse.Namespace: Pipeline configuration namespace
    :return: None
    """
    logger = logging.getLogger('root')
    logger.info("Running multimodal filter")
    genes = matrx.index.values
    _ = filter_gene_set(genes,
                        matrx,
                        CPU=args.CPU,
                        gene_mean_filter=args.min_mean_filter,
                        min_prob_filter=args.min_prob_filter,
                        output_dir=args.output_dir,
                        sensitive=args.sensitive)
    mm_genes = os.path.join(args.output_dir, 'MultiModalGenes')
    logger.info('Multimodal genes are located here:\n%s' % mm_genes)


def sweep(matrx, args):
    """
    Applies multivariate clustering of multimodally expressed genes across
    input gene set database. Use regex to select specific gene sets.

    :param matrx (pandas.DataFrame): Expression matrix
    :param args (argparse.Namespace):
    :return:
    """
    logger = logging.getLogger('root')
    # Determine which gene sets are included.
    if args.test:
        logger.info("Loading debug gene sets...")
        sets, _ = get_test_genesets(src)
        genesets = read_genesets(sets)
        args.gmt = 'TEST'

    elif args.gmt:
        genesets = get_genesets(args.gmt, args.gmt_regex)

        if len(genesets) == 0:
            raise ValueError("Need to specify gene sets for analysis.")

    else:
        logger.warn("No gene set database given.")
        genesets = {'ALL_GENES': set(matrx.index.values)}

    # Find overlap in alias space
    pth = os.path.join(src, 'data/alias-mapper.gz')
    alias_mapper = pd.read_csv(pth, sep='\t')
    logger.info("Looking for gene aliases...")
    genesets = find_aliases(genesets, alias_mapper, matrx.index)

    logger.info("Starting gene set clustering analysis:\n%s" % args.gmt)
    gene_set_clustering(genesets, matrx, args)


def enrich(matrx, args):
    """

    :param matrx:
    :param args:
    :return:
    """
    logger = logging.getLogger('root')
    if args.gmt is not None:
        logger.info("Starting unsupervised enrichment analysis:\n%s" % args.gmt)
        enrichment_analysis(matrx, args)

    elif args.gmt is None and args.go_enrichment:
        logger.info("Starting unsupervised enrichment analysis using clusterProfiler GO analysis")
        args.gmt = "GO"
        enrichment_analysis(matrx, args)

    else:
        raise ValueError("Either specify a GMT file or use the --go-enrichment flag!")


def run_notebook():

    sys.path.append(src)

    logger = logging.getLogger('root')

    regex = re.compile('token=(?P<token>\w*)')

    cmd = ['jupyter',
           'notebook',
           '--allow-root',
           '--ip', '0.0.0.0',
           '--no-browser']

    found = False

    try:
        logger.info("Starting Jupyter Notebook")
        logger.info('Ctrl-C to stop session...')
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)

        for stderr_line in iter(p.stderr.readline, ""):
            m = regex.search(stderr_line)

            if m and not found:
                logger.info("TOKEN: %s" % m.group('token'))
                found = True

        p.stderr.close()
        return_code = p.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    except KeyboardInterrupt:
        logger.info('\nClosing Jupyter Notebook')


def main():
    """
    Hydra pipeline for identifying multimodally expressed genes. Takes an expression
    matrix and outputs differentially expressed gene sets or identifies enrichment
    of gene sets for further analysis.
    """
    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('mode',
                        help=textwrap.dedent("""
                        filter - applied multimodal filter
                        enrich - performs unsupervised enrichment analysis
                        sweep - performs within gene set clustering
                        notebook - spins up jupyter notebook"""))

    parser.add_argument('-e', '--expression',
                        help='Gene symbol by sample matrix.\nDo not center the data.',
                        required=False)

    parser.add_argument('--gmt',
                        help='Gene set database in GMT format.',
                        required=False)

    parser.add_argument('-m', '--multimodal',
                        dest='mm_path',
                        help='Path to MultiModalGenes directory',
                        required=False)

    parser.add_argument('--go-enrichment',
                        help='Performs GO enrichment analysis using gene set database',
                        dest='go_enrichment',
                        action='store_true',
                        required=False)

    parser.add_argument('--all-genes',
                        help='Uses all multimodal genes in expression matrix',
                        dest='all_genes',
                        default=False,
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

    parser.add_argument('--min-gene-filter',
                        dest='min_num_filter',
                        help='Skips gene sets with fewer than X multimodal genes.',
                        type=int,
                        default=5)

    parser.add_argument('--min-effect-filter',
                        dest='min_effect_filter',
                        help='Removes genes with a minimum effect less than value.',
                        type=float,
                        default=None)

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

    parser.add_argument('-o', '--output-dir',
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

    parser.add_argument('--gmt-regex',
                        dest='gmt_regex',
                        help='Regex for subsetting gmt file',
                        required=False)

    args = parser.parse_args()

    # Set up logger
    level = logging.INFO
    if args.debug:
        level = logging.DEBUG

    # Make the output directory if it doesn't already exist
    mkdir_p(args.output_dir)

    logging.basicConfig(filename=os.path.join(args.output_dir, 'hydra.log'),
                        level=level)
    logging.getLogger().addHandler(logging.StreamHandler())

    logger = logging.getLogger('root')

    # Spin up notebook
    if args.mode == 'notebook':
        run_notebook()
        return

    logger.info("Started Hydra...")
    logger.info("Parameters:")
    for key, value in vars(args).items():
        logger.info('\t%s: %s' % (key, value))

    if args.min_mean_filter is not None:
        logger.info("Minimum gene expression mean: %0.2f" % args.min_mean_filter)

    if args.min_prob_filter is not None:
        assert args.min_prob_filter < 1.0, 'Probability filter must be expressed as a decimal!'
        logger.info("Minimum component probability: %0.2f" % args.min_prob_filter)

    if args.expression is None:
        raise ValueError('Path to expression file is required! (-e)')

    # Read in expression data
    logger.info("Reading in expression data:\n%s" % args.expression)
    matrx = pd.read_csv(args.expression,
                        sep='\t',
                        index_col=0)

    # Remove duplicates in index
    logger.info("Removing duplicate genes:")
    logger.info("Number of genes: %d" % matrx.shape[0])
    matrx = matrx[~matrx.index.duplicated(keep='first')]
    logger.info("Number of genes after removing duplicates: %d" % matrx.shape[0])

    drop_genes = set()
    for gene in matrx.index:
        if '/' in gene:
            logger.info("Gene names cannot contain forward slashes! %s" % gene)
            drop_genes.add(gene)

        if "'" in gene or "\"" in gene:
            logger.info("Gene names cannot contain quotation marks! %s" % gene)
            drop_genes.add(gene)

    matrx = matrx.drop(list(drop_genes), axis=0)
    logger.info("Number of genes after "
                "removing misformatted genes: %d" % matrx.shape[0])

    if args.mode == 'filter':
        filter(matrx, args)

    elif args.mode == 'sweep':
        sweep(matrx, args)

    elif args.mode == 'enrich':
        enrich(matrx, args)
        enrichment_notebook(args.expression,
                            args.output_dir)

if __name__ == '__main__':
    main()
