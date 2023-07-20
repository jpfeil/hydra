import bnpy
import logging
import os
import multiprocessing
import shutil
import numpy as np

from library.fit import subprocess_fit
from library.utils import mkdir_p

# Global variable for keeping track of multimodally expressed genes
manager = multiprocessing.Manager()
analyzed = manager.dict()

def fit_gene(gene, X, sensitive=False):
    # Default parameters
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
    return subprocess_fit(gene,
                          X,
                          gamma=5.0,
                          sF=sF,
                          K=K,
                          bstart=bstart,
                          mstart=mstart,
                          dstart=dstart)


def get_gene_fit(gene,
                 data,
                 min_prob_filter=None,
                 output_dir=None,
                 sensitive=False,
                 mm_path=None,
                 overwrite=False):
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
    logger.debug("\nFiltering %s" % gene)

    save = True

    converged = False
    params = None
    stdout = None

    if mm_path:
        model_pth = os.path.join(mm_path, gene)
        try:
            logger.info("Loading fit:\n%s" % model_pth)
            model = bnpy.ioutil.ModelReader.load_model_at_prefix(model_pth,
                                                                 prefix=gene)
            converged = True       # We don't save models that don't converge
            save = False

        except (IOError, ValueError):
            logger.info("Couldn't find model:\n%s" % model_pth)
            res = fit_gene(gene,
                                                                data,
                                                                sensitive=sensitive)
            if res is None:
                return gene, False 
            
            model, converged, params, stdout, stderr = res

    else:
        logger.debug("Fitting gene %s" % gene)
        res = fit_gene(gene,
                                                            data,
                                                            sensitive=sensitive)
        if res is None:
            return gene, False
        
        model, converged, params, stdout, stderr = res

    probs = model.allocModel.get_active_comp_probs()
    min_prob = np.min(probs)

    # Do not consider genes where the model did not converge
    if not converged:
        logger.info("Model did not converge for %s" % gene)
        return gene, False

    # Remove genes that have a low component frequency
    elif min_prob is None or min_prob < min_prob_filter:
        logger.info("Gene %s was removed by min prob filter." % gene)
        analyzed[gene] = (gene, False)
        return gene, False

    elif len(probs) > 1:
        if save and output_dir:
            _dir = os.path.join(output_dir, 'MultiModalGenes', gene)
            if os.path.exists(_dir) and overwrite:
                logger.info("Erasing previous fit %s" % gene)
                shutil.rmtree(_dir)

            # This line is causing problems
            # assert not os.path.exists(_dir), 'Tried to overwrite previous gene-level fit!'
            mkdir_p(_dir)
            bnpy.ioutil.ModelWriter.save_model(model,
                                               _dir,
                                               prefix=gene)

            pth = os.path.join(_dir, "PARAMS")
            with open(pth, "w") as f:
                f.write(str(params))

            pth = os.path.join(_dir, "STDOUT")
            with open(pth, "w") as f:
                f.write(str(stdout))

        logger.debug("Gene %s is multimodal" % gene)
        analyzed[gene] = (gene, True)
        return gene, True

    elif len(probs) == 1:
        logger.debug("Gene %s was removed because it is unimodal" % gene)
        analyzed[gene] = (gene, False)
        return gene, False

    else:
        raise ValueError("Couldn't place gene!")


def filter_gene_set(lst,
                    expression,
                    CPU=1,
                    gene_mean_filter=None,
                    min_prob_filter=None,
                    output_dir=None,
                    sensitive=False,
                    mm_path=None,
                    overwrite=False):
    """
    Applies non-parametric mixture model to expression data. Can optionally add
    a covariate (e.g. survival, IC50) to select genes that vary with a variable of
    interest.

    Loops over a list of genes and selects rows in the expression matrix. Creates
    a bivariate data set for analyzing genes that covary with a variable of interest.

    :param lst (list): list of genes
    :param covariate (np.array): covariate vector
    :param CPU (int): number of CPUs available
    :param gene_mean_filter (float): mean threshold for filtering genes
    :param min_prob_filter (float): min probability threshold for filtering genes
    :param output_dir (str): path to output directory for saving intermediate files
    :param sensitive (bool): Boolean to determine if mixture model is run in sensitive mode
    :return:
    """
    logger = logging.getLogger('root')
    logger.debug("Filtering genes...")

    if sensitive:
        logger.debug('SENSITIVE MODE...')

    pool = multiprocessing.Pool(CPU)
    results = []
    for gene in lst:
        # If we have analyzed this sample before,
        # then just take that value and save some time
        if gene in analyzed:
            if analyzed[gene][1]:
                logger.debug("Gene %s was previously found to be multimodal" % gene)
                yield gene

            else:
                logger.debug("Gene %s was previously removed" % gene)
                continue

        # Center the expression data. This data should not be
        # centered before this step
        data = expression.loc[gene, :].values
        assert len(data) > 1, 'Dataset needs at least 2 observations!'

        mean = np.mean(data)
        data = data - mean

        # Skip Genes that have a mean below the mean filter
        if mean < gene_mean_filter:
            logger.debug("Gene %s was removed by mean filter." % gene)
            analyzed[gene] = (gene, False)
            continue

        # Reshape the data so that it is compatible with the
        # mixture model
        data = data.reshape(len(data), 1)

        # Convert data to a bnpy XData object
        X = bnpy.data.XData(data)
        logger.debug("Gene: %s\nnDim %d\n Size %d\n" % (gene,
                                                        X.dim,
                                                        X.get_size()))

        res = pool.apply_async(get_gene_fit,
                               args=(gene,
                                     X,
                                     min_prob_filter,
                                     output_dir,
                                     sensitive,
                                     mm_path,
                                     overwrite,))
        results.append(res)

    output = [r.get() for r in results]
    pool.close()
    pool.join()

    for gene, multimodal in output:
        if multimodal:
            yield gene


def filter_gene_sets(matrx, genesets, args):
    """
    Takes gene sets and filters down to the multimodally expressed genes

    :param genesets (dict): Dictionary containing gene sets
    :param args (argparse.Namespace): Input parameters
    :return: Filtered gene sets
    :rtype: collections.defaultdict
    """
    logger = logging.getLogger('root')
    for gs, genes in genesets.items():
        logger.debug(gs)
        start = len(genes)
        res = filter_gene_set(list(genes),
                              matrx,
                              CPU=args.CPU,
                              gene_mean_filter=args.min_mean_filter,
                              min_prob_filter=args.min_prob_filter,
                              output_dir=args.output_dir,
                              sensitive=args.sensitive,
                              mm_path=args.mm_path,
                              overwrite=args.overwrite)
        res = list(set(res))
        end = len(res)

        logger.info("Filtering: {gs} went from {x} to {y} genes".format(gs=gs,
                                                                        x=start,
                                                                        y=end))
        if end < args.min_num_filter:
            logger.info("Skipping {gs} because there "
                        "are not enough genes to cluster".format(gs=gs))
            continue

        yield gs, res


def filter(matrx, args):
    """
    Runs the multimodal filter and stops

    :param args argparse.Namespace: Pipeline configuration namespace
    :return: None
    """
    logger = logging.getLogger('root')
    logger.info("Running multimodal filter")

    if args.mm_path:
        logger.info("Loading previous fits...")
        rec = os.path.join(os.path.dirname(args.mm_path), 'ANALYZED')
        with open(rec) as f:
            for line in f:
                gene, result = line.strip().split('\t')
                analyzed[gene] = (gene, bool(int(result)))

    genes = matrx.index.values
    res = filter_gene_set(list(genes),
                          matrx,
                          CPU=args.CPU,
                          gene_mean_filter=args.min_mean_filter,
                          min_prob_filter=args.min_prob_filter,
                          output_dir=args.output_dir,
                          sensitive=args.sensitive,
                          mm_path=args.mm_path,
                          overwrite=args.overwrite)
    res = list(set(res))
    logger.info("Found %d multimodally distributed genes" % len(res))
    with open(os.path.join(args.output_dir, 'ANALYZED'), 'w') as f:
        for key, (_, value) in analyzed.items():
            f.write("%s\t%d\n" % (key, int(value)))

    if args.mm_path:
        mm_genes = args.mm_path
    else:
        mm_genes = os.path.join(args.output_dir, 'MultiModalGenes')
    assert os.path.exists(mm_genes), 'No output generated!'
    logger.info('Multimodal genes are located here:\n%s' % mm_genes)

