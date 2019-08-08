#!/usr/bin/env python2.7
import bnpy
import collections
import itertools
import logging
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import pandas as pd
import seaborn as sns
import subprocess
import tempfile
import uuid

from scipy.stats import ttest_ind
from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster, dendrogram

from library.fit import run

src = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

y_genes = set()
pth = os.path.join(src, 'data/chrY-genes.txt')
with open(pth) as f:
    for line in f:
        y_genes.add(line.strip())


class EnrichmentAnalysis(object):
    def __init__(self,
                 mm_path=None,
                 exp_path=None,
                 min_prob_filter=0.2,
                 min_effect_filter=None,
                 gmt_path=None):
        """
        Performs hydra enrichment analysis.

        :param mm_path (str): Path to MultiModelGenes directory
        :param exp_path (str): Path to expression matrix
        :param gmt_path (str): Path to gene set database in GMT format
        :param min_prob_filter (float): Minimum cluster probability
        """

        if mm_path is None:
            raise ValueError('Need path to MultiModalGenes!')

        self.mm_path = mm_path

        if exp_path is None:
            raise ValueError('Need path to original expression matrix!')

        self.exp = pd.read_csv(exp_path, sep='\t', index_col=0)
        self.num_genes = self.exp.shape[0]
        self.num_samples = self.exp.shape[1]
        self.background = list(self.exp.index.values)
        self.gmt = gmt_path
        self.min_comp_filter = min_prob_filter
        self.min_effect_filter = min_effect_filter

        self.logger = logging.getLogger('root')
        self.logger.info("Startng hydra enrichment clustering analysis")
        self.logger.debug("Background: %d" % len(self.background))

        self.mm_genes = self.get_multimodal_genes()
        self.enrich = self.run_enrich()

    def get_multimodal_genes(self):
        """
        Reads in hydra gene-level fits. Filters genes
        that have low impact on population based on the
        minimum component probability parameter.

        :return:
        """
        assert os.path.exists(self.mm_path), "Can't find path to MultiModalGenes"
        genes = os.listdir(self.mm_path)
        mm = set()
        self.min_probs = []
        for gene in genes:
            if gene in y_genes:
                self.logger.debug('Skipping chrY genes: %s' % gene)
                continue

            model_pth = os.path.join(self.mm_path, gene)
            model = bnpy.ioutil.ModelReader.load_model_at_prefix(model_pth,
                                                                 prefix=gene)
            probs = model.allocModel.get_active_comp_probs()
            self.min_probs.append(min(probs))
            if min(probs) < self.min_comp_filter:
                self.logger.debug('Minimum Component Probability Filter: %s' % gene)
                continue

            effects = []
            for i, j in itertools.combinations(range(len(probs)), 2):
                if i == j:
                    continue

                mi = model.obsModel.get_mean_for_comp(i)
                mj = model.obsModel.get_mean_for_comp(j)
                effects.append(abs(mj - mi))

            if min(effects) < self.min_effect_filter:
                self.logger.debug('Minimum Effect Filter: %s' % gene)
                continue

            mm.add(gene)
        self.logger.info('Found %d multimodal genes' % len(mm))
        return list(mm)

    def run_enrich(self):
        """
        Runs clusterProfiler enrichment analysis

        :return: Enriched gene sets
        :rtype: pandas.DataFrame
        """
        mm_temp = os.path.join(tempfile.gettempdir(), 'MultiModal' + str(uuid.uuid4()))
        b_temp = os.path.join(tempfile.gettempdir(), 'BackGround' + str(uuid.uuid4()))
        res = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))

        if len(self.mm_genes) == 0:
            raise ValueError('No multimodal genes were identified.\n '
                             'Try using a lower minimum component probability threshold.')

        with open(mm_temp, 'w') as f, open(b_temp, 'w') as g:
            f.write('\n'.join(self.mm_genes))
            g.write('\n'.join(self.background))

        if os.path.exists(self.gmt):
            cmd = ['Rscript',
                   os.path.join(src, 'bin', 'enrich.R'),
                   mm_temp,
                   b_temp,
                   res,
                   self.gmt]

        elif self.gmt == 'GO':
            cmd = ['Rscript',
                   os.path.join(src, 'bin', 'go_enrich.R'),
                   mm_temp,
                   b_temp,
                   res]

        else:
            raise ValueError("Can't find GMT file:\n%s" % self.gmt)

        try:
            self.logger.debug(' '.join(cmd))
            stdout, stderr = run(cmd)

        except subprocess.CalledProcessError:
            self.logger.info(stdout)
            self.logger.info(stderr)
            raise

        except IOError:
            self.logger.info(stdout)
            self.logger.info(stderr)
            raise

        return pd.read_csv(res)

    def get_enriched_terms(self):
        """
        Helper function to return enriched term dataframe

        :return: Enriched Terms
        :rtype: pandas.DataFrame
        """
        return self.enrich

    def get_enriched_term_genes(self, regex=None):
        """
        Helper function to return all enriched go term genes

        :return:
        """
        if self.enrich.shape[0] == 0:
            raise ValueError("No enriched terms found.")

        t = self.enrich.copy()

        if regex:
            t = t[t['Description'].str.contains(regex)]

        genes = set()
        for g in t['geneID'].values:
            genes.update(g.strip().split('/'))

        return list(genes)


class MultivariateMixtureModel(object):
    def __init__(self,
                 data,
                 gamma=5.0,
                 variance=2.0,
                 K=5,
                 center=False,
                 verbose=False):

        self.gamma = gamma
        self.variance = variance
        self.K = K
        self.center = center
        self.verbose = verbose

        self.og_data = data.copy()
        self.hmodel = None
        self.clusters = None
        self.cluster_features = None

        self.hmodel = self.fit()

        if self.og_data.shape[0] > self.og_data.shape[1]:
            print 'WARNING: Number of genes outnumbers samples. ' \
                  'Consider more stringent filtering.'

    def fit(self, name='MultivariateAnalysis', verbose=False):
        """

        :param name:
        :return:
        """
        # This is a pandas dataframe: genes x samples
        data = self.og_data
        if self.center:
            if self.verbose:
                print 'centering data'
            data = data.apply(lambda x: x - x.mean(), axis=1)

        data = data.T.values
        xdata = bnpy.data.XData(data)

        workdir = tempfile.mkdtemp(prefix="%s_" % name)
        output_dir = 'K={K}-gamma={G}-ECovMat={Cov}-moves=birth,merge,delete,shuffle/'.format(K=self.K,
                                                                                              G=self.gamma,
                                                                                              Cov=self.variance)
        output_path = os.path.join(workdir, output_dir)
        hmodel, info_dict = bnpy.run(xdata,
                                     'DPMixtureModel',
                                     'Gauss',
                                     'memoVB',
                                     nLap=1000,
                                     nTask=1,
                                     nBatch=1,
                                     gamma0=self.gamma,
                                     sF=self.variance,
                                     ECovMat='eye',
                                     K=self.K,
                                     initname='randexamplesbydist',
                                     moves='birth,merge,delete,shuffle',
                                     b_startLap=0,
                                     m_startLap=2,
                                     d_startLap=2,
                                     output_path=output_path,
                                     doWriteStdOut=verbose)
        self.hmodel = hmodel
        self.clusters = collections.defaultdict(list)
        for sample, cluster in zip(self.og_data.columns, self.get_assignments(self.og_data)):
            self.clusters[cluster].append(sample)
        return self.hmodel

    def get_assignments(self, data):
        # This is a new pandas dataframe that may
        # not be the same as the one we trained on

        genes = self.og_data.index.values
        data = data.reindex(genes).dropna()
        if self.center:
            data = data.sub(self.og_data.mean(axis=1), axis=0)

        if data.ndim == 1:
            data = data.values.reshape(1, data.shape[0])

        else:
            data = data.T.values

        xdata = bnpy.data.XData(data)

        unclass = 1 - np.sum(self.hmodel.allocModel.get_active_comp_probs())
        LP = self.hmodel.calc_local_params(xdata)
        asnmts = []
        for row in range(LP['resp'].shape[0]):
            _max = np.max(LP['resp'][row, :])
            if _max < unclass:
                asnmts.append(-1)

            else:
                _arg = np.argmax(LP['resp'][row, :])
                asnmts.append(_arg)
        return asnmts

    def get_cluster_features(self, exp, gmt=None):
        """
        Returns dictionary with enriched gene set terms per cluster

        :param exp (pd.DataFrame): Original expression matrix with all genes
        :param gmt (str): Path to GMT file. Gene annotation must match expression matrix
        :return:
        """
        if self.clusters is None:
            raise ValueError("Need to fit model first!")

        if len(exp) < 1000:
            raise ValueError('Please pass the original expression matrix before multimodal filtering')

        self.cluster_features = {}
        for c, samples in self.clusters.items():
            ins = samples
            outs = []
            for _c, _samples in self.clusters.items():
                if _c == c:
                    continue
                outs.extend(_samples)

            res = ttest_ind(exp[ins].values,
                            exp[outs].values,
                            axis=1).statistic

            tstats = pd.DataFrame(index=exp.index,
                                  data=res).dropna()

            tstats = tstats.sort_values(0,
                                        ascending=False)

            gsea = n1(tstats, gmt=gmt)
            self.cluster_features[c] = gsea.sort_values('NES', ascending=False)
        return self.cluster_features

    def sub_cluster_gsea(self, data, constant=0.05, gmt=None, return_diff=False, alpha=0.05, debug=False):
        """
        Performs N-of-1 GSEA by normalizing the sample to the cluster background. (BETA)

        :param data (pd.Series): Expression Series for N-of-1 sample
        :param constant (float): PAM normalizing factor for ranking genes
        :param gmt (str): Path to gene set file in GMT format
        :param return_diff (bool): Removes enriched gene sets that could have been identified at the cohort level
        :param alpha (float): P-value threshold for assessing statistical significance.
        :return:
        """
        a = self.get_assignments(data).pop()
        if a == -1:
            raise ValueError("Sample did not place in a cluster.")

        background = self.og_data
        data = data.reindex(background.index)
        bsamples = self.clusters[a]
        zscore = (data - background[bsamples].mean(axis=1)) / (background[bsamples].std(axis=1) + constant)
        zscore = zscore.sort_values(ascending=False)
        if debug:
            print(zscore.head())
        fgsea = n1(zscore, gmt)
        if return_diff:
            zscore2 = (data - background.mean(axis=1)) / (background.std(axis=1) + constant)
            fgsea2 = n1(zscore2, gmt)
            sig1 = fgsea[fgsea['padj'] < alpha]
            sig2 = fgsea2[fgsea2['padj'] < alpha]
            index = set(sig1.index.values) - set(sig2.index.values)
            fgsea = fgsea.reindex(index)
        return a, fgsea

    def cohort_gsea(self, data):
        raise NotImplementedError('Sorry! Still working on this')


class HClust(object):
    def __init__(self, data, method='ward', metric='euclidean'):
        self.og_data = data.copy()
        self.cluster(method=method,
                     metric=metric)
        self.row_groups = None
        self.col_groups = None

    def cluster(self, method='ward', metric='euclidean'):
        zscore = self.og_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

        self.row_linkage = hierarchy.linkage(
            distance.pdist(zscore.values),
            method=method, metric=metric)

        self.col_linkage = hierarchy.linkage(
             distance.pdist(zscore.values.T),
             method=method, metric=metric)

    def plot_row_linkage(self, dist):
        """

        :param dist: Distance for determining groups
        :return:
        """

        fancy_dendrogram(
            self.row_linkage,
            truncate_mode='lastp',
            p=12,
            leaf_rotation=90.,
            leaf_font_size=12.,
            show_contracted=True,
            annotate_above=10,
            max_d=dist,
        )
        plt.show()

    def get_row_groups(self, dist):
        """

        :param dist: Distance for determining groups
        :return:
        """
        clusters = fcluster(self.row_linkage,
                            dist,
                            criterion='distance')

        groups = collections.defaultdict(list)

        for gene, cluster in zip(self.og_data.index.values,
                                   clusters):
            groups[cluster].append(gene)

        self.row_groups = groups

        return groups

    def plot_col_linkage(self, dist):
        """

        :param dist: Distance for determining groups
        :return:
        """
        fancy_dendrogram(
            self.col_linkage,
            truncate_mode='lastp',
            p=12,
            leaf_rotation=90.,
            leaf_font_size=12.,
            show_contracted=True,
            annotate_above=10,
            max_d=dist,
        )
        plt.show()

    def get_col_groups(self, dist):
        """

        :param dist: Distance for determining groups
        :return:
        """
        clusters = fcluster(self.col_linkage,
                            dist,
                            criterion='distance')
        groups = collections.defaultdict(list)
        for sample, cluster in zip(self.og_data.columns,
                                   clusters):
            groups[cluster].append(sample)
        self.col_groups = groups
        return groups

    def plot(self):
        return sns.clustermap(self.og_data,
                              z_score=0,
                              col_linkage=self.col_linkage,
                              row_linkage=self.row_linkage,
                              method='ward',
                              center=0,
                              cmap=sns.diverging_palette(240, 10, n=7),
                              figsize=(10, 10))


class ScanEnrichmentAnalysis(object):
    def __init__(self, mm_path, exp_path, gmt_path, min_prob_range=None, min_effect_filter=1.0, **kwargs):
        """
        Class to explore the function of multimodally expressed genes.

        :param min_prob_range (iterable): Iterable containing floats between 0 and 0.5
        """
        self.mm_path = mm_path
        self.exp_path = exp_path
        self.gmt_path = gmt_path
        self.K = kwargs['K'] if 'K' in kwargs else 1
        self.variance = kwargs['variance'] if 'variance' in kwargs else 2.0
        self.gamma = kwargs['gamma'] if 'gamma' in kwargs else 5.0
        self.CPU = kwargs['CPU'] if 'CPU' in kwargs else 1
        self.cluster = kwargs['cluster'] if 'cluster' in kwargs else True

        if min_prob_range is None:
            min_prob_range = [round(x, 2) for x in np.linspace(0.1, 0.4, 10)]
        self.results = pd.DataFrame(index=min_prob_range,
                                    columns=['num_genesets', 'gs_terms', 'gs_term_genes',
                                             'num_genes', 'num_clusters', 'num_samples'])
        self.results.index.name = 'min_prob_filter'
        self.min_prob_range = min_prob_range
        self.min_effect_filter = min_effect_filter

    def scan(self):
        pool = multiprocessing.Pool(self.CPU)
        results = []
        for min_prob in self.min_prob_range:
            args = (self.mm_path,
                    self.exp_path,
                    self.gmt_path,
                    min_prob,
                    self.min_effect_filter,
                    self.cluster,
                    self.gamma,
                    self.variance,
                    self.K,)
            res = pool.apply_async(_get_enrichment_analysis, args=args)
            results.append(res)

        results = [x.get() for x in results]

        for res in results:
            min_prob = res[0]
            self.results.loc[min_prob, :] = res[1:]
        self.plot()
        return self.results

    def plot(self):
        fig, ax = plt.subplots(1,
                               figsize=(6, 4))

        _scan = self.results.reset_index()
        _scan.plot(x='min_prob_filter',
                   y='num_genes',
                   ax=ax,
                   legend=False)
        ax2 = ax.twinx()
        _scan.plot(x='min_prob_filter',
                   y='num_clusters',
                   ax=ax2,
                   color='black',
                   legend=False)
        ax.axhline(_scan['num_samples'].unique()[0], color='red')
        ax.set_ylabel("Number of Samples/Number of Genes")
        ax2.set_ylabel("Number of Clusters")
        ax.set_xlabel("Minimum Component Probability Threshold")
        ax.legend(['#Genes', '#Samples'], loc=1)
        ax2.legend(["#Clusters"], loc=7)
        plt.show()


class SweepAnalysis(object):
    def __init__(self, path):
        """

        :param path (str): Path to MultivariateAnalysis directory
        """

        self.path = path

    def rank(self):
        hits = pd.DataFrame(columns=['gene-set', 'num_clusters', 'maxKL'])

def fancy_dendrogram(*args, **kwargs):
    """
    Code was adapted from:
    https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/

    :param args:
    :param kwargs:
    :return:
    """
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


def n1(zscore, gmt=None):
    """
    N-of-1 GSEA. GMT annotation must matched expression annotation

    :param zscore: Ranked dataframe
    :param gmt: Path to GSEA GMT file
    :return: GSEA Dataframe
    """
    rnk_temp = os.path.join(tempfile.gettempdir(), 'RNK' + str(uuid.uuid4()))
    fgsea_temp = os.path.join(tempfile.gettempdir(), 'FGSEA' + str(uuid.uuid4()))

    if gmt is None:
        gmt = os.path.join(src,
                           'gene-sets',
                           'Human_GOBP_AllPathways_no_GO_iea_December_01_2018_symbol.gmt')
    cmd = ['Rscript',
           os.path.join(src, 'bin', 'fgsea.R'),
           gmt,
           rnk_temp,
           fgsea_temp]
    zscore.to_csv(rnk_temp,
                  header=False,
                  sep='\t')
    subprocess.check_call(cmd)

    return pd.read_csv(fgsea_temp, index_col=0)


def _get_enrichment_analysis(mm_path, exp_path, gmt_path, min_prob_filter, min_effect_filter,
                             cluster=True, gamma=5.0, variance=2.0, K=1, center=True):

    res = EnrichmentAnalysis(mm_path=mm_path,
                             exp_path=exp_path,
                             gmt_path=gmt_path,
                             min_prob_filter=min_prob_filter,
                             min_effect_filter=min_effect_filter)

    try:
        enriched_terms = res.get_enriched_terms()
        terms = enriched_terms['Description'].values
        enriched_term_genes = res.get_enriched_term_genes()

        if cluster:
            exp = pd.read_csv(exp_path, sep='\t', index_col=0)
            clus = MultivariateMixtureModel(exp.reindex(enriched_term_genes),
                                            gamma,
                                            variance,
                                            K,
                                            center)
            num_clusters = len(clus.clusters)
        else:
            num_clusters = np.nan
        return [res.min_comp_filter,
                len(enriched_terms),
                '|'.join(terms),
                '|'.join(enriched_term_genes),
                len(res.get_enriched_term_genes()),
                num_clusters,
                res.num_samples]

    except ValueError:
        return [res.min_comp_filter,
                0,
                np.nan,
                np.nan,
                0,
                np.nan,
                res.num_samples]


def kl(m1, S1, m2, S2):
    """
    https://stats.stackexchange.com/questions/60680/kl-divergence-between-two-multivariate-gaussians

    :param m1:
    :param S1:
    :param m2:
    :param S2:
    :return:
    """
    m1 = m1.reshape(len(m1), 1)
    m2 = m2.reshape(len(m2), 1)
    detS2 = np.linalg.det(S2)
    detS1 = np.linalg.det(S1)
    d = len(m1)
    invS2 = np.linalg.inv(S2)
    trS2S1 = np.trace(invS2 * S1)
    mah = np.matmul(np.matmul((m2 - m1).T, invS2), (m2 - m1))
    return float(0.5 * (np.log(detS2 / detS1) - d + trS2S1 + mah))