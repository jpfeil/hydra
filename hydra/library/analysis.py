#!/usr/bin/env python2.7
import bnpy
import collections
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import subprocess
import tempfile
import uuid

from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster, cophenet, linkage, dendrogram
from scipy.spatial.distance import pdist

from library.fit import subprocess_fit

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
                 min_comp_filter=0.2,
                 gmt=None):

        if mm_path is None:
            raise ValueError('Need path to MultiModalGenes\nRequires --save-genes setting')

        if exp_path is None:
            raise ValueError('Need path to original expression matrix')

        self.exp = pd.read_csv(exp_path, sep='\t', index_col=0)
        self.background = list(self.exp.index.values)

        self.mm_path = mm_path
        self.min_comp_filter = min_comp_filter

        self.mm_genes = self.get_multimodal_genes()

        self.gmt = gmt
        self.enrich = self.run_enrich()

    def get_multimodal_genes(self, verbose=False):
        assert os.path.exists(self.mm_path), "Can't find path to MultiModalGenes"
        genes = os.listdir(self.mm_path)
        mm = set()
        self.min_probs = []
        for gene in genes:
            if gene in y_genes:
                if verbose:
                    print 'Skipping chrY genes: ', gene
                continue

            model_pth = os.path.join(self.mm_path, gene)
            model = bnpy.ioutil.ModelReader.load_model_at_prefix(model_pth,
                                                                 prefix=gene)
            probs = model.allocModel.get_active_comp_probs()
            self.min_probs.append(min(probs))
            if min(probs) < self.min_comp_filter:
                if verbose:
                    print 'Minimum Component Probability Filter: ', gene
                continue

            mm.add(gene)
        return list(mm)

    def _recommend_min_comp(self, variance=0.15):
        "In development. Not ready for use."

        if self.min_probs is None:
            raise ValueError('Need to init class first!')

        sns.distplot(self.min_probs)

        d = pd.DataFrame(index=range(len(self.min_probs)),
                         columns=[0])

        d.loc[:, 0] = self.min_probs
        m = MultivariateMixtureModel(d, center=True, variance=variance)
        a = m.get_assignments(d)
        groups = [[] for _ in m.hmodel.allocModel.get_active_comp_probs()]
        for v, c in zip(self.min_probs, a):
            groups[c].append(v)
        means = [np.mean(x) for x in groups]
        _max = np.argmax(means)
        return min(groups[_max])

    def run_enrich(self):
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
            subprocess.check_call(cmd)

        except subprocess.CalledProcessError:
            print ' '.join(cmd)
            raise

        return pd.read_csv(res)

    def get_enriched_terms(self):
        return self.enrich

    def get_enriched_term_genes(self):
        if self.enrich.shape[0] == 0:
            raise ValueError("No GO terms found.")

        genes = set()
        for g in self.enrich['geneID'].values:
            genes.update(g.strip().split('/'))

        return list(genes)


class MultivariateMixtureModel(object):
    def __init__(self,
                 data,
                 gamma = 5.0,
                 variance = 2.0,
                 K = 5,
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

        self.hmodel = self.fit()

        if self.og_data.shape[0] > self.og_data.shape[1]:
            print 'WARNING: Number of genes outnumbers samples. ' \
                  'Consider more stringent filtering.'

    def fit(self):
        # This is a pandas dataframe: genes x samples
        data = self.og_data
        if self.center:
            if self.verbose:
                print 'centering data'
            data = data.apply(lambda x: x - x.mean(), axis=1)

        data = data.T.values
        xdata = bnpy.data.XData(data)

        output = subprocess_fit('EnrichmentAnalysis',
                                xdata,
                                gamma=self.gamma,
                                sF=self.variance,
                                K=self.K)

        if output[1] is False:
            print 'WARNING: DPMM algorithm did not converge!'

        self.hmodel = output[0]
        self.clusters = collections.defaultdict(list)
        for sample, cluster in zip(self.og_data.columns, self.get_assignments(self.og_data)):
            self.clusters[cluster].append(sample)
        return self.hmodel

    def get_assignments(self, data):
        # This is a new pandas dataframe that may
        # not be the same as the one we trained on

        genes = self.og_data.index.values
        data = data.reindex(genes)
        if self.center:
            data = data.sub(self.og_data.mean(axis=1), axis=0)
        xdata = bnpy.data.XData(data.T.values)

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

    def n_of_1(self, data, constant=0.05):
        a = self.get_assignments(data)
        if a == -1:
            raise ValueError("Sample did not place in a cluster.")

        data = data.reindex(self.og_data.index)
        background = self.clusters[a]
        zscore = (data - self.og_data[background].mean(axis=1)) / (self.og_data[background].std(axis=1) + constant)
        return n1(zscore)


class HClust(object):
    def __init__(self, data, method='ward', metric='euclidean'):
        self.og_data = data.copy()
        self.cluster(method=method,
                     metric=metric)

    def cluster(self, method='ward', metric='euclidean'):
        zscore = self.og_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

        self.row_linkage = hierarchy.linkage(
            distance.pdist(zscore.values),
            method=method, metric=metric)

        self.col_linkage = hierarchy.linkage(
             distance.pdist(zscore.values.T),
             method=method, metric=metric);

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

        for sample, cluster in zip(self.og_data.columns,
                                   clusters):
            groups[cluster].append(sample)

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

class HydraUnsupervisedAnalysis(object):
    def __init__(self):
        raise NotImplementedError


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


def n1(zscore):
    rnk_temp = os.path.join(tempfile.gettempdir(), 'RNK' + str(uuid.uuid4()))
    fgsea_temp = os.path.join(tempfile.gettempdir(), 'FGSEA' + str(uuid.uuid4()))

    cmd = ['Rscript',
           os.path.join(src, 'bin', 'fgsea.R'),
           os.path.join(src, 'data', 'Human_GO_AllPathways_no_GO_iea_October_01_2018_symbol.gmt'),
           rnk_temp,
           fgsea_temp]

    zscore = zscore.sort_values(ascending=False)

    zscore.to_csv(rnk_temp,
                  header=None,
                  sep='\t')

    subprocess.check_call(cmd)

    return pd.read_csv(fgsea_temp, index_col=0)