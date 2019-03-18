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
                 min_comp_filter=0.2):

        if mm_path is None:
            raise ValueError('Need path to MultiModalGenes\nRequires --save-genes setting')

        if exp_path is None:
            raise ValueError('Need path to original expression matrix')

        self.exp = pd.read_csv(exp_path, sep='\t', index_col=0)
        self.background = list(self.exp.index.values)

        self.mm_path = mm_path
        self.min_comp_filter = min_comp_filter

        self.mm_genes = self.get_multimodal_genes()
        self.go = self.go_enrich()

    def get_multimodal_genes(self, verbose=False):
        assert os.path.exists(self.mm_path), "Can't find path to MultiModalGenes"
        genes = os.listdir(self.mm_path)
        mm = set()
        for gene in genes:
            if gene in y_genes:
                if verbose:
                    print 'Skipping chrY genes: ', gene
                continue

            model_pth = os.path.join(self.mm_path, gene)
            model = bnpy.ioutil.ModelReader.load_model_at_prefix(model_pth,
                                                                 prefix=gene)
            probs = model.allocModel.get_active_comp_probs()
            if min(probs) < self.min_comp_filter:
                if verbose:
                    print 'Minimum Component Probability Filter: ', gene
                continue

            mm.add(gene)
        return list(mm)

    def go_enrich(self):
        mm_temp = os.path.join(tempfile.gettempdir(), 'MultiModal' + str(uuid.uuid4()))
        b_temp = os.path.join(tempfile.gettempdir(), 'BackGround' + str(uuid.uuid4()))
        res = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))

        if len(self.mm_genes) == 0:
            raise ValueError('No multimodal genes were identified.\n '
                             'Try using a lower minimum component probability threshold.')

        with open(mm_temp, 'w') as f, open(b_temp, 'w') as g:
            f.write('\n'.join(self.mm_genes))
            g.write('\n'.join(self.background))

        cmd = ['Rscript',
               os.path.join(src, 'bin', 'go_enrich.R'),
               mm_temp,
               b_temp,
               res]

        subprocess.check_call(cmd)

        return pd.read_csv(res)

    def get_go_terms(self):
        return self.go

    def get_go_term_genes(self):
        if self.go.shape[0] == 0:
            raise ValueError("No GO terms found.")

        genes = set()
        for g in self.go['geneID'].values:
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

        pth = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        hmodel, info_dict = bnpy.run(xdata,
                                     'DPMixtureModel',
                                     'Gauss',
                                     'memoVB',
                                     output_path = pth,
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
                                     d_startLap=2)
        return hmodel


    def get_assignments(self, data):
        # This is a new pandas dataframe that may
        # not be the same as the one we trained on

        genes = self.og_data.index.values
        data = data.reindex(genes).sub(self.og_data.mean(axis=1), axis=0)
        xdata = bnpy.data.XData(data.T.values)

        unclass = 1 - np.sum(self.fit.allocModel.get_active_comp_probs())

        LP = self.fit.calc_local_params(xdata)
        asnmts = []
        for row in range(LP['resp'].shape[0]):
            _max = np.max(LP['resp'][row, :])
            if _max < unclass:
                asnmts.append(-1)

            else:
                _arg = np.argmax(LP['resp'][row, :])
                asnmts.append(_arg)
        return asnmts


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
