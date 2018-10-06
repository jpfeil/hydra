import nbformat as nbf
import os


def create_notebook(src, gs, output_dir):
    """
    Creates a jupyter notebook for analyzing the output of
    hydra.

    :param name:
    :param output_path:
    :return:
    """

    nb = nbf.v4.new_notebook()

    text = """\
Code for investigating gene and sample clusters"""

    import_statements = """\
%matplotlib inline
import bnpy
import collections
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster, cophenet, linkage, dendrogram
from scipy.spatial.distance import pdist;"""

    load_model = """\
gene_set = "{gs}"
pth = os.path.abspath(".")
   
hmodel = bnpy.ioutil.ModelReader.load_model_at_prefix(pth,
                                                      prefix="{gs}")
                                                          
means = []
for comp in range(len(hmodel.allocModel.get_active_comp_probs())):
    m = hmodel.obsModel.get_mean_for_comp(comp)
    means.append(m)
    
if len(means) == 1:
    raise ValueError("Only one component was identified!")
    
pth = os.path.join(os.path.abspath('.'), 'training-data.tsv')
exp = pd.read_csv(pth, sep='\\t', index_col=0)
exp.head()
       
mean_df  = pd.DataFrame(data=np.vstack(means).T,
                        index=exp.index.values, 
                        columns=['cluster_%d' % c for c in range(len(means))])
                        
method = 'ward'
metric = 'euclidean'

zscore_df = mean_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1) 

row_linkage = hierarchy.linkage(
    distance.pdist(zscore_df.values), 
    method=method, metric=metric)

col_linkage = hierarchy.linkage(
    distance.pdist(zscore_df.values.T), 
    method=method, metric=metric);""".format(gs=gs,
                                             output=os.getcwd())

    functions = """\
def fancy_dendrogram(*args, **kwargs):
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
    return ddata;
"""

    row_dendro = """\
dist = 15

fancy_dendrogram(
    row_linkage,
    truncate_mode='lastp',
    p=12,
    leaf_rotation=90.,
    leaf_font_size=12.,
    show_contracted=True,
    annotate_above=10,
    max_d=dist,
)
plt.show()

rclusters = fcluster(row_linkage, dist, criterion='distance')

rcmap = sns.color_palette("Set2", max(rclusters))

rcolors = [rcmap[i-1] for i in rclusters]

row_groups = collections.defaultdict(list)
for gene, cluster in zip(zscore_df.index, rclusters):
    row_groups[cluster].append(gene);"""

    cluster_map = """\
cmap = sns.diverging_palette(240, 10, n=9)

sns.clustermap(zscore_df,
               col_linkage=col_linkage,
               row_linkage=row_linkage,
               row_colors=rcolors,
               cmap=cmap);"""

    print_genes = """\
print "\\n".join(row_groups[1]);"""

    fgsea = """\
import subprocess
from scipy.stats import ttest_ind

pth = <PATH TO BACKGROUND EXPRESSION>
background = pd.read_csv(pth, sep='\\t', index_col=0)

assign = pd.read_csv('assignments.tsv', 
                     sep='\
                     \t', 
                     index_col=0, 
                     header=None)

cmd = ["/usr/bin/Rscript",
       "%s/bin/fgsea.R",
       <PATH TO GENE SET FILE (.gmt)>,
       "/tmp/fgsea-analysis.rnk",
       "/tmp/fgsea-analysis.fgsea"]

fgseas = {}
for cluster, rows in assign.groupby(1):
    ins = rows.index.values
    outs = [x for x in background.columns if x not in ins]
    
    res = ttest_ind(background[ins].values,
                    background[outs].values,
                    axis=1).statistic
                    
    tstats = pd.DataFrame(index=background.index, 
                          data=res).dropna()
                          
    tstats = tstats.sort_values(0, ascending=False).reset_index()
    
    tstats.to_csv('/tmp/fgsea-analysis.rnk',
                  header=None,
                  sep='\\t',
                  index=False)
                  
    subprocess.check_call(cmd)
    
    fgsea = pd.read_csv('/tmp/fgsea-analysis.fgsea')
    
    fgseas[cluster] = fgsea
    
    os.remove('/tmp/fgsea-analysis.rnk')
    os.remove('/tmp/fgsea-analysis.fgsea')""" % src

    nb['cells'] = [nbf.v4.new_markdown_cell(text)] + \
                  [nbf.v4.new_code_cell(x) for x in [import_statements,
                                                     load_model,
                                                     functions,
                                                     row_dendro,
                                                     cluster_map,
                                                     print_genes,
                                                     fgsea]]

    pth = os.path.join(output_dir, '%s.ipynb' % gs)
    nbf.write(nb, pth)

    return pth
