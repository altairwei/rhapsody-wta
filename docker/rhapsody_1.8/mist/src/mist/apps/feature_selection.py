from scipy.stats import ttest_ind
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from mist.apps.utils import *
from mist.apps.evaluation import *
import numpy as np
import matplotlib.pyplot as plt
import itertools
import rpy2.robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
import csv

# Some functions for scoring how good a split is
def sigclust(X,y):
    if len(np.unique(y)) != 2:
        print('ERROR: can only have 2 unique labels')
        return
    y = str_labels_to_ints(y)
    ro.r('library(sigclust)')
    ro.r('''g<-function(X,y){
                p = sigclust(X,100,labflag=1,label=y,icovest=3)
                return(p@pval)
            }''')
    g = ro.globalenv['g']
    p = g(numpy2ri(X),numpy2ri(y))
    return np.nan_to_num(-np.log10(p[0]))

# Feature selection using random forest
def skRandomForest(X,Y):
    clf = RandomForestClassifier()
    clf.fit(X,Y)
    return clf

def select_genes_using_RF(X,Y,return_score=False,verbose=False,return_comparisons=False,score_function=None):
    # find most important features using random forest
    n = 100 # number of times to repeat random forest
    m = int(np.round(X.shape[1]/2)) # number of features to output at the end
    if m > n: m = n
    # Save the number of times each gene appears in the top n most important features
    gene_scores = np.zeros(X.shape[1])
    # For each fitting, score how well the fitting classifies the data
    if return_score and score_function is None: scores = np.zeros(n)
    for i in range(0,n):
        clf = skRandomForest(X,Y)
        importances = clf.feature_importances_
        y_clf = clf.predict(X)
        if return_score and score_function is None: scores[i] = NMI(Y,y_clf)
        gene_ranks = np.flipud(np.argsort(importances))
        gene_scores[gene_ranks[0:m]] += 1
    # Rank genes by how often it appeared in the top n most important features
    gene_ranks = np.flipud(np.transpose(np.argsort(np.transpose(gene_scores))))[0:m]
    gene_scores /= float(n)
    outputs = [gene_ranks,gene_scores[gene_ranks]]
    if verbose: print('Mean prediction score: %.3f'%(np.mean(scores)))
    if return_score:
        if score_function is None: outputs.append(np.mean(scores))
        else: outputs.append(score_function(X,Y))
    if return_comparisons: outputs.append(compare_feature_means(X,Y)[gene_ranks])
    return tuple(outputs)

def log_select_genes_using_RF(X,Y,**kwargs):
    return select_genes_using_RF(np.log(1+X),Y,**kwargs)

def log_select_genes_using_RF_sigclust(X,Y,**kwargs):
    return select_genes_using_RF(np.log(1+X),Y,score_function=sigclust,**kwargs)

# Feature selection using L1-regularized logistic regression
def skLogisticRegression(X,Y):
    lr = LogisticRegression(penalty='l1')
    lr.fit(X,Y)
    return lr

def select_genes_using_LR(X,Y,return_score=False,verbose=False,return_comparisons=False,score_function=None):
    clf = skLogisticRegression(X,Y)
    y_clf = clf.predict(X)
    gene_scores = np.abs(clf.coef_)[0]
    gene_ranks = np.flipud(np.argsort(gene_scores))
    outputs = [gene_ranks,gene_scores[gene_ranks]]
    # if verbose: print 'Score: %.2f'%(score)
    if return_score:
        if score_function is None: score = compute_clustering_accuracy(Y,y_clf)
        else: score = score_function(X,Y)
        outputs.append(score)
    if return_comparisons: outputs.append(compare_feature_means(X,Y)[gene_ranks])
    return tuple(outputs)

def log_select_genes_using_LR(X,Y,**kwargs):
    return select_genes_using_LR(np.log(1+X),Y,**kwargs)

def log_select_genes_using_LR_sigclust(X,Y,**kwargs):
    return select_genes_using_LR(np.log(1+X),Y,score_function=sigclust,**kwargs)

# Feature selection using Welch's t-test
def select_genes_using_Welchs(X,Y,return_score=False,verbose=False,return_comparisons=False,score_function=None):
    if len(np.unique(Y)) > 2:
        print('ERROR: Y should only have 2 unique values')
    N,M = np.shape(X)
    y = str_labels_to_ints(Y)
    t,p = ttest_ind(X[y==1,:],X[y==0,:],equal_var=False)
    gene_inds = np.array(list(range(M)))
    finite_inds = np.isfinite(t)
    t,p,gene_inds = t[finite_inds],p[finite_inds],gene_inds[finite_inds]
    # Argsort first by p-value, then by -np.abs(t-value)
    keep_inds = np.array([j[2] for j in sorted([(p[i],-np.abs(t[i]),i) for i in range(len(t))])])
    if len(gene_inds) == 0 and len(keep_inds) == 0:
        outputs = [[], []]
        if return_score:
            if score_function is None: score = 0
            else: score = score_function(X,Y)
            outputs.append(score)
        if return_comparisons: outputs.append([])
    else:
        gene_ranks = gene_inds[keep_inds]
        gene_scores = np.nan_to_num(-np.log10(p[keep_inds]))
        outputs = [gene_ranks,gene_scores]
        # if verbose: print 'Score: %.2f'%(score)
        if return_score:
            if score_function is None: score = gene_scores[0]
            else: score = score_function(X,Y)
            outputs.append(score)
        if return_comparisons: outputs.append((t>0).astype(int)[keep_inds])
        # Note: For dendrosplit, each split involves exactly two clusters. 'L' will be mapped
        #   to 0, and 'R' will mapped to 1. A positive t statistic indicates that the 1-clust
        #   has a higher median than the 0-clust (i.e. 'L' has less of the feature than 'R' on
        #   average). Therefore return_comparisons will index the features that are more greater
        #   in 1-clust (or 'R'-clust). In the visualize_history function, 0's are red, 1's are
        #   green.
    return tuple(outputs)

def log_select_genes_using_Welchs(X,Y,**kwargs):
    return select_genes_using_Welchs(np.log(1+X),Y,**kwargs)

def log_select_genes_using_Welchs_sigclust(X,Y,**kwargs):
    return select_genes_using_Welchs(np.log(1+X),Y,score_function=sigclust,**kwargs)

# Visualization scripts
def one_from_rest_gene_selection(X,Y,feature_selector=log_select_genes_using_Welchs):
    selected_genes_for_each_cluster = []
    for ctype in np.unique(Y):
        selected_genes_for_each_cluster.append(feature_selector(X,Y == ctype,return_score=True))
    return selected_genes_for_each_cluster

def one_from_rest_visualize_genes(X,genes,x1,x2,y,num_genes=3,
                                  feature_selector=log_select_genes_using_Welchs):
    selected_genes_for_each_cluster = one_from_rest_gene_selection(X,y,feature_selector)
    unique_clusts = np.unique(y)
    for k,selected_genes in enumerate(selected_genes_for_each_cluster):
        plt.figure(figsize=(4*(num_genes+1),3))
        plt.subplot(1,num_genes+1,1)
        plt.scatter(x1,x2,c=y==unique_clusts[k],edgecolors='none')
        plt.title('Clust: %d, Acc: '%(unique_clusts[k])+sn(selected_genes[2]))
        _ = plt.axis('off')
        for i in range(num_genes):
            g = genes[selected_genes[0][i]]
            plt.subplot(1,num_genes+1,i+2)
            plt.scatter(x1,x2,c=X[:,genes == g],edgecolors='none')
            plt.title(g+' '+sn(selected_genes[1][i]))
            _ = plt.axis('off')

# Note that if you use save_name while show_plots is True, all the plots will be saved as well
def save_more_highly_expressed_genes_in_one_clust(X, genes, y, x1=None, x2=None, num_genes=3,
                                                  output_header=None, verbose=True,
                                                  feature_selector=log_select_genes_using_Welchs, save_name=None,
                                                  show_plots=True, pval_cutoff=10, skip_singleton=True):
    if show_plots:
        if x1 is None or x2 is None:
            print('NEED TO PASS IN x1, x2 FOR PLOTTING')

    if save_name is not None:
        with open(save_name+'_Cluster_Features.csv', 'w') as f:
            rb = csv.writer(f)
            for row in output_header:
                rb.writerow(row)
            rb.writerow(['Cluster', 'Gene', 'p-Value', 'Mean_of_Expression', 'Fold_Change_of_Expression'])

    for c in np.unique(y):

        if c == -1 and skip_singleton:
            continue

        c1 = X[y == c, :]
        c2 = X[y != c, :]
        keep_genes = np.mean(c1, 0) > np.mean(c2, 0)
        g_temp = genes[keep_genes]

        if len(g_temp) == 0:
            continue

        if verbose:
            print('Cluster %d: %d/%d genes kept'%(c, np.sum(keep_genes), len(keep_genes)))

        gene_ranks, gene_scores, score = feature_selector(X[:, keep_genes], y == c, return_score=True)
        num_genes_c = np.min([len(gene_ranks), num_genes])

        if show_plots and x1 is not None and x2 is not None:
            plt.figure(figsize=(4*(num_genes_c+1), 3))
            plt.subplot(1, num_genes_c+1, 1)
            plt.scatter(x1, x2, c=y == c, edgecolors='none')
            plt.title('Clust: %d, Score: ' % c + sn(score))
            _ = plt.axis('off')

        for i in range(num_genes_c):
            g = g_temp[gene_ranks[i]]
            g_ind = np.where(genes == g_temp[gene_ranks[i]])[0]

            if show_plots and x1 is not None and x2 is not None:
                plt.subplot(1,num_genes_c+1, i+2)
                plt.scatter(x1, x2, c=X[:, genes == g], edgecolors='none')
                plt.title(g+' '+sn(gene_scores[i]))
                _ = plt.axis('off')

            if save_name is not None and gene_scores[i] > pval_cutoff:
                g_mean_in = np.mean(X[y == c, g_ind])
                g_mean_out = np.mean(X[y != c, g_ind])
                with open(save_name + '_Cluster_Features.csv', 'a') as f:
                    rb = csv.writer(f)
                    rb.writerow([str(c), g_temp[gene_ranks[i]], sn(gene_scores[i], j=3), sn(g_mean_in, j=3),
                              sn(g_mean_in/g_mean_out, j=3)])

        if show_plots and x1 is not None and x2 is not None and save_name is not None:
            plt.savefig(save_name+'_cluster'+str(c)+'.png', format='png', dpi=300)


# Note that if you use save_name while show_plots is True, all the plots will be saved as well
def pairwise_cluster_comparison(X, genes, y, x1=None, x2=None, num_genes=20, output_header=None, verbose=True,
                                feature_selector=log_select_genes_using_Welchs, save_name=None,
                                show_plots=True, pval_cutoff=10, skip_singleton=True):
    if show_plots:
        if x1 is None or x2 is None:
            print('NEED TO PASS IN x1, x2 FOR PLOTTING')

    if save_name is not None:
        with open(save_name+'_Pairwise_Cluster_Features.csv', 'w') as f:
            rb = csv.writer(f)
            for row in output_header:
                rb.writerow(row)
            rb.writerow(['Comparison', 'Gene', 'p-Value', 'Larger_Cluster',
                         'Fold_Change_of_Expression_for_Larger_Cluster'])

    for (i,j) in itertools.combinations(np.unique(y), 2):
        if i == -1 and skip_singleton:
                continue
        labels_ij = np.logical_or(y==i,y==j)
        gene_ranks,gene_scores,score = feature_selector(X[labels_ij,:],y[labels_ij],return_score=True)

        if verbose:
            print('\nComparing cluster %s and %s'%(i, j))

        if show_plots and x1 is not None and x2 is not None:
            plt.figure(figsize=(4*(num_genes+1),3))
            plt.subplot(1,num_genes+1,1)
            plot_labels_legend(x1[labels_ij],x2[labels_ij],y[labels_ij],labels=[str(i),str(j)])
            plt.title('Cluster %s vs. Cluster %s, Score: '%(i, j)+sn(score))
            _ = plt.axis('off')
        num_effective_genes = min(num_genes, len(gene_ranks))
        for k in range(num_effective_genes):
            g_mean_i = np.mean(X[y == i, gene_ranks[k]])
            g_mean_j = np.mean(X[y == j, gene_ranks[k]])

            if g_mean_i > g_mean_j:
                c = i
                fold = g_mean_i/g_mean_j

            else:
                c = j
                fold = g_mean_j/g_mean_i

            if verbose:
                print('%s (more expressed in %s) '%(genes[gene_ranks[k]],c)+sn(gene_scores[k]))

            if show_plots and x1 is not None and x2 is not None:
                plt.subplot(1, num_genes+1, k+2)
                plt.scatter(x1[labels_ij], x2[labels_ij], c=X[labels_ij, gene_ranks[k]], edgecolors='none')
                plt.title(genes[gene_ranks[k]]+' '+sn(gene_scores[k]))
                _ = plt.axis('off')

                if save_name is not None:
                    plt.savefig()

            if save_name is not None and gene_scores[k] > pval_cutoff:
                with open(save_name+'_Pairwise_Cluster_Features.csv', 'a') as f:
                    rb = csv.writer(f)
                    rb.writerow(['Cluster'+str(i)+'_vs_Cluster'+str(j), genes[gene_ranks[k]], sn(gene_scores[k], j=3),
                                 str(c), sn(fold, j=3)])

        if show_plots and x1 is not None and x2 is not None and save_name is not None:
            plt.savefig(save_name+'_Cluster'+str(i)+'_vs_Cluster'+str(j)+'.png', format='png', dpi=300)
