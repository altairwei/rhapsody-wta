from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
import matplotlib.colors
import time
import logging
from mist.apps.preprocessing import *
from mist.apps.feature_selection import *
from mist.apps.utils import *
import sys
from functools import reduce
sys.setrecursionlimit(10000)

# Using the linkage matrix Z, iteratively perform the following from the top:
#    1. Split into two clusters
#    2. If one of the clusters are too small (< min_clust_size), mark all points in it as outliers/error. remove it
#    3. Else, use random forest to classify the two types against each other. If high enough accuracy, keep clusters
#    4. Repeat on subclusters that aren't marked as outliers/error

# Given an inner node, split the two subtrees into two groups 
def node_to_leaf_labels(Z,i,ltypes=['0','1'],defaultlabel='_'):
    # Z: hierarchical clustering result returned by scipy linkage function
    # i: row of Z to perform splitting on 
    # ltypes: 2x1 vector. labels given to subtrees
    n = int(Z[-1,3])
    labels = [defaultlabel for _ in range(n)]
    def recursion(j,label):
        # base case 1: leaf node
        if j < n: 
            labels[j] = label
        # base case 2: parent of two leaf nodes
        elif Z[j-n,3] == 2:
            labels[int(Z[j-n,0])] = label
            labels[int(Z[j-n,1])] = label
        # recursive case: 
        else:
            recursion(int(Z[j-n,0]),label)
            recursion(int(Z[j-n,1]),label)
    recursion(int(Z[i,0]),ltypes[0])
    recursion(int(Z[i,1]),ltypes[1])
    return np.array(labels)

# Label updating algorithm for dendrosplit
def update_labels(y,labels,background='_'):
    label = None
    for i in range(len(y)):
        if y[i] != background:
            if label == None: label = labels[i]
            labels[i] += y[i] # Update labels 
    return label,labels

# Dendrosplit algorithm
def dendrosplit(X, preprocessing=log_correlation, min_clust_size=2, score_threshold=10, method='complete',
                verbose=False, split_evaluator=log_select_genes_using_Welchs, disband_percentile=100):
    '''
        Perform iterative hierarchical clustering. Obtain a dendrogram using a distance matrix. 
        Starting from the top, iteratively split the tree into two subtrees. Perform a test to 
        see if the two subtrees are different enough.
        
        Inputs:
            X: Matrix of counts 
            preprocessing: function to convert X to a distance matrix. This argument
                can be set to 'precomputed', in which case X should be a tuple (D,X)
                where D is the already-computed distance matrix
            min_clust_size: min size of cluster in order for a cluster split to be considered
            score_threshold: min score achieved by split required to keep the split
            method: linkage method for hierarchical clustering
            verbose: whether or not to turn on print statements while the algorithm runs
            split_evaluator: function to score how good a split is
            disband_percentile: if both subtrees have less than min_clust_size samples or if 
                a candidate split does not achieve high enough a score, look at the pairwise 
                distances amongst samples in this final cluster. if all are greater than 
                this percentile of the distances, then mark all of the points as singletons
            
        Outputs:
            y: cluster labels 
            history: list of tuples of information generated at each split. Use the functions
                'print_history' and 'visualize_history' with this output
    '''
    
    start_time = time.time()
    
    # Catch edge cases
    if preprocessing != 'precomputed': 
        if np.sum(np.sum(X,0) == 0) > 1:
            logging.info('Please remove genes that sum to 0 (try split.filter_genes())')
            return
#         if not check_if_all_entries_are_whole(X):
#             print 'Please make sure X is a valid distance matrix'
#             return
        D = preprocessing(X)

        if verbose: logging.info('Preprocessing took %.3f s'%(time.time()-start_time))

    elif type(X) != tuple:
        logging.info('preprocessing = precomputed. Feed in (D,X) tuple as first argument.')
        return

    elif not np.all(np.isclose(X[0],X[0].T, equal_nan=True)):
        logging.info('Need a valid distance matrix.')
        return

    else:
        D,X = X[0],X[1]
        if np.sum(np.sum(X,0) == 0) > 1:
            logging.info('Please remove columns that sum to 0 (try split.filter_genes())')
            return
#         if not check_if_all_entries_are_whole(X):
#             print 'Please make sure X is a valid distance matrix'
#             return

    # Perform the hierarchical clustering
#    D = (D+D.T)/2
    Ds = squareform(D)
    Z = linkage(Ds,method=method)
    N = len(X)

    # Compute disband percentile
    disband_threshold = np.percentile(flatten_distance_matrix(D),disband_percentile)
    
    # start with the label 'r' (for 'root') for all samples
    labels = ['r' for i in range(N)]

    # Stuff we may want to remember about each stage of splitting:
    #  - labels before the split and # of labels
    #  - labels after the split and # of each labels
    #  - splitting score
    #  - important features for distingushing the two subtypes
    #  - which cluster has greater expression of each feature
    history,scount = [],[0]
    
    def recursion(j,split_depth,labels):

        j = int(j)

        # Check number of leaves below node exceeds required threshold
        # Or if we're at a leaf node
        if j < 0: return
        if Z[j,3] < min_clust_size: return
        
        # Split subtrees and get labels
        y_node = node_to_leaf_labels(Z,j,ltypes=['L','R'])
      
        lsize = np.sum(y_node == 'L')
        rsize = np.sum(y_node == 'R')
        if verbose: logging.info('Potential split result: %d and %d'%(lsize,rsize))
        
        # If the split produces a cluster that's too small, ignore the smaller cluster
        if lsize < min_clust_size and rsize < min_clust_size: 
            # if both clusters are too small, check distances. If all distances are above
            # the disband_percentile overall distance, then make all samples singletons
            if np.sum(flatten_distance_matrix(D,y_node!='_') < disband_threshold) == 0:
                if verbose: logging.info('Disbanding (points in cluster too far from each other)')
                for i in range(N):
                    # guarantee that each sample has a unique label
                    if y_node[i] != '_': labels[i] += y_node[i]+str(i)
            return

        if lsize < min_clust_size: 
            label,labels = update_labels(y_node,labels)
            history.append(({label: np.sum(y_node != '_'),
                             label+'L': np.sum(y_node == 'L'),
                             label+'R': np.sum(y_node == 'R')},
                            y_node,Z[j,2],split_depth))
            recursion(Z[j,1]-N,split_depth,labels)
            return
        if rsize < min_clust_size: 
            label,labels = update_labels(y_node,labels)
            history.append(({label: np.sum(y_node != '_'),
                             label+'L':np.sum(y_node == 'L'),
                             label+'R':np.sum(y_node == 'R')},
                            y_node,Z[j,2],split_depth))
            recursion(Z[j,0]-N,split_depth,labels)
            return
        
        # Given a data matrix and some labels, fit a classifier. Evaluate classification accuracy
        scount[0] += 1
        genes,rankings,score,comparisons = split_evaluator(X[y_node!='_',:],y_node[y_node!='_'],
                                                           return_score=True,return_comparisons=True)
        #if verbose and type(score) is not list: print ' Split score '+sn(score)
        if verbose: logging.info(' Split score '+sn(score))

        l = str_labels_to_ints(y_node)
        intra1 =median_cdist_corr(D,0,0,l)
        intra2 =median_cdist_corr(D,1,1,l)
        inter = median_cdist_corr(D, 0, 1, l)
        weighted_diff = (np.sum(l==0)**2*intra1+np.sum(l==1)**2*intra2)/(np.sum(l==0)**2+np.sum(l==1)**2)-inter
        #print intra1, intra2, inter, weighted_diff
        # Decide if we should keep the split based on the score
        if score > score_threshold and (history == [] or weighted_diff < 0):
            label,labels = update_labels(y_node,labels)
            history.append(({label: np.sum(y_node != '_'),
                             label+'L': np.sum(y_node == 'L'),
                             label+'R': np.sum(y_node == 'R')},
                             y_node, Z[j,2], split_depth+1,
                             genes, rankings, score, comparisons))
            recursion(Z[j,0]-N,split_depth+1,labels)
            recursion(Z[j,1]-N,split_depth+1,labels)
        else: 
            if np.sum(flatten_distance_matrix(D,y_node!='_') < disband_threshold) == 0:
                if verbose: logging.info('Disbanding (points in cluster too far from each other)')
                for i in range(N):
                    if y_node[i] != '_': labels[i] += y_node[i]+str(i)
            return
        
    recursion(len(Z)-1,0,labels)
    if verbose: 
        logging.info('# of times score function was called: %d'%(scount[0]))
    logging.info('Total computational time was %.3f s'%(time.time()-start_time))

    return np.array(labels),history


# calculate full dendrogram
def plot_dendro(D, sample_names=None, p=None, labels=None, save_name=None, method='complete', lr=90,
                return_cell_order=False, fig_dim=(100,10), font_size=7):
    '''
    Calculate the full dendrogram based on a distance matrix D
    
    Inputs: 
    D: distance matrix (NxN)
    sample_names: name for each sample (uses the index by default, length-N vector)
    p: depth to cut the tree (goes to leaves by default)
    labels: a label for each sample (also a length-N vector, good for indicating cluster IDs)
    method: linkage methods for scipy's hierarchical clustering algorithm
    lr: leaf rotation
    return_cell_order: boolean. If true, returns the indices of the samples in the order of the dendrogram
    fig_dim: tuple indicating the size of the figure
    font_size: font size of leaves (7 by default)
    '''

    if p is None: p = len(D)
    if sample_names is None: sample_names = list(range(len(D)))
    if labels is not None:
        singleton_col = '#808080' # gray
        Ncolors = len(np.unique(labels))
        for i,c in enumerate(np.unique(labels)):
            if np.sum(labels == c) == 1:
                labels[labels == c] = -1
                Ncolors -= 1
        HSVs = [(x*1.0/Ncolors, 0.8, 0.9) for x in range(Ncolors)]
        RGBs = [colorsys.hsv_to_rgb(*x) for x in HSVs]
        color_map = {i:matplotlib.colors.rgb2hex(RGBs[j]) for j,i in enumerate(range(Ncolors))}
        color_map[-1] = singleton_col
        colors = [color_map[x] for x in labels]
    # Convert n-by-n D to (n choose 2)-by-1 vector for scipy's linkage function
    Ds = squareform(D)
    # 'complete' merges clusters based on farthest points in two clusters
    # 'single' merges clusters based on closest points in two clusters
    # 'average' merges clusters based on the average pairwise distance between all points in two clusters
    # 'weighted' merges clusters A and B based on the average distance of the two clusters that were merged to become A
    Z = linkage(Ds,method=method)
    plt.figure(figsize=fig_dim)
    link_cols = {}

    for i, i12 in enumerate(Z[:,:2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(Z) else color_map[labels[x]] for x in i12)
        if c1 == c2:
            link_cols[i+1+len(Z)] = c1
        else:
            link_cols[i+1+len(Z)] = singleton_col
    D = dendrogram(Z,color_threshold=None,p=p,leaf_rotation=lr,
                   leaf_font_size=font_size,labels=sample_names,link_color_func=lambda x: link_cols[x])
    ax = plt.gca()

    if labels is not None:
        xlbls = ax.get_xmajorticklabels()
        for i in range(len(xlbls)):
            xlbls[i].set_color(colors[int(float(xlbls[i].get_text()[0:]))])

    if save_name is not None:
        plt.savefig(save_name+'.png', format='png', dpi=300)

    if return_cell_order:
        return np.array([int(i.get_text()) for i in ax.get_xmajorticklabels()])

        
# Visualize how the splitting progresses
def visualize_history(X, x1, x2, genes, history, score_threshold=10, save_name=None):
    prefix = 'Split'
    ii = 1

    if type(history[0]) != tuple:
        history = history[2:]
        prefix = 'Merge'

    for i in history:
        if len(i) < 8: continue
        de_gene = [[genes[i[4][j]],i[7][j],sn(i[5][j])]  for j in range(len(i[7])) if i[5][j]>score_threshold]
        if len(de_gene) == 0:
            logging.warning("Historically, first column tried access undeclared variable j for genes[i[4][j]]. "
                            "This was likely a typo! Accessing genes[i[4][0]]")
            de_gene = [[genes[i[4][0]], i[7][0], sn(i[5][0])]] # merging
        if len(de_gene) < 10:
            num_plot = len(de_gene)+1
        else:
            num_plot = 10
        plt.figure(figsize=(8*num_plot,6))
        plt.subplot(1,num_plot,1)
        plot_labels_legend(x1,x2,str_labels_to_ints(i[1]))
        history_l = str_labels_to_ints(i[1])
        cor_all, p = pearsonr(np.mean(X[history_l==0, :], axis=0), np.mean(X[history_l==1, :], axis=0))
        if prefix == 'Merge':
            plt.title('%s %d: '%(prefix,ii)+'score = '+de_gene[0][2])
        else:
            plt.title('%s %d: '%(prefix,ii)+'%d genes express differently\ncorrelation %.2f'%(len(de_gene), cor_all))
        for j in range(0, num_plot-1):
            plt.subplot(1,num_plot,j+2)
            plt.subplots_adjust(left=0.1, right=0.9, top=0.85, bottom=0.05)
            plt.scatter(x1,x2,edgecolors='none',c=X[:,i[4][j]])
            g_mean_l = np.mean(X[history_l==i[7][j],i[4][j]])
            g_mean_s = np.mean(X[history_l==1,i[4][j]]) if i[7][j] == 0 else np.mean(X[history_l==0,i[4][j]])
            fold = g_mean_l/g_mean_s
            dr0 = np.mean(X[history_l==0,i[4][j]].astype(bool))
            dr1 = np.mean(X[history_l==1,i[4][j]].astype(bool))
            er0=np.where(np.argsort(np.mean(X[history_l==0,:], axis=0))[::-1]==i[4][j])[0]
            er1=np.where(np.argsort(np.mean(X[history_l==1,:], axis=0))[::-1]==i[4][j])[0]
            plt.title('%s: score = '%(de_gene[j][0])+de_gene[j][2]+
                      '\nc%d: %.2f (%.2f fold)\ndetection rate %.2f %.2f\nranked: %d %d'%(de_gene[j][1], g_mean_l, fold,
                                                                                  dr0, dr1, er0, er1))
            plt.colorbar()
            _ = plt.axis('off')
        if save_name is not None:
            plt.savefig(save_name+'_'+str(ii)+'.png', format='png', dpi=300)

        ii += 1


# Print the history in an interpretable way
def print_history(genes,history):
    prefix = 'Pre-split'
    if type(history[0]) != tuple:
        logging.info('%d of %d samples are singletons'%(np.sum(np.bincount(str_labels_to_ints(history[0]))==1),len(history[0])))
        print_singleton_merging_result(history[0],history[1])
        history = history[2:]
        prefix = 'Post-merge'

    for i in history:
        if len(i) < 8:
            continue

        if type(i[0]) == dict:
            nP = max(i[0].values())
            for key in i[0]:
                if i[0][key] < nP and key[-1] == 'L': nL = i[0][key]
                if i[0][key] < nP and key[-1] == 'R': nR = i[0][key]
            logging.info(prefix+": %-5s  L: %-5s  R: %-5s  Score: "%(nP,nL,nR)+ \
                  sn(i[6])+"  Top Gene: %-10.10s  Top Gene Score: "%(genes[i[4][0]])+sn(i[5][0]))


# Save full splitting results
def save_history(genes, history, score_threshold=10, num_genes=50, output_header=None, save_name=''):
    prefix = '_Split'
    ii = 1

    if type(history[0]) != tuple:
        history = history[2:]
        prefix = '_Merge'

    with open(save_name + prefix + '_history.csv', 'w') as f:
        rb = csv.writer(f)
        for row in output_header:
            rb.writerow(row)
        rb.writerow(['Split', 'Gene', 'p-Value', 'Larger_Cluster'])

        for i in history:
            if len(i) < 8: continue
            idx = 0
            for j in range(len(i[7])):
                if i[5][j] > score_threshold:
                    idx += 1
                    if idx > num_genes:
                        break
                    rb.writerow(list(map(str, [ii, genes[i[4][j]], sn(i[5][j]), i[7][j]])))
            ii += 1


# Print more information on how singletons were handled
def print_singleton_merging_result(L1, L2):
    # Make sure L1 is the result of merging singletons (has less unique clusters)
    if len(np.unique(L1)) > len(np.unique(L2)): L1,L2 = L2,L1
    outliers = [str(i) for i in np.where(L1 == -1)[0]]
    logging.info('Singleton(s) '+', '.join(outliers)+' marked as outliers (N = %d)'%(len(outliers)))
    for c in np.unique(L1):
        if c != -1:
            s = ''
            large_clust = None
            Ltemp = L2[L1 == c]
            if len(np.unique(Ltemp)) > 1:
                for j in np.unique(Ltemp): 
                    if np.sum(Ltemp == j) == 1: s += str(j)+', ' 
                    else: large_clust = j
                if large_clust is not None:
                    logging.info('Singleton(s) '+s[:-2]+' merged with cluster %d (N = %d) to form cluster %d (N = %d)' \
                          %(large_clust,np.sum(Ltemp == large_clust),c,np.sum(L1==c)))
                else:
                    logging.info('Singleton(s) '+s[:-2]+' merged to form cluster %d (N = %d)' \
                          %(c,np.sum(L1==c)))


# Analyze a split
def analyze_split(X, x1, x2, genes, history, split_num,num_genes=12, clust=None, show_background=True):
    prefix = 'Split'
    if type(history[0]) != tuple: 
        history = history[2:]
        prefix = 'Merge'

    ii = 0
    for k,i in enumerate(history):
        if len(i) == 8: ii += 1
        if ii == split_num: break

    i = history[k]
    y = str_labels_to_ints(i[1])
    important_gene_indices = i[4]
    scores = i[5]

    comparisons= i[7]

    if clust is not None:
        if num_genes > np.sum(comparisons == clust): 
            logging.info('num_genes is greater than the number of genes more highly expressed in clust %d'%(clust))
            num_genes = np.sum(comparisons == clust)
        important_gene_indices = important_gene_indices[comparisons == clust]
        scores = scores[comparisons == clust]
        comparisons = comparisons[comparisons == clust]

    plot_labels_legend(x1,x2,y)

    if show_background is False:
        x1 = x1[y != 2]
        x2 = x2[y != 2]
        X = X[y != 2,:]
        
    plt.title('%s %d: '%(prefix,split_num))
    for j in range(num_genes):
        if j%4 == 0: plt.figure(figsize=(16,4))
        plt.subplot(1,4,j%4+1)
        plt.scatter(x1,x2,edgecolors='none',c=X[:,important_gene_indices[j]])
        plt.title(genes[important_gene_indices[j]]+' (%d): '%(comparisons[j])+sn(scores[j]))
        _ = plt.axis('off')


# Select features from history
def feature_extraction_via_split_depth(history):
    feature_lists = []
    for i in history:
        feature_lists.append(i[1][:np.max([20-i[-1],1])])
    return reduce(np.union1d,feature_lists)


# Remove entries from shistory that would not be there if the splitting step was performed 
# with the specified threshold
def filter_out_extraneous_steps(shistory, threshold):
    shistory_filt = [i for i in shistory if len(i) < 8 or i[6] > threshold]
    a = [sorted(list(i[0].keys()), key=lambda x:len(x))[0] for i in shistory_filt]
    flags = [False for i in range(len(a))]
    keep_entries = {i:None for i in a}
    # Given a list of strings, for each string s of length N, remove it if s[:-1] is not in the list
    prev_len = len(keep_entries)
    while True:
        for j, i in enumerate(list(keep_entries.keys())):
            if len(i) > 1 and i[:-1] not in keep_entries:
                keep_entries.pop(i)
        if prev_len == len(keep_entries): break
        prev_len = len(keep_entries)
    
    for j, i in enumerate(a):
        if i not in keep_entries: flags[j] = True
        
    return [i for j,i in enumerate(shistory_filt) if not flags[j]]


# Perform parameter sweeping using the above function
def get_clusters_from_history(D, shistory, threshold, disband_percentile):
    shistory = filter_out_extraneous_steps(shistory,threshold)
    
    # For a given threshold, recover the clustering from D
    N = len(D)
    labels = ['r' for i in range(N)]
    for q,i in enumerate(shistory): _,labels = update_labels(i[1],labels)
    labels = np.array(labels).astype('S%d'%(max([len(i) for i in labels])+len(str(N))))

    # Remove singletons
    disband_threshold = np.percentile(flatten_distance_matrix(D), disband_percentile)
    for i in np.unique(labels):
        inds = labels == i
        if np.sum(inds) > 1 and np.sum(flatten_distance_matrix(D, inds) < disband_threshold) == 0:
            for j in range(N):
                if inds[j]:
                    labels[j] = b"".join((labels[j], str(j).encode('utf-8')))
                    
    return labels
