from collections import namedtuple, deque, defaultdict
from itertools import combinations
import cython


ReducedCluster = namedtuple('ReducedCluster', ['counts', 'child_node'])


def seq_error_correction(umi_counter, cell_index, gene):
    """
    Args:
        umi_counter: dictionary that relates the UMI to the count of its reads
        cell_gene_entry: tuple representing the cell-gene combination to which the UMI belong

    Returns: A molecule annotation in the format:
    [cell-label, barcode, gene, read_count, rsec_corrected_read_count, parent_barcode]
    If the barcode is a parent -- that is to say, it is not an error -- then parent_barcode is empty.
    """

    adj_list = get_adjacency_list(umi_counter, 1)
    clusters = get_connected_components(adj_list, umi_counter)
    reduced_clusters, child_to_parent_umi = reduce_clusters(clusters, umi_counter)

    for umi, umi_count in umi_counter.items():
        try:  # Child node
            parent_umi = child_to_parent_umi[umi]
            rsec_corrected_count = 0
        except KeyError:  # Parent node
            parent_umi = ''
            rsec_corrected_count = reduced_clusters[umi].counts
        finally:
            yield [cell_index, umi, gene, umi_count, rsec_corrected_count, parent_umi]


def reduce_clusters(clusters, umi_counts):
    """
    Args:
        clusters: the set of reads that
        umi_counts: dictionary that relates the UMI to the count of its reads

    Returns: A smaller set of clusters labeled with a representative UMI
    """
    reduced_clusters = {}
    child_to_parent_umi = {}
    for cluster in clusters:
        if cluster:
            sorted_nodes = sorted(cluster, key=lambda x: umi_counts[x], reverse=True)
            parent_umi = sorted_nodes[0]
            child_umis = sorted_nodes[1:]
            collapsed_reads = sum(umi_counts[umi] for umi in cluster)
            reduced_clusters[parent_umi] = ReducedCluster(counts=collapsed_reads, child_node=child_umis)
            for child in child_umis:
                child_to_parent_umi[child] = parent_umi
    return reduced_clusters, child_to_parent_umi


@cython.cfunc
@cython.locals(
    adj_list=dict,
    umis=list,
    count_a=cython.int, count_b=cython.int, umi_length=cython.int,
)
@cython.returns(dict)
def get_adjacency_list(umi_counts, threshold=1):
    """
    Args:
        umi_counts: dictionary that relates the UMI to the count of its reads
        threshold: hamming distance above which two reads are considered different UMIs

    Returns: An adjacency list where the parent UMI has at least twice the read count of its children UMIs
    """

    umis = list(umi_counts.keys())
    adj_list = {umi: [] for umi in umis}

    if len(umis) > 25:  # If there are more than 25 distinct UMIs, use hashed comparisons
        umi_length = len(umis[0])
        substr_idx = build_substr_idx(umis, umi_length, threshold)
        iter_umi_pairs = iter_nearest_neighbours(umi_counts, substr_idx)
    else:  # With fewer (or equal to) 25 UMIs, it is faster to consider all pair-wise comparisons without hashing
        iter_umi_pairs = combinations(umi_counts.items(), r=2)

    for (umi_a, count_a), (umi_b, count_b) in iter_umi_pairs:
        if count_a >= (count_b * 2):
            if edit_distance(umi_a, umi_b, threshold):
                adj_list[umi_a].append(umi_b)
        elif count_b >= (count_a * 2):
            if edit_distance(umi_b, umi_a, threshold):
                adj_list[umi_b].append(umi_a)
    return adj_list


@cython.cfunc
@cython.locals(
    first=str, second=str,
    threshold=cython.int, num_diff=cython.int, seq_len=cython.int, i=cython.int
)
@cython.returns(cython.bint)
def edit_distance(first, second, threshold):
    """
    Args:
        first, second: UMIs being compared
        threshold: hamming distance above which two reads are considered different UMIs

    Returns: True if the edit distance/hamming distances between its two arguments is greater than the given threshold
    """
    num_diff = 0
    seq_len = len(first)
    for i in range(seq_len):
        if first[i] != second[i]:
            num_diff += 1
        if num_diff > threshold:
            return False
    return True


@cython.cfunc
@cython.locals(
    adj_list=dict,
    found=set, component=set,
    node=str
)
@cython.returns(list)
def get_connected_components(adj_list, umi_counts):
    """
    Args:
        adj_list: graph that connects reads to their erroneous copies
        umi_counts: dictionary that relates the UMI to the count of its reads

    Returns: unordered clusters of connected UMIs
    """
    found = set()
    components = list()

    for node in sorted(adj_list, key=umi_counts.get, reverse=True):
        if node not in found:
            component = breadth_first_search(node, adj_list)
            for umi_key in component:
                adj_list.pop(umi_key, None)
            found.update(component)
            components.append(component)
    return components


def breadth_first_search(node, adj_list):
    """
    Args:
        node: node with which to begin the search
        adj_list: graph that connects reads to their erroneous copies

    Returns: unordered set of all children of the root node

    """
    found, queue = set(), deque([node])
    while queue:
        node = queue.pop()
        if node not in found:
            if node in adj_list:  # UMI node has more than one parent
                queue.extend(adj_list[node])
                found.add(node)
    return found


# UMI tools' hashed network construction algorithm
# lightly modified from https://github.com/CGATOxford/UMI-tools


def get_substr_slices(umi_length, idx_size):
    """
    Create slices to split a UMI into approximately equal size substrings
    Returns a list of tuples that can be passed to slice function
    """
    cs, r = divmod(umi_length, idx_size)
    sub_sizes = [cs + 1] * r + [cs] * (idx_size - r)
    offset = 0
    slices = []
    for s in sub_sizes:
        slices.append((offset, offset + s))
        offset += s
    return slices


def build_substr_idx(umis, umi_length, min_edit):
    """
    Build a dictionary of nearest neighbours using substrings, can be used to reduce the number of pairwise comparisons.
    """
    substr_idx = defaultdict(lambda: defaultdict(set))
    slices = get_substr_slices(umi_length, min_edit + 1)
    for idx in slices:
        for u in umis:
            u_sub = u[slice(*idx)]
            substr_idx[idx][u_sub].add(u)
    return substr_idx


def iter_nearest_neighbours(umi_counts, substr_idx):
    """use substring dict to get (approximately) all the nearest neighbours to each in a set of umis."""
    for umi, umi_count in umi_counts.items():
        neighbours = set()
        for idx, substr_map in substr_idx.items():
            u_sub = umi[slice(*idx)]
            neighbours = neighbours.union(substr_map[u_sub])
        neighbours.remove(umi)
        for nbr in neighbours:
            yield (umi, umi_count), (nbr, umi_counts[nbr])