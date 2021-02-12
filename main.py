import random
import copy
import sys
import itertools


############################################################################################################
# WEEK 1
############################################################################################################


def string_composition(k, text):
    """ Takes an integer k representing the length of k-mer, and text (String) which is the given DNA
    sequence. It then outputs the composition of the sequence, i.e. the collection of k-mers that make
    up the sequence. """

    coll = []
    for i in range(len(text) - k + 1):
        coll.append(text[i: i + k])

    return sorted(coll)


def string_spelled_by_genome_path(genome_path):
    """ Takes in a collection of strings, genome path (List), and reconstructs the sequence. The genome
    path is the collection of k-mers that when placed in sequence and overlap with one another, spell
    the correct genome """

    reconstructed_genome = genome_path[0]
    for i in range(1, len(genome_path)):
        if genome_path[i][0:-1] == genome_path[i - 1][1:]:
            reconstructed_genome = reconstructed_genome + genome_path[i][-1]

    return reconstructed_genome


def generate_unique_key(base_pattern, collection):
    """ Takes in a base pattern that is to be made unique and a collection to check against """

    if base_pattern in collection.keys():
        pattern = str(base_pattern.split('_')[0])
        key_num = int(base_pattern.split('_')[1])
        key_num = key_num + 1
        new_key = pattern + '_' + str(key_num)
    else:
        new_key = base_pattern

    return new_key


def overlap_graph(kmers):
    """ Takes in a collection of strings, kmers (List), and generates the overlap graph in the form of an
    adjacency list. """

    overlap_dict = {}
    for kmer_outer in kmers:
        if kmer_outer not in overlap_dict:
            unique_key = kmer_outer
            overlap_dict[unique_key] = []
        else:
            test_key = kmer_outer + '_2'
            unique_key = generate_unique_key(test_key, overlap_dict)
            overlap_dict[unique_key] = []
        for kmer_inner in kmers:
            if kmer_outer[1:] == kmer_inner[0:-1] and \
                    kmer_inner not in overlap_dict[unique_key]:
                overlap_dict[unique_key].append(kmer_inner)

    return overlap_dict


def debruijn_graph(k, text):
    """ Takes in an integer k and a DNA sequence, txt (String), and outputs the De Bruijn Graph. The De Bruijn Graph
    for a sequence is one that has edges that contain each k-mer and nodes that contain each (k-1)-mer. The graph is
    then collapsed such that repeats of nodes are joined together and their adjacent nodes are combined. """

    node_length = k - 1
    nodes = [text[i: i + node_length] for i in range(len(text) - node_length + 1)]

    al = {}
    for i in range(len(nodes) - 1):
        if nodes[i] not in al:
            al[nodes[i]] = [nodes[i + 1]]
        else:
            al[nodes[i]].append(nodes[i + 1])
    return al


def debruijn_graph_from_kmers(kmers):
    """ Takes in a collection of kmers (List), and outputs the De Bruijn graph as a Dictionary. It works by generating
    isolated nodes and edges where each edge is a directed edge that is labelled by a given kmer. The nodes that are
    on either side of the edge are made from the kmer's prefix and suffix (in that order). Identical nodes are then
    glued together, generating the final De Bruijn graph. """

    dbg = {}
    for kmer in kmers:
        if kmer[0:-1] not in dbg:
            dbg[kmer[0:-1]] = [kmer[1:]]
        else:
            dbg[kmer[0:-1]].append(kmer[1:])

    return dbg


#######################################################################################################
# WEEK 2
#######################################################################################################
def unexplored_edges(graph):
    """ Helper function for eulerian_cycle() function. Takes in a graph (Dict) and checks to see if there are still
    unexplored edges """

    for k, v in graph.items():
        if len(v) > 0:
            return True

    return False


def recycle(new_start, cycle):
    """ Helper function for eulerian_cycle() function. Takes in the newly found starting position (String) and the
    most current cycle (List). It then rearranges the previous cycle to begin at the new start and returns the new
    cycle. """

    del cycle[0]
    new_start_index = cycle.index(new_start)
    for i in range(new_start_index):
        node = cycle.pop(0)
        cycle.append(node)
    cycle.append(new_start)
    return cycle


def get_new_start(cycle, graph):
    """ Helper function for eulerian_cycle function(). This takes in a cycle (List) and a graph (Dict) and outputs
    the new node from which to begin traversal in the next iteration of the algorithm. """

    for n in cycle:
        if len(graph[n]) > 0:
            return n


def walk(start, tracker, cycle):
    """ Helper function for eulerian_cycle() function. Takes in a starting point (String), a tracker graph (Dict) that
    tracks the edges that have been visited, and a cycle (List) which maintains the most current Eulerian cycle. This
    is a recursive function that traverses the given graph and builds the current Eulerian cycle.
    *Note* for very large datasets, the Max Recursion Limit (1000) is exceeded. In practical implementations,
    an iterative approach for walking may be better. """

    if len(tracker[start]) > 0:
        step = tracker[start].pop(0)
        cycle.append(step)
        cycle = walk(step, tracker, cycle)
    else:
        return cycle
    return cycle


def eulerian_cycle(graph):
    """ Takes in a graph (Dict) that contains an adjacency list for an Eulerian directed graph. It outputs the Eulerian
    cycle in the graph. This works by first finding a cycle by traversing the given graph while only visiting each
     edge once. If a cycle is found and there are still unvisited edges in the given graph, then it finds a node
     from the newly generated cycle that still contains unvisited edges. It traverses the graph again, this time starting
     at the node with unused edges, walks through the edges of the previously generated cycle and continues onward
     from there randomly. Each time the traversal gets 'stuck', it repeats this procedure until the entire graph is
     traversed. """

    tracker_graph = copy.deepcopy(graph)
    k = random.choice(list(graph.keys()))
    cycle = [k]
    cycle = walk(k, tracker_graph, cycle)
    while unexplored_edges(tracker_graph):
        new_start = get_new_start(cycle, tracker_graph)
        cycle = recycle(new_start, cycle)
        cycle = walk(new_start, tracker_graph, cycle)

    return cycle


def find_new_endpoint(missing_in, missing_out, cycle):
    """  Helper function for Eulerian Path function. Takes in the missing edge (Strings) as well as the cycle (List)
    and outputs the index of where the cycle should end. The new endpoint for the cycle should be the node that did
    not contain an outgoing edge in the original graph. """

    for i in range(len(cycle) - 2 + 1):
        if cycle[i] == missing_out and cycle[i + 1] == missing_in:
            return i


def transform(missing_in, missing_out, cycle):
    """  Helper function for Eulerian Path function. Takes in the missing edge (Strings) as well as the cycle (List)
    and outputs the transformed cycle in which the missing edge is no longer there. It first finds the new endpoint
    and re-arranges the sequence of the cycle such that there is no edge where there shouldn't be. """

    # endpoint = cycle.index(missing_out)
    endpoint = find_new_endpoint(missing_in, missing_out, cycle)
    del cycle[0]
    for i in range(endpoint):
        node = cycle.pop(0)
        cycle.append(node)
    return cycle


def get_in_counts(graph):
    """  Helper function for Eulerian Path function. Takes in a graph (Dict) and outputs another Dict containing
    the count of incoming edges for all nodes. """

    in_counts = {}
    for k, v in graph.items():
        for n in v:
            if n not in in_counts:
                in_counts[n] = 1
            else:
                in_counts[n] = in_counts[n] + 1

    return in_counts


def get_out_counts(graph):
    """  Helper function for Eulerian Path function. Takes in a graph (Dict) and outputs another Dict containing the
    count of outgoing edges for all nodes. """

    out_counts = {}
    for k, v in graph.items():
        if k not in out_counts:
            out_counts[k] = len(v)
        else:
            out_counts[k] = out_counts[k] + len(v)

    return out_counts


def get_missing_edge(in_counts, out_counts):
    """ Helper function for Eulerian Path function. Takes in in_counts and out_counts (Dict's) and outputs a tuple
     representing the missing edge that would make the graph balanced. This missing edge is added to the graph,
     the Eulerian cycle for this balanced graph is generated and the Eulerian Path is derived from that through a
     simple re-arrangement of the cycle (transform()). """

    for k, v in in_counts.items():
        if k not in out_counts:
            missing_out = k
            break
        elif out_counts[k] < in_counts[k]:
            missing_out = k
            break
        else:
            missing_out = ''

    for k, v in out_counts.items():
        if k not in in_counts:
            missing_in = k
            break
        elif out_counts[k] > in_counts[k]:
            missing_in = k
            break
        else:
            missing_in = ''

    return missing_in, missing_out


def get_balanced_graph(graph, missing_in, missing_out):
    """ Helper function for Eulerian Path function. Takes in a graph (Dict) and the missing edge (Strings) and outputs
    the balanced graph. The balanced graph simply contains the missing edge that would allow the possibility of an
    Eulerian cycle. """

    if missing_out not in graph:
        graph[missing_out] = [missing_in]
    else:
        graph[missing_out].append(missing_in)

    return graph


def eulerian_path(graph):
    """ Takes in a graph (Dict) and outputs the Eulerian Path. The graph is assumed to be nearly balanced, meaning that
     there is only one missing edge that needs to be added in order for an Eulerian cycle to exist in the graph. Once
     the missing edge is added, the Eulerian cycle is found. The cycle is then re-arranged such that the missing edge
     is removed, i.e. the node that originally did not contain an outgoing edge is at the end of the cycle and the node
     that originally did not contain an incoming edge is at the beginning of the cycle. """

    in_counts = get_in_counts(graph)
    out_counts = get_out_counts(graph)

    missing_in, missing_out = get_missing_edge(in_counts, out_counts)

    if missing_in != '':
        balanced_graph = get_balanced_graph(graph, missing_in, missing_out)
        cycle = eulerian_cycle(balanced_graph)
        path = transform(missing_in, missing_out, cycle)
    else:
        path = eulerian_cycle(graph)

    return path


def string_reconstruction(k, patterns):
    """ Takes in an integer k and a list of k-mers and assembles the genome. It first builds the de Bruijn graph which
    is a Dict containing each k-mer comprising the edges and each (k-1)-mer comprising the nodes. From this graph
    an Eulerian path is found. The Eulerian path for a given set of k-mers represents a sequence of nucleotide
    snippets that when combined, the ending nucleotides (suffix) of a snippet match the beginning nucleotides (prefix)
    of the subsequent snippet. The final genome is then built by connecting the snippets in the proper sequence and
    removing overlapping portions of the snippets. """

    db = debruijn_graph_from_kmers(patterns)
    path = eulerian_path(db)
    s = string_spelled_by_genome_path(path)

    return s


def generate_binary_kmers(k):
    """ Helper function for k_universal_circular_string() function. This function simply creates the collection of
     k-mer patterns comprised of all values from 0 to 2^k in binary. This collection is then used to build the
     k-universal circular string which is a string that contains all binary values from 0 to 2^k in one string. This
     string is also 'circular' in that the ending nucleotides match the beginning nucleotides."""

    patterns = []
    for s in range(2 ** k):
        patterns.append('{0:0{1}b}'.format(s, k))

    return patterns


def k_universal_circular_string(k):
    """ Takes in an integer k, and generates the k-universal circular string which is a string that contains all binary
     values from 0 to 2^k in one string. This string is also 'circular' in that the ending nucleotides match the
     beginning nucleotides. """

    patterns = generate_binary_kmers(k)
    db = debruijn_graph_from_kmers(patterns)
    path = eulerian_cycle(db)
    s = string_spelled_by_genome_path(path)
    end_index = (k - 1) * -1

    return s[0:end_index]


def string_spelled_by_gapped_patterns(gapped_patterns, k, d):
    """ Takes in a List of gapped_patterns and integers k and d and outputs the string that contains perfect overlap
    between the prefixes and suffixes of the gapped patterns. """

    first_patterns = []
    second_patterns = []
    for kdmer in gapped_patterns:
        beg_end = kdmer.split('|')
        first_patterns.append(beg_end[0])
        second_patterns.append(beg_end[1])

    prefix_string = string_spelled_by_genome_path(first_patterns)
    suffix_string = string_spelled_by_genome_path(second_patterns)

    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return "There is no string spelled by the gapped patterns"

    s = prefix_string + suffix_string[-(k + d):]

    return s


def debruijn_graph_from_read_pairs(kdmers):
    """ Similar to the dbruijn_graph_from_kmers() function except that it takes in a collection of (k,d)-mers and
     generates the debuijn graph. A (k,d)-mer read pair is a pair of k-mers in the genome that are a distance d
      apart from each other. """

    dbg = {}
    for kdmer in kdmers:
        beg_end = kdmer.split('|')
        prefix = beg_end[0][0:-1] + '|' + beg_end[1][0:-1]
        suffix = beg_end[0][1:] + '|' + beg_end[1][1:]
        if prefix not in dbg:
            dbg[prefix] = [suffix]
        else:
            dbg[prefix].append(suffix)

    return dbg


def string_reconstruction_from_read_pairs(k, d, gapped_patterns):
    """ Takes in integers k-mer length k and gap distance d, as well as a List of read-pairs, gapped_patterns. It
    outputs a valid genome sequence built from the given (k,d)-mer composition. """

    dbg = debruijn_graph_from_read_pairs(gapped_patterns)
    path = eulerian_path(dbg)
    s = string_spelled_by_gapped_patterns(path, k, d)

    return s


def in_found_cycle(paths, node, edge):
    """ Helper function for find_isolated_cycles() function. It takes in a List of paths, a node, and it's adjecent node,
     edge, and checks to see whether the node-edge pair has already been found. This prevents isolated cycles to be
     captured redundantly in the final contigs collection. """

    for path in paths:
        if node in path and edge in path:
            return True

    return False


def find_isolated_cycles(graph, in_counts, out_counts):
    """ Helper function for maximal_nonbranching_paths() function. Takes in a graph (Dict), in_counts, and out_counts
    (Dict's) and outputs a List of isolated cycles. An isolated cycle is a path that contains only 1-in-1-out nodes
    and begins and ends at the same node. Isolated cycles are not captured by the maximal_non_branching_paths
    algorithm because it only captures cycles that begin with branching nodes. """

    paths = []
    for k, v in graph.items():
        if (k in in_counts and k in out_counts) and (in_counts[k] == 1 and out_counts[k] == 1):
            cycles = [k, v[0]]
            w = v[0]
            if in_found_cycle(paths, k, w):
                break
            while (w in in_counts and w in out_counts) and (in_counts[w] == 1 and out_counts[w] == 1):
                cycles.append(graph[w][0])
                if graph[w][0] == cycles[0]:
                    break
                w = graph[w][0]
            if cycles[0] == cycles[len(cycles) - 1]:
                paths.append(cycles)

    return paths


def maximal_nonbranching_paths(graph):
    """ Helper function for the contig_generation() function. Takes in a de Bruijn graph (Dict) and outputs all maximal
    nonbranching paths (List). A maximal nonbranching path is one where all intermediate nodes of the path are
    1-in-1-out nodes. These represent portions of the final genome that are likely to be contiguous. """

    in_counts = get_in_counts(graph)
    out_counts = get_out_counts(graph)
    paths = []
    for node in graph:
        if node not in out_counts:
            continue
        if (node not in in_counts) or (out_counts[node] != 1 or in_counts[node] != 1):
            if out_counts[node] > 0:
                for edge in graph[node]:
                    nonbranching_path = [node, edge]
                    w = edge
                    # if w in in_counts and w in out_counts:
                    while (w in in_counts and w in out_counts) and (in_counts[w] == 1 and out_counts[w] == 1):
                        nonbranching_path.append(graph[w][0])
                        w = graph[w][0]

                    paths.append(nonbranching_path)

    isolated_cycles = find_isolated_cycles(graph, in_counts, out_counts)
    for cycle in isolated_cycles:
        paths.append(cycle)

    return paths


def contig_generation(kmers):
    """ Takes in a collection of kmers (List) and outputs all contigs given by all maximal non-branching paths within
     the De Bruijn graph. """

    strings = []
    dbg = debruijn_graph_from_kmers(kmers)
    paths = maximal_nonbranching_paths(dbg)
    for path in paths:
        s = string_spelled_by_genome_path(path)
        strings.append(s)

    return strings


if __name__ == '__main__':
    ######################################################################################################
    # WEEK 1
    ######################################################################################################

    #######################################################################################################
    # String Composition Challenge - Generate the k-mer composition of a given DNA sequence. The k-mer composition
    # is the collection of all k-mer substrings in the DNA sequence ordered lexicographically
    #######################################################################################################

    # k_arg = 100
    # t_arg = 'TGCTTACTTCGAGGCTTGTTGGGGCGTCTGTACTATTGTATGCGTGACCGGCCGCACGTGGAACAGGCGGCCAATTGTCCATTCAGTCTAGCTTGGTCCCGGTATGATCCGCGTCGAAATGTTGGTGCACCACAATCCCGGGAAGCCTACAGAAATATAACCATCCCGCGGCGCGGCCCGATTTGGTCCACAGGGACACTCACACAATGATATTTAGTCGGGCCACGGGTGCGTGTAATGATTCTCCGCGGAACTGCGAATCCATAGGTACAGTCCTCAGATGCACTGATCGAAAAATAGGATTGCTGAGCGAAGAAGGAACTACCCGAGATGTCTGCTGTTGAGACGGCTTTCGTACCTGTCTGTGACGGCATGAGCCCGTAGCAGTCGCCTAATTCGGAGGCTATGAGGATGGAATGATCGCTGAACGGGTCTCGCGATCGACCGATTACCACCGACCCGTGGGCCTTATATTAGAATAATTACACTATACGGCCGTGAGCGTGTTAATTAGGCTGCAAATACCCAATGGTCTAGACCATATGTCCTCAGTATCAGGATGCTGGCCGGGTTGGGTCTGATCCTAACAGCTCCGCGCCACCGTTCTAGCCCTTACGCAGCTGTGGTAGGGGACATGGCAGGAATGCATTGTAACGCCTATTCACTATCGTTCTATGTCGTCTGCTTAGTTCTAGCTGTACTCGACATTCCAGAGAGCAAAGAAGCAACTACCGCGGATCAACTCCGGGACGGGTACGAACTTCCGGGTGAAGGACTGGTGCTAAACACTATTCTACTCTCCACATAGGCGAGGCTTCTATTATCGATTTTCAGCATTAATGCCAATAAATGGTGTTCGGCGACGTATCACGAAGGTCTGTTCCTCAGGCGAACGTGAACACTCGCACATGCATGTTGCTGCTAAAAGTCAGATTCAGGTATAGTGGTACTGGTGGAGGCCCCGCGTCTGCGGGATGTCTCGGAGTTTCTCATATTCCATCCGAGGAGTTCTGTTCCCTGTACGTTGTAAAGATAGGGAAGAGGTATTGGGATATACCTTAACGACAAAAACTGTTTTATTAAAGGACGCACCGTATCTTAAAGCTGATTTCTCATTGATGACCAGTCATCCAGAGGTACCCGACTATCTGGCCGTGGGTCTTGACATCTATACCCAGTCCAGAGCCACAGAGGTTCCTAGGATGGCTCTCAGGTTGTATGATAGGTTCTGGGTAACTAAAACTGCTCTGTTCATTCCACCTTGTACACCTAACACCCGGTAGCTACCCACGCGAGTCTGCGTATACGGGTTACACGAAAATAAACGTCCCTACGCTTCGTTGTGCGCAATCCTCTCTACAAAGGACTGGCAAAATACCGCCCGTATACCGGGTTAGAGGCAGGCCAAAGCACAGTTACGCGTTTGTTTGCCTTTGGAATAAGAGATGCATTAACCTATCCTGCAGAGGGGCGTGAGTCTGTACACACTGCAAGTTCGCCCTAGAACGGACTCTCAGTTAGGCATATGCAGGGTTTGCTGAACACTCATAGCAGCCTGATTTGCCCCCGCGTCGATGATTTCTTACCCAAGGGAAACACAGATCCGATGAATTGGCGCATAGAATGTTTTTTCGCGCTAGCACTTTAGAATAGATTACCGCGAGTGGCCAATCTGAACGCTAGAAGCACACCTGTACTCCGCGACGGGAGCCGGGGGCAGCGTCCAGCGGACCATCCTCGACGTCTCAACTAGGCTCAAGGCGCGCTCCCGGCGCGTGAAATTCCATAGCAGTGGTTACAGATTGAGTGTTTCGTGCCTCCTAAGCCGTTACGAGGCCATTATAGGGCTTACGCCGTCATTACGAGTGTCTTCACCGGGCTACGGCCGCTAGCGTCTCCCACTTGACCAGCATGATCTCATAAAAGAACGACGTCTCTGTTAGTACTCTGATGTCACTAATACAGTGAAAGCGTGACATACCGGACAAACCTCATGGCGAATGTATGCTAGTTTCTAGCCGCGCCCGTGGGAACCAAATTAGCATGCTAGGAAGTGACTATTACCGCACTGTGGTTGAACTGCAAAGCTCAGAAATGGTGCCCTGGATCTTCAGATTTTTTGTGACAGATCTAGTGGGCTAACCCCCTTGGTTTATGTAATTACGATGTGTCGTACGGGCGACGGGAAATTCATCTCTTTATGGCTGACAATGCGGAAAGGGAGCACTACCAACCACACTTTCCCGGACGGGGCAACGACAAGTCATAACTTAACCAGATTATTAGTAAAAGGCCGGTGTGTGCCGATCTCTACGGGCAACAGCGGCCACTACCGCGTGAAGTATCGATTCCTCCGCCTAGTCAACATTCGGGCGCTTGCGTATAGTACGAGCCCTGTCCTATTGACAAAGCCCACGCCCTGTGGCCTAATTTTCCTCCTGATAATGGAAAGGTCGGGTGAAAAATAGTCACCACAGTTAGATTACCTCGCCTTTTCAGCCAGACATATATCATTAGATGAACGTCACGTAGCGTGCCTCTGGGCCGTGACCGTCATCAAGCACCGTGCAAGAATCCTCCCTCCGTCCTTTTAAAGGGTGTCACGTAGGATAAACATGGATTAGGATGAATGATTGCAAATGAAGAGCCTTTAGGGCTGTCAGCGACCCACCGACCCGGCAGGAACTCCGGCACCATTCATATCATGGACGTGCCTCAGCGTCTTACTTTGGGTGAGGGCTAATGAAATCGGCTCCGCCTGCAGGTCGCGCTATACCATAGTGGACCAGATGGCCGCAAGAGACTGTACGTTGCATGAGGGCACCGTGCAGAATCGTACAGCGATTAAGAGATTTGCACCTGCCGAGTCAGACGCATAGCTACGCGGCCTCTGGCGGTAGCTGGTCCCCGCTTTTTTCCGCCACCGAAGGCGTAGAGTGCCTATAAGATCCGCTCACGCAACGAGCCACCTGGGATTGGTTTGAGTAGGCTAGGTGCCTGTCTGCGTGGGAGCTTTAGTTGACTGCAGTTCTCTGTTAGTGCCACATCCGCGCCACGTTGCACCCTTCGTTATCACGCGTTACTCTAACGCAAGTCCCATGTAATACTCTAGCACTCTACTAGAATCTTAGACATTGTCAGACCCCCTACCCGTAATCCCTACATCGTTTGTTACGGTGTATGGTTTCAAGGCCTTTAGAGTAGGATCGACATGCGAATTCCTAGTTATTATTGAACTTTGCCTAATCGGATTATCTCGCATCTTGCGGGTTTACAGAGATCTTCTCGTTGCGGAACTCGGCCCCGGTCCGCACCACAGAACAAAATGCTGCTCCAAAAGAATATCTCCTGAGGGTGAATCTACGTGCGCGTGGTGGTGAGAGTACTCCGAATAGTTACATGGAAAGCCTCCCCAGGGAAAAAATATTTCCGCGATTAGCTCTCGAAACCAATTAACATACTGCGATCCACGAGCAGTGCTCCGCTGTAACCATCTTACCGCGGGGCTGGCGCTAGTCCGCATACGTCGTATGGGTGGCTAACAATCGATGGAGTCCACGACGCCCGCACTTCCGGCTTAACGGCTTCCAGTCTGCCCTTAAATGCGGCGTCCAGACCGCGATATCTAGGAGTAGGGTAGGCCGGACTCATAGGAATACCCGTCTAAACTATTTCTGGATGTTAAGTGGACTTCCTCGTTTGTGATCAAGTTGTATATTGCTTGTTTACTAGGGCAACGCATTATGGTAAACAGTTCAAGTAAAACAGCCGTATTGACTCAGTCCTCGGTCTAGGTCAACCAGAAATTGAAATTCACTTTATAAAAGCCAGGTGTCGCGAGGGTAACGAGGTCTGCTCATGCCACGTGAACCTGATGTACATGGTGGGCCCCTAATACATGGTCGTATCATTGGTTTTAAATCAGTCGTAGATAACGTCTACAGGATTGGTGTCACATTGCACCGGACCAGCTAGCACAAGGATACCTCATGACTGAGCTTATCGACAATTTGACTTTTGCCACGGTTACACCATGTTGAGATCATTTTATCGAGAAGTATCCACCCATCTGTCCTCCGCCGGAAGCCAGGCTGCTCTGTCCATCACATCATGAAATACCATGATAGTAGCGGACGAAACCGTACCGGCATATGTCGTGCAAAACGCAGCTGGTAATTACATGCGCTGTCGTAGCAGTTTTCCGTCCCTGACATTCTTCTCTTCAAAAAACCGTCCTTACAATGGTCGACGCAGTGGATGCGCTGATATAATGAAGGGTGATCTGATTGAAGATGTTGGATTAATATCCCTGTGACTAGTCAGTGTTGGGCCAAATCCCATTGCGAGTTCGTGCATAGCTCGACGTATCTCAAATTCGCCAGTAAGAGCAGCCGGACGACTTCAACAAGTTTACCCGAGGTACCGGGTCGGTTTATCCGTATGCCCTACCAGGCTCTCCATGTACAAAGCGCACAAAGTAGGATCGTCGAGAAACCGAAGGCAATTCTCATTATTCAGTTTGACACCATTTATATAGGCAGGAACGGTGCTGTAAGGCTGCTCAGATGGAGATTCCCTGTCGTCTCTGCGGCTAGAGCCCATCATCTATCTACCGTATTTGTCGTTGACTTCCCGGGACGGCTGCCTATCCAACAACTTATTACGTATAGGAATTGGGGCGTCGTGAGTAACCCGACCGTACTCATATAACTTCCATACTTACAGTCTTATTAACCGTCTGCCCATGAGTTGTCACACCCTCGCGGAAAGGGGACTGAATCCCGTCAGGATCCTAAACCGCAGGGTGGGTGGACTCTGAAAATTACAATCGTGGGGATGGCGCGCCCTCGGAAGAAAGAACACCACGAGAGAGTCTCCAACCTCCCACGTAATAATCTAACTTATAGCAAAGTTGTTAACCCCTATCCATGTGATTGCTAGCCCGACAGAACCGTTGGTCAACGGATCAATTATGCGCAGGTTGCGTGCGGATTCCGCGC'

    # c_out = string_composition(k_arg, t_arg)
    # [print(s) for s in c_out]

    #######################################################################################################
    # String Spelled by a Genome Path Challenge - Reconstruct a genome given the k-mers that comprise it.
    # The k-mers are in sequence and overlap, i.e. the genome path is given.
    #######################################################################################################

    #     path_arg = """
    # """

    # path_arg = path_arg.strip().split('\n')
    # r_out = string_spelled_by_genome_path(path_arg)
    # print(r_out)

    #######################################################################################################
    # Overlap Graph Challenge - Given a collection of k-mers, generate the overlap graph in the form of an
    # adjacency list. An adjacency list shows each node together with every other node that it connects to
    #######################################################################################################

    #     k_arg = """
    # """
    #     k_arg = k_arg.strip().split('\n')
    #     o_out = overlap_graph(k_arg)
    #     for k, v in o_out.items():
    #         if not len(v) == 0:
    #             print(k + ' -> ' + ','.join(v))

    #######################################################################################################
    # Generate a 4-universal binary string - A 4-universal binary string is a string that contains each of
    # the 16 possible binary strings of length 4 from 0000 to 1111
    #######################################################################################################

    # k_universal_string(16)

    #######################################################################################################
    # De Bruijn Graph from String Challenge - Find the De Bruijn Graph in the form of an adjacency list given an integer
    # k and DNA Sequence. The genome path is given. Repeats of nodes are 'glued' together such that the overall number
    # of nodes is reduced but the original DNA sequence can still be obtained.
    #######################################################################################################

    # k_arg = 12
    # t_arg = 'CTCTCGAGATCTACTCGTCACAACATCATCCCGATTACTTATAGTATGGCAGGGCGTGATTGTGCACCGGTCCTTCACCACGCTCCCTACGCTACACACCGCCGACTAGGGCTGTGGGCTCCCCGGGGGACTGTATAAGTTAAAGCCAAGGTTAGAACATCCCGAGAGGGCTTCCCGACGCTGCCTGTGTGCCGCGTATCCAGTGAACAAACAGAGGTTGCCCCGCCCTTTGTGGCATGAACAAAGAATGCTAGCCTCTCGCTCCTAATCCCTTACCGTCGGTACCGGGCTGAGTCCGTACTGTTGGCAGTTTGTAACATCCACTATCCCCACAGTTGAAAAACCACAAGGTTTATCGTCAAATGTGGCAACCCATGCTAGCTACGCGCCCGTCTTCGGGATTTCTTCCTGATCCGCCAGAGCCTTATGTAATCGGTTACCCGGGGGGTCGGGTGGACGTGCCCCGTCCCTGAAGACCTCTATGTCACGGCCCTCAAAATGTGAGATTGTAGGTTCCTTGCCGCCGTTAGAGACTCAGGTCAGCCAGGCATATCATGATACCAGGCTTCACCCTCTGAGTGGCTCAAACCGCAGCCATCATTGACCACAGGGCGCTTACGTCGTTCGACGCAACGCAGACATGGTCCTCTCTCACGATGTCTTGATACTTCCGGGGGATACCCCAGCTTGCGCTCTCAGCGCAACTTTCCATAAGCTGTGGTACCAATTACGGTAGCCTACACTGTGGTAGTTGAGTAACTCGCCTCGTCTTAAATTCACCTAAACTAATTGGAGGACTATCGGATCCGGCATCGCTTTCGTAACACTTATGGGGAGATAAACTTTCTTGAATAGCCGTGCCGCAATCTTATACGATGCCCAGTGTAGGCGCTTCGTAAGGTCTTTTTTGGCCCTGACGCCCCAATGTGACAACTCTTCATTCGTGCGAACGGTCCCGGCGTGCTACATCACAGAGCATAGTGGGCACAGGCCCAAACTTTGAGTCTCTCTGGTTCCGTCTACTTCGGCGGAGAAAGTTGCGAAACGAAATAAACAGAACTTTCCAGTCCTCAGAATTTATGCCTCCTTCGAAGCCCTATGTCATTCTTATTGTACGTCAATAAGACGCGAGGCCCCAGATGGTCCAGATGCACTCCGAAGCGCGGTGCTGTACGACAATAGATCATAGAGCAGCGATGGAAGTCCCAGTCCTGTTTCCAGGCTCAAGTCGCGCTGTGGATCCTAGTTACGAGCGTCTGTCTTTTACATTTTTACGAGTCTGCTTAATTGTCGCACTAAAGTAACTCTAATTGTCCGCCAGGGGAACGCAGAGCTTCTCGTCCCGAACAACCTATATGCATGGCAACTAAGCTCGTGTTCCACGGATGCCTGGATTGATTCTGATTCGGTGTTGTAATTCCTTACACGCATAAGGGCATGGGACGGTATTATTTAGACCCGAGAGGCGAACTATGGCTTGAGGCGCCGATATAATAGTTCCCCACTCTGCCTGGTTGCCTTGTCGCCGGAACTCTCGTCTTTACTAGTTACTCGCATTGATCGGTTGGACCGGTGACCTGACCTTCTCTGAGACCACCCGATTACGCTTCGGCTGTGTTCCGCCCGAAACTCCACGTAGTCGGCTTCATCGGAATCTCAACCTGAATATACCTACTAGGGCGACAGGAAGGAGCATCCGTTTTACTACCTGTGCCATGACTGTTTCGTGGTGCAAGGATGTGTCAGGAGTCCGGATAACTCGACCTCGGCATCGAGTCTCATTCCATGGCTATTGTCTGCCGACAAGCATGGCTTTGGTACAGGTAGAGTAGACGTGTTCCTAGGGGTCAGCTGACCACATTCGCGCCAAGGCGAAGGCAATGATCGCCCTGTTAATGTCCAGAGAGTTGGAAGACTAAGTCCTTGCCTGCGTCCCCGAAGCTGGCGAGGCGTGTGCAACGTGACCAAAGTTCGCTTTCCAAATGGAGAAAACAGAACTT'
    # al_out = debruijn_graph(k_arg, t_arg)
    # for k, v in al_out.items():
    #     print(k + ' -> ' + ','.join(v))

    #######################################################################################################
    # De Bruijn Graph from kmers Challenge - Given a set of k-mers, generate the De Bruijn Graph as an
    # adjacency list
    #######################################################################################################

    # kmers_arg = """"""

    # kmers_arg = kmers_arg.strip().split('\n')
    # dbg_out = debruijn_graph_from_kmers(kmers_arg)
    # for k, v in dbg_out.items():
    #     print(k + ' -> ' + ','.join(v))

    #######################################################################################################
    # WEEK 2
    #######################################################################################################

    #######################################################################################################
    # Eulerian Cycle Challenge - Generate the Eulerian cycle given the adjacency list of an Eulerian
    # directed graph (Challenge Data Set can be found in 'test_datasets' directory)
    #######################################################################################################

    # al_arg = """[dataset removed from here because it's too large]"""

    # graph_arg = {}
    # al_arg = al_arg.strip().split('\n')
    # for item in al_arg:
    #     node_edges = item.split(' -> ')
    #     graph_arg[node_edges[0]] = node_edges[1].split(',')

    # c_out = eulerian_cycle(graph_arg)
    # c_print = ''
    # for i in range(len(c_out)):
    #     if i < len(c_out) - 1:
    #         c_print = c_print + c_out[i] + '->'
    #     else:
    #         c_print = c_print + c_out[i]

    # print(c_print)

    #######################################################################################################
    # Eulerian Path Challenge - Generate the Eulerian path given the adjacency list of an unbalanced
    # graph (Challenge Data Set can be found in 'test_datasets' directory)
    #######################################################################################################

    # graph = """[dataset removed from here because it's too large]"""

    # graph = graph.strip().split('\n')
    # graph_arg = {}
    # for edges in graph:
    #     nodes = edges.split(' -> ')
    #     graph_arg[nodes[0]] = nodes[1].split(',')

    # path_out = eulerian_path(graph_arg)
    # path_print = ''
    # for i in range(len(path_out)):
    #     path_print = path_print + path_out[i]
    #     if i < len(path_out) - 1:
    #         path_print = path_print + '->'

    # print(path_print)

    #######################################################################################################
    # String Reconstruction Problem - Reconstruct the Genome given a list of k-mers
    # (Challenge Data Set can be found in 'test_datasets' directory)
    #######################################################################################################

    # k_arg = 25
    # p_arg = """[dataset removed from here because it's too large]"""

    # p_arg = p_arg.strip().split('\n')

    # sys.setrecursionlimit(2500)
    # s_out = string_reconstruction(k_arg, p_arg)
    # sys.setrecursionlimit(1000)

    # print(s_out)

    #######################################################################################################
    # k-Universal Circular String Problem - given an integer k, Generate the string that contains all possible
    # binary k-mers
    #######################################################################################################

    # k_arg = 9
    # universal_string_out = k_universal_circular_string(k_arg)
    # print(universal_string_out)

    #######################################################################################################
    # String Reconstruction from Read-Pairs Challenge - Given a collection of (k,d)-mers, assemble the
    # genome (Challenge Data Set can be found in 'test_datasets' directory)
    #######################################################################################################

    # k_arg = 50
    # d_arg = 200
    # gapped_patterns_arg = "[dataset removed from here because it's too large]"
    # gapped_patterns_arg = gapped_patterns_arg.strip().split('\n')

    # sys.setrecursionlimit(6000) # need to increase max recursion limit for large datasets
    # s_out = string_reconstruction_from_read_pairs(k_arg, d_arg, gapped_patterns_arg)
    # sys.setrecursionlimit(1000)
    # print(s_out)

    #######################################################################################################
    # Contig Generation Challenge - Given a collection of k-mers, output all contigs that can be constructed
    # (Challenge Data Set can be found in 'test_datasets' directory)
    #######################################################################################################

    # kmers_arg = """[dataset removed from here because it's too large]"""
    # kmers_arg = kmers_arg.strip().split('\n')

    # c = contig_generation(kmers_arg)
    # [print(s) for s in c]
