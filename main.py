import random
import copy
import sys
import itertools
import pprint as pp


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


#######################################################################################################
# WEEK 3
#######################################################################################################

def reverse_complement(pattern):
    """ Helper function for peptide_encoding() function. Takes in a pattern (String) and returns its reverse
     complement (String). """

    complements = {
        'A': 'U',
        'U': 'A',
        'G': 'C',
        'C': 'G'
    }
    pattern_rc = ""
    pattern_rev = list(reversed(pattern))
    for c in pattern_rev:
        pattern_rc = pattern_rc + complements[c]

    return pattern_rc


def translate(pattern, code):
    """ Helper function for peptide_encoding() function. Takes in a pattern (String) and genetic code (Dict) and returns
     the corresponding amino acid sequence (String). """

    aa = ""
    for i in range(0, len(pattern), 3):
        aa = aa + code[pattern[i:i + 3]]

    return aa


def peptide_encoding(dna, aa, code):
    """ Takes in a dna sequence, dna (String), an amino acid sequence, aa (String), and the genetic code, code (Dict)
     and outputs all substrings of the given dna sequence that code for the given amino acid sequence. It works by
     moving through the transcripted DNA sequence, i.e. the RNA sequence, using a sliding window of size 3 times
     the given amino acid sequence length (since 3 nucleotides encode every 1 amino acid). It then translates the RNA
     sequence in each window along with its reverse complements and checks the resulting amino acid sequence against
     the given amino acid sequence, capturing any matches. """

    dna_sequences = []
    rna = dna.replace('T', 'U')
    aa_len = len(aa)
    window = aa_len * 3
    for i in range(len(rna) - window + 1):
        rna_current = rna[i: i + window]
        rna_rc = reverse_complement(rna_current)
        peptide0 = translate(rna_current, code)
        peptide1 = translate(rna_rc, code)
        if peptide0 == aa:
            dna_sequences.append(rna_current.replace('U', 'T'))
        if peptide1 == aa:
            dna_sequences.append(rna_current.replace('U', 'T'))

    return dna_sequences


def prefix_masses(peptide, alphabet):
    """ Helper function for linear_spectrum() and cyclospectrum() functions. It takes in a peptide (String) and outputs
    a collection of prefix masses (List of integers). A prefix mass is a mass that can be used to build up a target mass that is
    used in a mass spectrum. For example, to arrive at the mass for the 'QE' portion of the peptide 'NQEL', one could
    use the mass for 'NQE' (371) and subtract the mass for 'N' alone (114). Here, 371 and 114 make up the prefix masses.
    These prefix masses are used to build a peptide's theoretical mass spectrum in this way.
    """

    # alphabet = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    # alphabet = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    # alphabet = [n for n in range(57, 201)]
    # aa_mass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    #            'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186
    #            }
    # aa_mass = {'AA' + str(i): i for i in range(57, 201)}
    # print(alphabet)
    # print(aa_mass)

    prefix_mass = [0]
    for i in range(len(peptide)):
        for s in alphabet:
            if s == peptide[i]:
                # el = prefix_mass[i] + aa_mass[s]
                el = prefix_mass[i] + s
                prefix_mass.append(el)
                break
    return prefix_mass


def linear_spectrum(peptide, alphabet):
    """ Takes in a peptide amino acid sequence (String) and outputs the linear theoretical spectrum (List of
    integers). It works by building a list of prefix masses, using this to generate masses for all possible
    combinations of the amino acids that comprise the peptide sequence, and capturing each in a numerically ordered
    list. """

    prefix_mass_list = prefix_masses(peptide, alphabet)
    linear_spectrum_list = [0]
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            linear_spectrum_list.append(prefix_mass_list[j] - prefix_mass_list[i])

    return sorted(linear_spectrum_list)


def cyclospectrum(peptide, alphabet):
    """ Takes in a peptide amino acid sequence (String) and outputs the theoretical cyclospectrum (List of integers).
    This works in the same way as the linear_spectrum() function except that it takes into account combinations of
    amino acids that occur when the sequence "wraps around". For example, "LN" occurs in the cyclospectrum of "NQEL"
    but doesn't occur in its linear spectrum. """

    prefix_mass = prefix_masses(peptide, alphabet)
    peptide_mass = prefix_mass[len(peptide)]
    cyclic_spectrum_list = [0]
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            cyclic_spectrum_list.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(peptide):
                cyclic_spectrum_list.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))

    return sorted(cyclic_spectrum_list)


def peptides_to_masses(peptide):
    """ Helper function for cyclopeptide_sequencing() function. It takes in a peptide (String) and outputs the
    sequence of masses (List of integers) that correspond to the constituent amino acids. """

    # aa_mass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    #            'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186
    #            }
    aa_mass = {'AA' + str(i): i for i in range(57, 201)}

    masses = []
    for nuc in peptide:
        masses.append(aa_mass[nuc])

    return masses


def peptide_inconsistent(peptide, spectrum):
    """ Helper function for cyclopeptide_sequencing() function. It takes in a peptide amino acid sequence (String)
    as well as a given experimental mass cyclospectrum (List of integers) and checks to see whether there is any
    mass in the theoretical linear spectrum of the given peptide that is not in the given experimental spectrum,
    returning True or False accordingly. """

    cs = linear_spectrum(peptide)
    for aa in cs:
        if aa not in spectrum:
            return True

    return False


def expand_peptides(candidate_peptides, alphabet):
    """ Helper function for cyclopeptide_sequencing() function. Takes in a collection of candidate peptides (List)
    that represent the peptides that are consistent with the given experimental spectrum. It outputs a collection of
    expanded peptides (List) which is obtained by appending every possible amino acid to each of the peptides in the
    list of candidates. This discards candidate peptides that have been retained from the previous iteration. At
    first I didn't discard these at each iteration and ended up with a never ending while loop in
    cyclopeptide_sequencing() because the list of candidate peptides was never empty. """

    # alphabet = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    # alphabet = [[57], [71], [87], [97], [99], [101], [103], [113], [114], [115], [128], [129], [131], [137], [147], [156], [163], [186]]
    # alphabet = [[n] for n in range(57, 201)]
    alphabet = [[n] for n in alphabet]

    if len(candidate_peptides) == 1:
        return candidate_peptides + alphabet

    relevant_candidates = []
    k = len(candidate_peptides[len(candidate_peptides) - 1])
    for candidate in candidate_peptides:
        if len(candidate) == k and candidate[0] != '':
            relevant_candidates.append(candidate)

    new_candidates = []

    for mass in relevant_candidates:
        for nuc in alphabet:
            new_peptide = mass + nuc
            new_candidates.append(new_peptide)

    return new_candidates


def cyclopeptide_sequencing(spectrum):
    """ Takes in an experimental mass spectrum (sorted List of integers) and outputs all peptide amino acid sequences
    that satisfy the given spectrum. It works by employing a branch-and-bound algorithm which starts looking for
    peptides at very short sequences. It expands these sequences (branch phase) by exploring all possible extensions
    of the current sequences. It then checks these expanded sequences against the given experimental spectrum and
    removes any sequences that are inconsistent (there is an amino acid in the expanded sequence that is not found
    in the experimental spectrum) (bound phase). If not for the bound phase, the list of candidates to explore would
    grow exponentially. *Note* my implementation explores all possible expansions of each candidate peptide. A better,
    faster implementation would explore only those amino acid masses that are found in the experimental spectrum to
    begin with. """

    alphabet = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    candidate_peptides = [""]
    final_peptides = []
    while len(candidate_peptides) > 0:
        candidate_peptides = expand_peptides(candidate_peptides, alphabet)
        retained_candidates = []
        for i in range(len(candidate_peptides)):
            current_peptide = candidate_peptides[i]
            if len(current_peptide) == 0:
                continue
            current_peptide_masses = peptides_to_masses(current_peptide)
            if sum(current_peptide_masses) == max(spectrum):
                if cyclospectrum(current_peptide) == spectrum and current_peptide not in final_peptides:
                    final_peptides.append(current_peptide)
                continue
            elif peptide_inconsistent(current_peptide, spectrum):
                continue
            retained_candidates.append(current_peptide)
        candidate_peptides = retained_candidates
    return [peptides_to_masses(peptide) for peptide in final_peptides]


#######################################################################################################
# WEEK 4
#######################################################################################################

def linear_peptide_scoring(peptide, spectrum, alphabet):
    """ Takes in a peptide (List of integers), spectrum (List of integers), and alphabet (List of integers) and returns
    how consistent the peptide is with the given spectrum. The scoring process involves counting the number of
    occurrences of masses in the spectrum, generating the theoretical linear spectrum for the given peptide, and
    comparing the two. For every mass that is shared between the two spectra, the score increases by one. In order
    for this to work properly, the given spectrum must be sorted from smallest mass to largest mass. """

    score = 0
    theoretical_counts = {}
    experimental_counts = {}
    for mass in spectrum:
        if str(mass) not in experimental_counts:
            experimental_counts[str(mass)] = 1
        else:
            experimental_counts[str(mass)] = experimental_counts[str(mass)] + 1

    cs = linear_spectrum(peptide, alphabet)
    for mass in cs:
        if str(mass) not in theoretical_counts:
            theoretical_counts[str(mass)] = 1
        else:
            theoretical_counts[str(mass)] = theoretical_counts[str(mass)] + 1

    for k in theoretical_counts.keys():
        if k in experimental_counts:
            score = score + min([theoretical_counts[k], experimental_counts[k]])

    return score


def trim(leaderboard, spectrum, n, alphabet):
    """ Takes in leaderboard (List of Lists of integers), spectrum (List of integers), n (integer), and alphabet
     (List of integers) and returns the top n peptides that are most consistent with the given spectrum, including
     ties. It involves sorting the leaderboard from highest score to lowest score, checking past nth item in the
     collection to see if there are peptides with equally high scores and retaining. It then returns a List of Lists of
     integers as the new leaderboard.
     """

    leaderboard_dict = {}
    for j in range(len(leaderboard)):
        current_peptide = leaderboard[j]
        peptide_index = '-'.join(str(mass) for mass in current_peptide)
        leaderboard_dict[peptide_index] = linear_peptide_scoring(current_peptide, spectrum, alphabet)

    last_index = len(list(leaderboard_dict.items()))
    leaderboard_dict_sorted = {k: v for k, v in
                               sorted(leaderboard_dict.items(), key=lambda item: item[1], reverse=True)}
    leaderboard_dict_sorted_items_list = list(leaderboard_dict_sorted.items())

    for i in range(n, len(leaderboard_dict_sorted_items_list)):
        current_item = leaderboard_dict_sorted_items_list[i]
        if current_item[0] != '':
            current_peptide = [int(x) for x in current_item[0].split('-')]
        else:
            current_peptide = []
        nth_item = leaderboard_dict_sorted_items_list[n-1]
        if nth_item[0] != '':
            nth_peptide = [int(x) for x in nth_item[0].split('-')]
        if linear_peptide_scoring(current_peptide, spectrum, alphabet) < linear_peptide_scoring(nth_peptide, spectrum, alphabet):
            last_index = i
            break
    leaderboard_dict_sorted = {k: v for k, v in list(leaderboard_dict_sorted.items())[:last_index]}
    new_leaderboard = []
    for k, v in leaderboard_dict_sorted.items():
        if k != '':
            peptide = [int(x) for x in k.split('-')]
        else:
            peptide = []

        new_leaderboard.append(peptide)

    # return list(leaderboard_dict_sorted.keys())
    return new_leaderboard


def cyclopeptide_scoring(peptide, spectrum, alphabet):
    """ Takes in an amino acid sequence, peptide (List of integers), and a theoretical spectrum, spectrum (List of
    integers), and outputs the score of peptide against spectrum (how consistent the peptide is with the spectrum).
    It works in a similar way to linear_peptide_scoring() except that the theoretrical cyclospectrum is generated
    instead of the theoretical linear spectrum. For example given a peptide NQEL, the masses for LN are included
    in the cyclospectrum and subsequently scored by this function, whereas LN would not be generated and scored by this
    function. This is because LN is ontly obtained by wrapping around NQEL."""

    score = 0
    theoretical_counts = {}
    experimental_counts = {}
    for mass in spectrum:
        if str(mass) not in experimental_counts:
            experimental_counts[str(mass)] = 1
        else:
            experimental_counts[str(mass)] = experimental_counts[str(mass)] + 1

    cs = cyclospectrum(peptide, alphabet)
    for mass in cs:
        if str(mass) not in theoretical_counts:
            theoretical_counts[str(mass)] = 1
        else:
            theoretical_counts[str(mass)] = theoretical_counts[str(mass)] + 1

    for k in theoretical_counts.keys():
        if k in experimental_counts:
            score = score + min([theoretical_counts[k], experimental_counts[k]])

    return score


def leaderboard_cyclopeptide_sequencing(spectrum, n, alphabet):
    """ Takes in spectrum (List of integers), n (integer), and alphabet (List of integers) and returns the peptide
    that is most consistent with the given spectrum. It works by using a branch and bound algorithm where the branch
    step involves expanding the candidate peptides by all of the peptides given in the alphabet. The alphabet
    represents all of the acceptable amino acids. Initially this was all possible amino acids in the mass range of 57
    Da to 200 Da but was later modified to include only those that were most frequently appewaring in the convolution
    of the given experimental spectrum (see below). The bound step involves scoring each of the peptides retained at
    each iteration. The top n scoring peptides are kept for the subsequent interation while the rest are discarded (see
    trim function). **Note** commented code was part of a minor exercise that aksed to find all peptides with max score.
    """

    # best_peptides = []
    leaderboard_peptides = [[]]
    leader_peptide = []
    while len(leaderboard_peptides) > 0:
        leaderboard_peptides = expand_peptides(leaderboard_peptides, alphabet)
        retained_candidates = []
        for i in range(len(leaderboard_peptides)):
            current_peptide = leaderboard_peptides[i]
            if sum(current_peptide) == max(spectrum):
                cscore = cyclopeptide_scoring(current_peptide, spectrum, alphabet)
                # if cscore == 82:
                #     best_peptides.append(current_peptide)
                if cscore > cyclopeptide_scoring(leader_peptide, spectrum, alphabet):
                    leader_peptide = current_peptide
                retained_candidates.append(current_peptide)
            elif sum(current_peptide) < max(spectrum):
                retained_candidates.append(current_peptide)
        leaderboard_peptides = retained_candidates
        leaderboard_peptides = trim(leaderboard_peptides, spectrum, n, alphabet)

    # This code was used as part of a minor exercise that asked to find all peptides with max score
    # print(len(best_peptides))
    # c_out = ''
    # for i in range(len(best_peptides)):
    #     masses = best_peptides[i]
    #     c_intermediate = ''
    #     for j in range(len(masses)):
    #         c_intermediate = c_intermediate + str(masses[j])
    #         if j < len(masses) - 1:
    #             c_intermediate = c_intermediate + '-'
    #     c_out = c_out + c_intermediate
    #     if i < len(best_peptides) - 1:
    #         c_out = c_out + ' '
    # print(c_out)

    # max_score = cyclopeptide_scoring(leader_peptide, spectrum)
    # print(max_score)

    return leader_peptide


def get_top_elements(p, collection):
    """ Helper function for convolution_cyclopeptide_sequencing() function. Takes the top p elements of the given
    SORTED collection (Dict) and returns them in a new Dict. The sorted collection in this case was a convolution (Dict
    key=mass (Str), val=count (integer)) sorted from most frequently occurring mass to least frequently occurring mass.
    """

    last_index = len(list(collection.items())) - 1
    for i in range(p, len(list(collection.items()))):
        current_aa = list(collection.items())[i]
        if current_aa[1] < list(collection.items())[p-1][1]:
            last_index = i
            break

    top_elements = dict(list(collection.items())[0:last_index])
    return top_elements


def spectral_convolution(spectrum):
    """ Takes in spectrum (List of integers) and outputs the convolution of the spectrum (List of integers). The
    convolution of a spectrum is the collection of masses that are obtained by subtracting all masses from one
    another in the given spectrum and retaining all resulting masses greater than 0. THis allows us to determine the
    amino acid composition of the underlying peptide. The convolution is then used as the collection of acceptable
    amino acids that should be explored by the branch and bound algorithm. In effect, this simultaneously reduces the
    amount of candidates that need to be explored by the algorithm and expands the candidate pool to include valid
    masses. """

    spectrum = sorted(spectrum)

    convolution = []
    for i in range(len(spectrum) - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            res = spectrum[i] - spectrum[j]
            if res != 0:
                convolution.append(spectrum[i] - spectrum[j])

    return convolution


def rank_convolution(convolution):
    """ Helper function for convolution_cyclopeptide_sequencing() function. Takes in a convolution (Dict key=mass (int)
     val=count (int)) and returns the convolution sorted by most frequently occurring masses to least frequently
     occurring masses. It also restricts the convolution to only include masses that are in the range 57 to 200
     (inclusive). Each mass represents an amino acid. """

    counts = {}
    for aa in convolution:
        if 56 < aa < 201:
            if str(aa) not in counts:
                counts[str(aa)] = 1
            else:
                counts[str(aa)] += 1

    counts_sorted = {k: v for k, v in sorted(counts.items(), key=lambda item: item[1], reverse=True)}

    return counts_sorted


def convolution_cyclopeptide_sequencing(m, n, spectrum):
    """ Takes in integers m and n and spectrum (List of integers) and returns the peptide (List of integers) that is
    most consistent with the given spectrum. It works in a similar way to leaderboard_cyclopeptide_sequencing() and
    in fact calls this function. However, in this method, the top m elements of the convolution of the given spectrum
    is used to generate a shorter list of acceptable masses to explore during leaderboard_cyclopeptide_sequencing(). """

    convolution = spectral_convolution(spectrum)
    ranked_convolution = rank_convolution(convolution)
    top_aa = get_top_elements(m, ranked_convolution)
    alphabet = []
    for aa in top_aa.items():
        alphabet.append(int(aa[0]))

    leader_peptide = leaderboard_cyclopeptide_sequencing(spectrum, n, alphabet)

    return leader_peptide


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

    #######################################################################################################
    # WEEK 3
    #######################################################################################################

    #######################################################################################################
    # Translation Challenge - Given an RNA pattern, output the corresponding amino acid sequence
    #######################################################################################################

    genetic_code = """
AAA K
AAC N
AAG K
AAU N
ACA T
ACC T
ACG T
ACU T
AGA R
AGC S
AGG R
AGU S
AUA I
AUC I
AUG M
AUU I
CAA Q
CAC H
CAG Q
CAU H
CCA P
CCC P
CCG P
CCU P
CGA R
CGC R
CGG R
CGU R
CUA L
CUC L
CUG L
CUU L
GAA E
GAC D
GAG E
GAU D
GCA A
GCC A
GCG A
GCU A
GGA G
GGC G
GGG G
GGU G
GUA V
GUC V
GUG V
GUU V
UAA 
UAC Y
UAG 
UAU Y
UCA S
UCC S
UCG S
UCU S
UGA 
UGC C
UGG W
UGU C
UUA L
UUC F
UUG L
UUU F 
"""
    genetic_code_dict = {}
    genetic_code_list = genetic_code.strip().split('\n')
    for codon in genetic_code_list:
        current_codon = codon.split(' ')
        genetic_code_dict[current_codon[0]] = current_codon[1]

    # p_arg = "AUGAGCCAACGCACGGCACCUGACGUCACGUCUACUAAACUCUUACCAGUUCGGUUGUAUACGACCAGGCUUUUAAACAAGACUUUAUUGAAUGGCCGCAACUAUAAGCGGACGAAGAAAGGUCGGAGUCCGAACACUGCGUCAAAUCGCCGGUACCCGCGAGUAGGUAUGUAUCAUGGACCAUCCAAAUGCGCCGACCCAAAUAGGGUUAAUGACACCUAUGUCAAAACUGAUCUUGGCUUAUUCCCGCAGACAGGUGACACGUUUUUGACAGACGCACAGGUACAGACUGGAAUCACGGGCGUGGGUGCCACAAGCUGCUGGUGCAAGCUUGUUCGAUAUUAUAAUCUGUACCAAUAUACCGUGAUUGCUCGCCGUACUAUGCCUCUACCUCUGAAUCAAAGAUCAAUUCAACCCAGUAAAGUCUGUUUGCACGGGAAGCAUAUUAGCACCGGGCCGACGCAGUCUUGCUUGCCAGGACUACAUGGACGCAGGUAUUGGGCCAGACGCUUAGCAUAUACCCAGCUCCGGACCCCGCCCCUAGUUAAUUUUAGAGACACGUCGAGGGACGACUUGGGCUGGCCCCCAUUACUUUGUAGCACACUGCGCUGGCGGCGUAUAGGGUGUAACGCGAGAGAGAUCAUUUCCGUGGAGCUGCUAGAAUUAGAGCACGUUGUACACCGCCCGCCUAACAUCAAGCCCCGUGGGAAUCGCCUUGCAGUGAUGUACAUCCAUGAUGAGAUUCAUCUUACACCCCAACGAUCUAAUAGAAGUUUCAGGCAAAGCGUAGCGGUCAAGCUCUGGGGUAGACGGGUCCGCGUUGCUCGCACUCUAUUGCGCGAACCGUUCUUCCUGGACAUGGAGUUGUUGGGAAGUAGAUACAUGGUACUUUCUGAGACAUUACGAACCAGCGUCGUGAGCUUAAGCGGAAGGGAAGAGUAUCGGCACCAGCUUUUCAAGAGCCUGACUGUAUUCCGUGCUAGGCGUUUGGCCGGCUCGCAUGCAAGUAUUAACGGGUCAGAUGCCUACUGUACAUCUUCACUGAUGCGAUGGUUUAACCGAAAUCGAAUAAUCAUUAGGCGUUAUGAUUUGCGGAGCACAAUCCAAAUAACUAUUAUUCAAGAGGGGCGCAAAAUACUCUGUUAUCAACAUUACGGGUGGUCUAGGAACGGCUGUGUUGCUUCCUCUGGUGAUUCCACCAAGAAAUUCCGCCGUUAUCAGCGAUGGGUUGGGUCCAAAAUCGUUGUCGGCGCGAUUCAUCCAAUCGCGGGACCCCGGGUUACCAAAUUUGUAUUGCACAAUGGAUGCGGUAGCCGCGCUGAGUACCCAAUGGGCCUGUUUUUGCGUCCUCCUUUUUACUGUCUGAUAAUGCAACACCGACAGUAUGACAGACGAUAUGCGCAAGAGGAUAAUUGUUGGACUCGGCUUACGAGGGGUAGCACGCUUUUUGAAGGUAGCACAGUGGAAACCGCGUACGCGCUGGACGGGCAGCUUUCUACACUAGCUGCACAGGGGCUAGCGCAGCUUCGGGAAGCUGUCUGGGCAAAUCGGAAAUGCAUGUUAAGGGGACACAACUGCGAUGAAGCAAUCCAUUUUAAGAGCGCAUCAUUUCUGGACACGGUUCCCGGAGGAUGUAGUCGAGGACACACACGAAGCACAGCACGUAUUCAGAAACGAUAUUUACAUGUACUAUGUCCUGCAGUUAGUGAUGGAUGUGUCACUGGGGACACGACGGUCUACUUCCAAUCCCAGGGCAACCUGACGGCUACCUUCCACUUUGGCGGUGGUAAUCGUACUAAAACUCGUUGGGUCCACUGCUAUACAAGUUUCAACAAUUCUACAGACUCUGCGUACCCUGUCCUGUACGCGGCUUGGUCCGUGCUAAGCUAUUCGAAUAAGGCUGAUCGUAUCACGAAAAGUUCACUGGCAGGAUCGUACGGCACUGGAUUAACCCGCGGUAGUGACCGCCAUCGAAGUAGCCACUCCAGCUUGGGGGCACGUCUCGAAUUUCACAGCAGCCCCAGGCGGAGUCGUUUUGAUCGGUGCCCUAGCCUGAUGUUUCACUCGAUCAAACCAUUCCCGGCAAGUUCCAUCGUCGUAUUCUUCUUCGCGCGGAUCGCUCUAAAUAUCAUCGCUUUGCUACUGGCCCGUGAUUAUGUCGACAUAGGACCAUGUAAUGCCAGUAGGUUAGGCGAGCCUACGAUUUUUAUUGCCACUCAGUGCGACAGGAACUGUCCAACGGGCUUCUUUUCGAGCGUGAUUCGCUUGGAGCAUUGCUCAGUGGGCUUCCUAUUUAAGGAGGAUACUAUGCGUCUGGAUCGUAGCUUCGGGGGACCACCCUGUUUAGUUCAUAAGCAAGAGACGAUCCCCCCACCAAGAUGGGAGAUGAAGCUAACCUUACUCUACCGUGCGGGUCCGGAGGCGUAUUCCCUUUUCGUCAAAGAGCCCCGGUGCUCGAUCUGCUACAAACUCUACCACUCUUGCGACACUGAAGGCAAAUACAAUAGAGGUGUGAAACCGCUUUAUGGUUCUUCCGUGUGGUCUCAAACGGGUGACUUUGAAGUCCUACAGUCUAAACCUCAAAGGACACAGACCCUGGCUCAGGCACCCACCCGCGCUCAGAGUGCUAUAUUCAUGUGGAUUUACGUCUGCCAAAUAGUGGGGUGGUUGAACCAAUACCCGCCCGAAAGCAGCUGCGCAUGCUUCGGAAGCUCAAUAGGCAGAGUGCGUCUUAUCCUCCUGUACUCGGAUGUCUUGCCCAAACAGAAUUUCAAAAAGCAGAUCUAUCAGCUAUCCAAUCUUGAUGUGCGUGUUGUAGCACGGCAUGUAGUCCACCACUUUCAUCGAUAUUACAGAACAUGCAGGGUUUACAAACCACGCACGUGUCCCUAUACCGCAUUCAGACGCUGUUUCGCUCCCCAAAUCGGGACAGUAUCUGCAUGUAAUCUGCGAAUUCCGGAACUGCAGACAAGAUGGGGCCUGCGCCCAGAUCUACCUUUGCUAGCGCCAAUGUGCCCAGGCAGCUAUUGGAGGAAUGGGAUCGUCCUCGUAGACACCCCCCAAUAUGCUAGGCCGGCGCGAGUAGCUGUAAGCAUUGGUGGUUUCCGUGCGCCCCGUCUUCGGGGACCAAGAUUGGGAUCAGAGGGUGACCGACGCAAAGGCAUUAACUACGGGAUGCAAUUCAGUAUAGGAUGUGCUUCGCGCUGGUCCGUGCGCCUUCCUUCAAUCCUUGCGAUUUUCCUGGUCAUACUUCUGCGCGAGUUUAAAAGCAAUCUACGUGUCCACCACCGCUACGGUAAUACUACCCCAUCAAACGAGCACCUUGUAAUGAAGAAUGAGGGUGCUAGAACACCUCCUAUGCUCAUAGCAAAUUGGCCAAUCGUUGCAUGCUCAUUCGGCUCAAUCAGCCCAGAUGCUGUUGGGGGUAUUAGAGGCACUUUCGGGUUGGCACAUGCUAGCUCUUGGCGAAUGCGCGUCGCGGUCAGCUGUCGUGUCGGAGGUAUUAUCGUGCUUUUAAGGAAAUUUGUCCUAGCUGGGCUUGAUUUGUGUGCCCGAAAACAUAGAGUCGAUGCUCAACCGGGAACAUCGGCGAUGAUAUGUCUGGGCACCCGUGACUUUCGUACUGCAAAAUGUUAUGGGUGCCGAUCCCUAGUGAAGCACUACACGCGCCCUUUAACAGAUCAUUAUCCUUUUGAGGUUCGAAGGCUUGUGUACGUUAGUCCCAACACCAUGUUCAGACUCAACUGCGAGAUGGUUUCUGAAGGGAGGAGGACGUGCCAAAGGCCCAAUUACCACGGGCUGAUGCACCGUCCACGUCCUUACAUGCAGCCCUCGACUGGGUCCAGUGUGUCUACCUUACCUCAGAAUUUAGCUUUGCAUCGUAGUGUGCGUCUCGAGGAGACAGUGCUGCGAAUUAUUUGCAAUUUGCUCGCACUGUCACCAUGGUAUGCGUGGUUCGGACCCACAGCAGUUCUGACUGAGAUACCGUUGCCUUGGAUUAACUUGAAGGUGCCGCUCGAGAUUAGAUCAAACGUUAUUGACGCGCGCUAUCCUACGCAUGCAGCGUUGAGACUCAAGAAUUCAAGUGUACUAAAGGCCAAGAUCGAUUCACAAUACAAGCUCUGUACGCUAAGAAUACAAAAACUUACCACGCUCACUGCUACAGGGGUGCUGGGCGGAGCUACGAUCAUGUCGGCUCGAAUGCUCUACCCAUGUCUAUCUAUACCACUUUGCCGAGAAAAUACUCCCUUGUUUUUCACCUUGGGGAAUAUAGUGUCACGACUCUAUCGGUCCAGCAACAGCCGUGAUUUGCCCGUUUCUGUGGACGCGAAAGCAGCUUACGUUUGGGCUACGCUUCCCCCGCCGAUGAAUUCCCGUCUUAGGCCCAUAGUGGUUUUUGAGCAGCGCGUGUGCCAGCACACCACUCUCUCGAUACGCGGUAGAACAUGUAUCCGAACAUGGAAUCCACAUGGUGGCAAAAUGCCACCCCUUUAUUCCAGGAUAAAGCCUUUCGCGUUUGGAUUAUCGGGUGGGCGCUCUUCUGAAUGCGACAAGUGCCGUCAGAAUAUGAGCCUUCCUGGCACCGUCCGCGGGGCAGGACUGCCACCGAAGAUCGAGAAGCUAAGACAAGAAGGCUGCCCUUUUAUAUCUCUCCCUACUGUGGGUGUAUGUGUGAUACUUGUAAACUGGUACGCUACCUGCUGGAUAUGGUGCGCACAGACGCGGGCCGUCGGAGGGAGCGACUCGCAACUAGUUAGCUGUGCCACAAACGUUUUUGGGGCAGUACCCCGUAAAUGCGAGAAUCGGGAGGUGUUAUUACAACCGAACCCUCGCACGCCUAUCUGUAGGGAGAGUGCCAGAUGGGUCGAUCGAUGUGUGAAGCGCCACCUGCGAUUCGUCCCGGCAAUUUGUAAGGCUAGGAAUAGCCAUAGGGGCCACAUAGCUGGUGCACAGCUAUUGGCGCGUGUAUGGCUGGACGAUAAUGGGUGCGGAUCUCAUACGGCCAAAACCUCGGAGAUAGAAUGGUUGCCUUCCUGCAAUCACUUUCCCACAGCAUUACGUAACGUAGUAUACCAGGCUCAACCAACGCCGGCGCCAUGUUCGGGGGACCUAUAUUUCAUGCCGGCGGAUUUCGGCCAUCCACUGCUGGAGAUACAACACAGUAUUUCGACUCCUUCUUACCAAUGGCUGUUUGUCCUGAGUCCCGCGAGGGUUAGGAAACUCGAUAGCUUCCAGAACCGUCGGGCCCAAUCGCCGCAAUGCCGGUCCCUCUCGUUUCGCUUGCCUGUGUCCUAUAAGAGCACCAGCUCGAGUGCGCACAGUCACAUACACUGCUCGGUAACGCGUGCGCGUGGCACCAGUAGCCUCAACAGACCUCCAGCAGAUUCCUGGGAUCUUUCCGAUUCGGCAACUAACGAGCUCAGUUUGUCCUCCAGAGGUGGGAUCCGUUCUCUAUGGCCACGUUUUUUAGGAGCUCCGGUAAAGCUACGCGAUAAUCGCACAUGCCCCGCCCCUCACAUAUGGUCGCCGAGAGAGACGGUUCCUCCCUGGCAUGAUGUUACUCCUCGGGCGGGAACUAGUAGCGAGGCGCCGGCCGCCGCUAAAUGGAGGUGUCUACCCGUUCUCCAGCUAAGGAAUGAAAGUAACACCGUUAGUAAUUCCGUCGGACCCCCCGGACGUCCGGAUGUGAGAUCAUUCGCCUCAAUACUCUCGAGGACGCACACGCUGCCGGCAAUGUUUCAUUUCAGCAACUGCACUGGAAUAACCUCAGAAGCAGGGUCUGUUAACUCUGGCACUCUUCAACUUACAUACCGAAUUGUCAUUGCACGUCAAGACGCGCCUAUGUCCACAAAGUCUAAGUUGCACACAAUCCGCAAGUACAAUCUGAUCGGACUUCGGAGCCGGACGCGUUCUGCGAACCAGCUGGGUACACAAACCACGAAAGCUAAUAAAGGAGUGGUCCCCGACUCUAACUCAAAGCAGGUGUGCGGAUCCCUGUUCCUGCGGGGCGGAGAUUUACGCAUUUGCGCAAUGAUGAAUACGAAGGUGCAGAUACAAUCACAGCAUAAUGCUACCGAAUUCCCCCAGCUGUUCAUUGGCCGCGGCAACAUGUUUUGUUACAGUAGAGAACCCAUAAGAGCCCCCCGUGAUGAAAUCUCGCUAAUGACUAAGGUGACGGGCGGGUCGCCGCCAGGGGGUGCACAACACAAUCGGAUAGGAUUCCACGCUAAUGCAGUGUUCGGCUCGAUUGGGGACCGAUCCCUGGCCUGCGCUGACGGGUUUGCUAAAGGAUGGGGCACAGAGAUGGAAAAACUUUCUCUUCGAAUUUUAGCAAAAUCUCGGACUCUCGUUGGGCACAUUCCUUUGUCAAAAGCCCAAGGUGGAUUCUCCACAUCUUGUAUUCCGUCUGUCCCGUCGCGUAGAUUGUUAGCUAGCUGUUCAUGCGGAUCUGCCUGCCUGGGAACCCUACGCUCUCGGCAUGCUCAGGGCCUUCGAUAUCCCCUCCAUCCCCCAAUCUGUGGGCGGUCAACAUUGGACGAGGUGACUGUGUUGGACGGAGUAGUGGUCUCCAGACCCCACGCGUACAAACGCAAUCAAUCAAAUACUCUGCCCCUAUCACCCGUUCGUGGAGUGCCUUAUCGGCCUCAUCGCGGAGGCCGAAGACUCUUAGUGAGCCAAAUUCAUAUGACUCGGUUUGCUACUCAUACACCCUGUAGCAUCGCUAGUCGGCCUGAGCCCAAGCUUCAGGUGGGUCUUAGAGAACCAGCCUCCUGCUUGAAUUUUACGGUAACGAAGGUCGUGUGCUCAGAAGCCACCCGUACGCACCCUAACUGCUCGGUACAAACUCUAAUAGCAGCCACUGCCAUGGACUGUCGUCCUAAUCCGUUAGCGGUGUUUCUCGUAUACUUAAAUCCUACUUCAAGUAGGGUAGCUUUCAUAUAUCUUAUACAAUGCGUCAGAUGCUGGCGUCGUUUUGCAUGGGCAUGCCUCCAUCACGUGUUCAUAAGAGCACUGUUGUUUUAUGUAAAGAGCAUGGAAGUUGGUAGAAUUACGACACUCUAUCCCUGUAACAUUCUGCGUCCCGCGCGCUCGAUUACCCGCGCUACCGCGCUCGAAUUUACUGGCCCUCGUGAACAUCCCGACUCGGAACACUCGCGGAGAAGAAUCGCACAAUUCUAUUUUAGGGAUGCCAUACGAUGCAGGGACCUGCAUCAAAUCAGUACAAAUUUAAACGCUACACCGAUAAUACUCUAUGAGAUACGCUUGACUAAAGUCAUACUUGGAGUCCGCUAUGGAGAAAUCUUCUCCCUACCUAACCUGCCAGGUGCCAGUAGCGCCAAAUUUUUAAGCGGAUACGAAGUACCGGUAGGGGAGGAAUAUAAAACGCCUAGUUACAAAGUCAGCACCUGCAUGGCUACCAUCUUGUUAUGCACACUGUACCAGCGAAUAGGACGCAGAAGUCCGCCCACCAUAACGCGUACUCUUGCAGCUCACAUUCUUCUCAUUCCGGGGCAUCCACGGUGGUCGGCCUAUAGUUCUCGUCAAAGGGGUAAUUUUUCCUCCAGCUAUAGUGUCGGUGAGACCGCAAACUCCCUGGUUGAGCUAUCCGCUCCAGUGUGUUGUUACGGUCGUCUAUACAUAAGCGUACCUCGACAGGGUACUUACGGCCACAUGCAAAAGCCGCAAUCGGAGCUAUGUCGUAUAUGUCGUUUCUGUUGUAUCAUUUACCUCUCCUUCUUGCAAGAUGCGUCUAGGAGCAUCUAUAAUCUAGAACGGAAUAAGGGACAGACGAUGAUUGAUGCCAUGGUAUGCCCUUGCCUGUACCCGGGAGGUGGUUCAGCGAUGGCGACUCCCGACCACAUGCCGAUGUUUACGCGCCUGGAGCGUACACUCCGCGACGUACUCUCUGGAACUAGCAAGCAACCUCUUAGCCCCAGACAUGGCACUGGAACGUGUAAUUGUUAUCUUUGCCCUUCUUUUUACCUACGCACGCCAAAAGCAUUUUUUCAAUCCCGACCCAAUGAUGGUAUAUGUGCAACGACGUAG"
    # aa_out = translate(p_arg, genetic_code_dict)
    # print(aa_out)

    #######################################################################################################
    # Find substrings of a genome encoding a given amino acid sequence
    #######################################################################################################

    # dna_arg = "CGCCACATCAGGAATCGCGCACGAAACGGCCTGAAACTCCGCAATGCTAGATTGGATCAAATATGAAAAGGTTCCAGGCGCAATTCTAGATTACCGTGTTCTGGCCATAGAAAACAACGGCTGCTTAGCGCTACACCGGGTTTAATGTCGATGCGTCGTCATTGAGTGAAAACATCATAGGCGCAGACTTGATACTCCCCTGCCCCGGACTTGAGAACGGGAGACTGCGACTTAACACCTTGCGCTAGGTCATTTCCAATCAAGGTCGTCTGGGGACCTGGAGCACTTAGGAGGTCCGCAGCCTAACTACATTCGCTGGGCAATCCCACATAATCTTTGTATATCGCCCCGATATTATACCACGAGAAGCTCCTCTCACTCACCGTCAGACACGTATTGAGAAGCAGGAGCCCGTGGTGAGTTAATTAGTTGAGGACTTTACAGTTGGCTCCTTTACCGTGCGACGAAGGCGTTTGGCCCCTGACAAGATCCCGCACAATTTATGCATCAGCCCACGACCTACAAATGGAACCGTAACTGTTCGCCTCCCATAATCAGGCTAACTCTGTTCGTCCACTAGGATTTGTTTTTATCAATTTTGTATTCGTGCTCGAGGCAGTGAACCCATAAGAGTATAGTGCGACATGGCTGATCCCCTACTCGGCTCACGATCTTGAGAGCAAAATTCAGGCGTATGTCCGCTTGAGCATCACTCTACGGTCAGATATTGGTACCAGTGTTGGTACTGTTTACGAATTCAGTCGTCGATCTCCCGAATCGTCAGAAAATCTCAGCGTATGTGGAGGGTACACCATCATCGTCTAGTGGCTGCAGGTTTATGCAGCTCGTAGATCACCGGTTGGGCATGGTCGCTGTACATATAACGCTCGATAAGGTCCTCTTACTGATGACTATGCAAGTTTCAATTCGCATGCCATGGTGTGGGGGCCCTACGTGCTACGATTGCATTTGATCTTGCCGTGATACGTTGAAAAACGGAACGATACATAACGATTAATATGGCACGCACCGGTAACCGAGGCATACGAGGGCTGATGCACATCGAGGGCTGATACAGAGATTATGAGGGATTAGCCTTTTGCGCTACATTCCGGCTCTGGCCAAAAGCAGAACAGCATTATTAGGTAAGATGTGCCGCCGGCACCTGCCGCAAGACCCGGGACGTGTACGATTTAACTCCTCCGCGGTATGTCATCCAGACAGGTGATTTCCGCCAGTGTATTATTCGGCCACTTACGAGCTGCCGGCCTCCCTCTATATACGGCGGACTCGGTAGCTCCGTCGCAGAAGTGCTAACAAGGGACGCGGTGCTAGCAGTGACCTAATGCGCATTCATCACACAAAGGCTCCAAGATCTTACAGCAACCTACCGTAGTTCTCCATGGGACAGAGGACAGATCATTAACTTATCCATTTGCAAACGACGCAAATCGCCTAGGCTACGCATTCGTAAGGGACGGGTTCAAAACGGTTCTATTCATCGAAGCAACATCGCGTAGGCTGTCGGGGAACCGCATAATGATAAAGTTAATTAGGCGGCATCCCTCATAACTTATGTATCTCCCCAAGAAACGACACCGTCACGGAATCCGGTCTGTAAATAGGAGGCGCAGTTCTCGCCGGCAACTACGTCGCTGTTTCGCAGACCAGAGGTCTTACCATTGATCGCGCAATACCTCACAACTTGTGTATTTCACCGCGGAACCTGAGAAATATCCCAGACCAGGGCCAATTGGCGGTTTGTTACTTTGAGGCGAAAATCGCACGTTCGTCCTAATTACCGCTTAGCACGTCACAGATCTCGTAATTACACGCAATCGCCTCGGGGCTGCCTTTGGTTAGGGTGAGGCCTGCGGCTTCGAGACGCAGGGGATGCTTCCCGAATTGTTCTACAGGGGGCACCCGGTTGGCATGTGCTATGTAATTGACGATAGCGTCTGGTATATGCACAGTTGTACCGTGGTCAAGCCAGTAACGAACTCCGATTTTGGTGCAGCCTATTGCCATCCACTATTTTCATCATAATCAGTCTGCTGGGTATGCACTAGGTAGTGAGGTCACGCCCCGTGCTGGTCGTTGTCCTTGTATCTGATCGTTTAACATGAAGGAAATACGACAGTTGCAACTCCAAAGCCTTGGAGAAATGTAGGCACTCTCTCTGCCGGGAGACTCATATAACGACCATAGCGCTGAGCGTTCTAAAACACGGTGTCCGTGAAACGCCCCCGCATGAACATTAATCATTCATTAGGCTAATGCATACAACGTCGCTCTTCTACGCATACGACTTCTCTACTGACAGTTGTCTATATCGAGCCTAGTTACCGCTACGTGAAGTACAGGTCCTGCCCTTGGGGAAGGGTGCTACGAGCTATGCCTCTTCGAGTAATGTGGTCAATGTAGTGGTTACTCAGGAATGCCAATGGCGCATGTTAGGGGCTCGTCCGGTGGTAAAAGTATGGGTACTGACCGAAGTCGGAGTCCGAGAGCTCGGAGTAGTGCGTGTTCGTGCAGGTTCTCGGTGAGATGCACAAGTTGTGAGGGATGCGGCTTCCCAGGATTCCACACAACCTGTGCATATCCCCTCGGCTGTCAAACACCGCACTCGGGACCTATTTGAGACCTCTGGAACCTGTACAAGAACAAAGACAGACGTACGCCTGTTTGCGTAAACCTGACGCTTGATGGGTGGTGGGTTGACTACATTATAGAGCAACCGCGCAGAATTCCAGTTTGGCTTATGTCAGGCCCTTTCACACCGAGTGGCCGATTTGAGGACAAATTTTGAATATACCTCATAATCTATGTATCTCCCCTAGAACGGTTAATCTCCGATTACGTTTCGCGGATTGCCAAAGCCAGTATGCCAGCGGTCAGTTCCGTGTAACCAAACCTAGCAATTGCTGTGCACGATGGAGAAGTTGACGTGTCAGTTCGTCCGACCCCATCAGGTGTAGGGATCATTCGGACAATGAGCGGTGCGGCTTAAGCAATGAGATTGCTTCCATAATCATCTCCTCCTTCACAGTCCACGAAACAACTCGTGTACGCCTGTCCAATATTAAAGCAATGTGGGTGTGAATTATCATTTAGAACACATGACATTATGTGATAGAGTTCTACAGCTAATCGGCTCTAACCACTAGAATATCTGATACATTAGTCTTGTAAAGTAAAAGTCGCAGATAAAGCGCCTAAGTGCGGCAACCACTGCAACTTTGAAAGCGTGAATGAAGTCCATCAGAGATGTCGTTGACGGGGGGAAATACATAAGTTGTGCGGGATCACAAACTAAATTTTCGACCATACACGACAGAATGGGGGGAGTTGTCTACAACAGAAAGTGCGATATTATAGGACGGACGTCTAACTCGAGGGCTGATGCACAAGTTGATCCCTCATAATCTGTGTATTTCACCCCGCTGAGGCAAGTCTCAAACAGGTGGCACGACAACTACAATCGTCATTAACATGATGACAACGTAGGATGCACCGAAGTGGATAGCTTCGACTGACTTGAGACCCCTGCCAGCCCTTGCTGCATGGAATCGCAATACGCGTACGCAGGGTCTAGCAATTGCGAACATGGGTGTTCTACCCGAGAACTTAATAGGTCATGGCATACCAGAAGAGCAGCGGTAGATGGGTGCTAGGGCCCACAGTAGCCTAGGGTCGTGAGATTGGGCTTTCGCAGAATCGCGCTATTCCTCACAATTTGTGCATTAGCCCGCGGGTGATAACCACGTGCGGACTGGATTGCTTTGTACTCGTAGTCTCTATAGAGCGTGCTAGTAACAGTATGCTGAAGTGTTGGGAGAACTAAACTCACTATTCGTCTCTCGTTAGTCGACCATCCCTCACAATCTTTGTATATCGCCGCGCGCCGACTTAAGTTACTGTGTATGGATGGTGCCTACTCGCCAAGTCCAAACGCCGCCTGGTCACGCCTCTTAACTGCAAAGATTGCCCCGGGAGTGCGGCATAAGTGAGCTGACAGTGGAGCAAAAGATCCTGGTGGCCGAAGGCGCCAGCCCTCGTTTCGTAACAACAATAATCGTAAGGGCTTTGTAGGTTGAAGCTTAATCTCTGCCGGAGAGATGCCGGGGGCGTGGTGGCAGATCCTGTCAGCTAGTCCACGTAACCCGTGTTCTGACCGAGACATAAATGCTCATCCTTGCGTGGACTTATACACAAGTTGTGAGGAATCGCCACAGTTGGCCGCCGCGACTTCTACGTGTCGTTGAACCTAACTTGGCTGTACTAAGCCCGCTAACAACGCCCCTCGAGCCGAATATCTTTGGAGATTTACGTAATGTCTGCTCGGGGCGGTCTAATGCTATGATGCGACACGTAGTTCCGAGCAATTAAAAACCATTGTACAAGCGGCTCTCACCTGTCGGCTCAAACTTCTTAATTCCTCCGCTTTGCACTTGTCGGCGGCAATGGCAGCCTGTAGGATGGCCGTATAGAGATCAGAAAAGGGGAATCTTATCGTGGCTTCTTCTGACAGAGGGAAGCCAAGTTACAAATTCCCTTTCTGTGAGCGAGATCGTTAACATTCAAAGACGATAGATGTCCACCCCTACTGGAATCTATTCTATTTGGTGGCGAGCGAATCCACGGACTGCGAATCGCCTTCGTACCGTATGAAGCTAGCTGTCACGTTTTATTAGCCCCCCGCATAGTGTTATATTTACGCCATTTTGGGTTGGATCTGGGCGCGGACCATAGTGGTGAAAAGAAATCGCTCTCAATCAATCAAATTGGGGTTCAGAAAGCATACGATGTGCTCAGCGCTCACTGGACAATATTCCGGTGAGAAGCGAAGATGATGTGCACGCCCCCTCCTGCGCAAGCCACGGTGACCCATCCCACATACCAGCACCTTAACGATCTCGTTTTTAACGACACGTGGCGTCATGTACCAATATGCCGTACGGCCTCGTCTAGCTCCACGTCTATACCTCACAACCTCTGCATCAGCCCGCGGGTGCGCCCCTGTGGCAGACGTTGAGGGATTCTTTTCAAATACAAACGAAGTCAATTAACTTTTGAGCCTTTGACACCTAACATTCTTCTTGTTGGGCCTAGGGACTAGTCGTAAGGGTTCACTCCCTATAACTAAGACGCATCATAGGACTCCTCTCCGGTTCTCTGTGCGGTCGAGCTCCCTGCTTCTCCGGAAATAGAGTGTAATTCTATGCCAGACAGTAGAGCTTCCGATGGGTGGCTCTCACCGCGCGGGATCCTAATGTTAGCGGCCGTGAGCCTTCGGAACCGGAGTGCATTCCACACAATCTCTGCATAAGTCCCCGCTTTCGTCATTATGGCCACCGCATCCGGCACTGCTCCGATGTCTTAGTTGCAGAGTGTCTGATCCTCACGCCGTAACCACCGGAGGACGCCGCGCCGCGTGTCTATCGGCACGAGGTCCGTCAGCCGTTGAATGTCGTGTCAGGCGAATGGTGACCATACCTCGCGGTTGATACAAGACAGACTAGATGTTATCACAAGTTCCGGTTATGACCCGGCCTTGTGTGCCCAAAGTTACAGAGTCGGTGCAGTGTCTACCTGACCCCGAGGCCAGGATGCCTTATCGCGTTAGCCGCTCCAGGATAGCTCAGGCTGAGGGGAAACCTAAACGAAAATAAGCACTATACTGTCAGATCGGTCGCAGGAGTCCGACCGTACTTAGACGTCTTACTATAACTCGCTAATGTAACGCATCACCCGGTCATGGCAAAGGTTTCCTTTGTACCCCCTAGTAAATGCATAAACGAAGTTAACCGGTACAAAAATGTGTGCCCGTGCCATCTAACTCAGGACCAGCTTTTTTCACTAACAAAAGCAAACACGCGTTATCACAGGAAACAGGTATGTCGCGGTAAGTCCACGGTGCGGTGGTGGCTAAACGCCGAGTATCACAGCAGTATGCCTTTGTACTTCGCCGTGCCTCTTTTCAGCGTAGGTCTAGTTGACTGCATACATGGGACCCGTTCCCTGAGATTTCGTGGGATATCTACTCTCCGCAACGAAGACCATCTCTACAAAGGCACAGTTGTGCGCCGTCCTCAAAGTCCGTTGAGCCAACCGGTAATGTTTTGTTGTTCATAATCTGTATATGGCTATCATCAGGTCCAATCCTCAAGCATGTAGTATAAGAGTATGGCGGCTACGTTGACCGGATGATGGTGTGATTTCCAGCCGGACCGAGGGTGGCCCTTAATCCGATCGCTCATTAGTAATACCACCCTCACGTGCGAGAGTTTAGCTCTGCTAGATGGTCTCGCGCATATGACCAATTGTTGCGGTAACTACATATAGTCTACCCAGGCCTCGCCAGCGTTATAATAAAAGTCGTACATGCTAGAAAGTTCTTTCTTGGAAGCCACTAGATAGAGCCCGATTGGGTACGTGCATCTTTACGTTTGACCTTCGCCTTGGGACAACTCAAGGAGAGTGATATTTCGCGACAGATCGCACTGTCGCGGCCATCAGAGCGCCCGCGGCGTCTCCAATCGGTAAATTGTGTCCCCTTCTCTTCTTCACCCATATCACTTAAAACATGCACCACTGTCAGGTCTGGCCAGCCAAGTAAGACACCAAATTTGGCGGGTTGTGAAGATGGGTGTTGGACCGTGATCGGAGCGCCGATGACGCACCTCATTGCAAATTCACTGATACTCAATTTATGGCGTCCCAGCGAGTTACATAGGACGAGTCTCAGAAGGTTGGAAGGATTCTATCCTAGTTGGTGATCAGGGAGTGTTGTGACGGTACCTGCGCGGGCGGACTGGTTGCAACCAGTCCACGTACAGACATGTGCTATTGTTCTGCGCCCCCAAGTACACCCGTTTTATGTGGCAGGCTAAGCCTGATAGGTTATCGGTTCCTGCCGCCGCGGCGATATGCAAAGGTTGTGCGGAATCGAGATTCCTCCACAAAGATCACCACAAATATGAGAGCGCTTACGATCTCATAAACTTGCAGTGTTTTATGAGCCTGAGGGTAGGCGACGTATAATGCGATTTTTTAGGTTGTTTTAAGAAGTTCTGCTAGAGGGTGGGCTACAGCTAACGCTGCATAGACCCCCCTCCGGCCTTTAGTTCCTAAACCATTTGTTTTTGTTGATAATTCAGATCTCCATGTCGTGCTCACATATAGCAGTTGCGGCACACCATATCTTTCAATCGTACCTGTCGGAGCCTCTTCTCGACCCCCCTTGTGCCATTGACAATTAACGAGACTGGGTAACGATATTTTCCCCAACGACGTCTGCTTACACTTGTTATTATAACTGAAGAGCGATTGAAGTCGCCATACGACAGTCATTGGTTTTCTCGGCTAAGGCCGCAGGGAATAAATCTGGGCAGCGTAAAACTAGCACCAACTCCCATGAAGTATCCCAACGGGCGACAGAGTCACACGCTACGTCTAAAAGCCCGGGCCCTAATAACACAGTAAAAGGCCTCCCTGGAACTCGACGCTGAGTACACAGCTTAAACCATTAGTGATCTGTCTATCTTCGCTTCGCCTTGCAAGTCTCGTCAAGCGCAGAATCTTTACAAACTACGTGTCAACATACATGATCATTATAACCCTCTACCATACTCACTTTTGGGTTTACGTGACGTCGGGATTCTAATTACTTGCACGGCCGTTGCGCTCAGCAGTGCCCATCTTGTTTAATGGCGCACGCTCCGGCAACCTGGACCGGAGCTAATACCCCTGGGCTTAAAGGATATACGAATAAACTAAGCCACATAGCCATATTGACCCCCCCACTCCTAAGCGTTATAGCGGCGCCTTCATCGCGCTAACGAAAGTTTTAGGTCGAAGGGCATAGGGGGATAAATGCCATCGTAATCGATGAGAACGGCCGACAGCCAGTACCCCAAATTCGCGAAGCCGCGGCCTTCGCGAGCTTAAGAACGAATGCATTGTCGAGTGCACTCTAGGGGAGATGCAGAGATTGTGAGGTATGCGGTCTCCTGCAGGCAAATGTCTTGATTATGTGTGCAAGAGAGCAACTGGTAGAGTTGAATGTCCCCAGACCCAATTCCTCTTGTAACACTCCGAGCATGTTATCTATAAGCAAAGAGCGAAGCGGGGACGGAGCCGAACAGTCATGGTGAGTGATAAGGAGCTGGGGGAACAATTTAGTACCGCTTCGTGTTTAGCGAGGTAAACAGTCAGACCTTGTGTGTCATATACTTTGGATACCACATAACTTGTGCATCTCCCCACGGTATTTGTACCAAACCGTTTCACGGACCTATCTTCCCCTTGTCTCTCATTTCCCTTTATCCAGGAGGCAGGATATCTTTAAGGTCTCGT"
    # aa_arg = "IPHNLCISPR"

    # seq_out = peptide_encoding(dna_arg, aa_arg, genetic_code_dict)
    # [print(seq) for seq in seq_out]

    #######################################################################################################
    # Generate the theoretical spectrum of a cyclic peptide
    #######################################################################################################

    # peptide = "LHPHRGWEWIPFKYI"
    # l_out = cyclospectrum(peptide)
    # masses = ''
    # for el in l_out:
    #     masses = masses + str(el) + ' '
    # print(masses)

    # n_arg = 4
    # peptide = "NQEL"
    # l_out = cyclospectrum(peptide)
    # print(l_out)
    # print(len(l_out))

    #######################################################################################################
    # Cyclopeptide Sequencing Challenge - Given an experimental spectrum, find the peptide sequences that
    # are consistent with that spectrum
    #######################################################################################################

    # spectrum_arg = "0 71 97 101 113 131 131 131 156 168 198 202 232 244 262 269 269 287 299 329 333 375 388 400 400 400 400 430 446 485 501 531 531 531 531 543 556 598 602 632 644 662 662 669 687 699 729 733 763 775 800 800 800 818 830 834 860 931"
    # spectrum_arg = spectrum_arg.split(" ")
    # spectrum_arg = [int(s) for s in spectrum_arg]
    # spectrum_arg = sorted(spectrum_arg)
    # p_out = cyclopeptide_sequencing(spectrum_arg)
    # peptides_string = ''
    # for k in range(len(p_out)):
    #     peptide_string = ''
    #     current_peptide = p_out[k]
    #     for i in range(len(current_peptide)):
    #         peptide_string = peptide_string + str(current_peptide[i])
    #         if i < len(current_peptide) - 1:
    #             peptide_string = peptide_string + '-'
    #     peptides_string = peptides_string + peptide_string
    #     if k < len(p_out) - 1:
    #         peptides_string = peptides_string + ' '
    # print(peptides_string)

    #######################################################################################################
    # WEEK 4
    #######################################################################################################

    #######################################################################################################
    # Cyclopeptide Scoring Challenge - Compute the score of a cyclic peptide against a spectrum
    #######################################################################################################

    # p_arg = "QALYDFNCFDYAATTIHHRSCKTNFEFASTADDPISYFT"
    # s_arg = "0 71 71 71 71 87 87 97 97 101 101 101 101 101 103 103 113 113 113 113 113 114 114 115 115 128 128 129 129 131 137 137 142 147 147 147 147 147 147 156 158 163 172 184 188 190 199 202 210 210 214 215 217 218 226 228 229 229 231 232 234 242 243 243 244 244 246 248 250 250 259 260 261 261 262 273 274 276 276 276 278 293 300 305 305 312 315 318 319 323 323 325 331 332 343 344 346 347 347 349 351 355 357 358 359 362 364 364 365 376 377 380 386 387 389 390 390 391 406 406 413 419 420 423 425 428 430 434 434 436 438 446 447 452 452 456 457 459 460 461 474 479 483 488 490 490 491 494 496 502 504 505 505 505 507 511 517 521 523 528 533 535 537 537 543 547 551 556 557 560 560 565 567 569 575 575 576 581 589 593 594 599 599 603 608 608 608 611 618 618 618 619 620 620 622 626 630 638 642 644 644 652 657 660 666 668 670 670 680 680 682 682 689 689 689 689 695 700 703 707 709 712 712 713 715 721 722 723 731 731 731 733 735 745 748 755 757 757 765 766 769 771 781 783 784 789 792 796 796 802 802 804 804 809 811 813 813 816 817 826 827 828 828 832 834 836 836 837 849 860 861 862 868 869 870 870 872 872 882 885 886 887 894 894 897 899 903 912 915 918 920 924 927 928 928 931 933 935 939 940 941 941 949 951 956 957 958 962 962 963 965 973 974 975 983 985 986 991 999 1006 1009 1012 1012 1019 1021 1025 1027 1027 1028 1028 1032 1033 1041 1042 1046 1050 1054 1056 1059 1059 1062 1063 1063 1065 1067 1070 1076 1077 1080 1088 1099 1102 1104 1104 1110 1112 1114 1122 1125 1127 1128 1129 1133 1134 1137 1138 1143 1146 1146 1155 1156 1156 1159 1159 1160 1164 1165 1168 1168 1175 1177 1182 1183 1187 1193 1194 1205 1215 1217 1217 1217 1217 1217 1223 1230 1235 1236 1239 1240 1240 1246 1247 1249 1251 1252 1252 1256 1256 1258 1259 1259 1270 1271 1272 1276 1278 1283 1288 1288 1288 1297 1306 1306 1307 1312 1315 1318 1320 1324 1343 1345 1346 1346 1349 1352 1353 1355 1359 1359 1364 1364 1365 1368 1369 1371 1371 1373 1374 1375 1377 1383 1384 1386 1387 1389 1399 1403 1407 1410 1415 1416 1420 1420 1425 1430 1444 1446 1446 1446 1453 1457 1460 1461 1465 1466 1469 1472 1478 1480 1481 1483 1484 1487 1487 1490 1490 1493 1493 1496 1499 1502 1502 1502 1508 1516 1517 1517 1517 1518 1520 1523 1529 1544 1554 1557 1559 1559 1561 1567 1570 1574 1579 1583 1584 1584 1588 1588 1593 1593 1594 1594 1600 1603 1605 1607 1615 1616 1617 1617 1619 1620 1625 1627 1630 1630 1630 1631 1639 1640 1645 1657 1664 1664 1670 1671 1671 1672 1674 1676 1689 1689 1696 1697 1698 1701 1706 1707 1708 1716 1719 1722 1722 1730 1730 1731 1733 1733 1735 1740 1741 1745 1754 1754 1758 1758 1762 1763 1764 1767 1772 1772 1773 1776 1777 1777 1784 1785 1790 1793 1793 1802 1803 1810 1811 1827 1834 1835 1836 1836 1843 1843 1845 1845 1847 1848 1855 1855 1859 1859 1859 1860 1863 1864 1866 1867 1869 1874 1876 1877 1880 1882 1889 1891 1897 1903 1906 1906 1914 1920 1924 1930 1937 1939 1939 1940 1940 1948 1948 1948 1948 1950 1956 1960 1963 1964 1965 1974 1977 1977 1979 1983 1983 1990 1992 1994 1995 1995 2001 2004 2004 2006 2007 2019 2021 2021 2026 2031 2033 2040 2045 2050 2053 2053 2053 2054 2062 2064 2066 2076 2077 2077 2078 2080 2084 2084 2086 2090 2091 2091 2095 2095 2101 2102 2104 2105 2110 2110 2116 2117 2120 2121 2121 2126 2135 2150 2155 2158 2162 2164 2167 2167 2168 2173 2177 2179 2179 2182 2187 2187 2187 2190 2191 2191 2192 2192 2200 2206 2209 2209 2212 2214 2214 2218 2223 2224 2227 2229 2232 2233 2238 2238 2241 2265 2268 2268 2273 2274 2277 2279 2282 2283 2288 2292 2292 2294 2297 2297 2300 2306 2314 2314 2315 2315 2316 2319 2319 2319 2324 2327 2327 2329 2333 2338 2339 2339 2342 2344 2348 2351 2356 2371 2380 2385 2385 2386 2389 2390 2396 2396 2401 2402 2404 2405 2411 2411 2415 2415 2416 2420 2422 2422 2426 2428 2429 2429 2430 2440 2442 2444 2452 2453 2453 2453 2456 2461 2466 2473 2475 2480 2485 2485 2487 2499 2500 2502 2502 2505 2511 2511 2512 2514 2516 2523 2523 2527 2529 2529 2532 2541 2542 2543 2546 2550 2556 2558 2558 2558 2558 2566 2566 2567 2567 2569 2576 2582 2586 2592 2600 2600 2603 2609 2615 2617 2624 2626 2629 2630 2632 2637 2639 2640 2642 2643 2646 2647 2647 2647 2651 2651 2658 2659 2661 2661 2663 2663 2670 2670 2671 2672 2679 2695 2696 2703 2704 2713 2713 2716 2721 2722 2729 2729 2730 2733 2734 2734 2739 2742 2743 2744 2748 2748 2752 2752 2761 2765 2766 2771 2773 2773 2775 2776 2776 2784 2784 2787 2790 2798 2799 2800 2805 2808 2809 2810 2817 2817 2830 2832 2834 2835 2835 2836 2842 2842 2849 2861 2866 2867 2875 2876 2876 2876 2879 2881 2886 2887 2889 2889 2890 2891 2899 2901 2903 2906 2912 2912 2913 2913 2918 2918 2922 2922 2923 2927 2932 2936 2939 2945 2947 2947 2949 2952 2962 2977 2983 2986 2988 2989 2989 2989 2990 2998 3004 3004 3004 3007 3010 3013 3013 3016 3016 3019 3019 3022 3023 3025 3026 3028 3034 3037 3040 3041 3045 3046 3049 3053 3060 3060 3060 3062 3076 3081 3086 3086 3090 3091 3096 3099 3103 3107 3117 3119 3120 3122 3123 3129 3131 3132 3133 3135 3135 3137 3138 3141 3142 3142 3147 3147 3151 3153 3154 3157 3160 3160 3161 3163 3182 3186 3188 3191 3194 3199 3200 3200 3209 3218 3218 3218 3223 3228 3230 3234 3235 3236 3247 3247 3248 3250 3250 3254 3254 3255 3257 3259 3260 3266 3266 3267 3270 3271 3276 3283 3289 3289 3289 3289 3289 3291 3301 3312 3313 3319 3323 3324 3329 3331 3338 3338 3341 3342 3346 3347 3347 3350 3350 3351 3360 3360 3363 3368 3369 3372 3373 3377 3378 3379 3381 3384 3392 3394 3396 3402 3402 3404 3407 3418 3426 3429 3430 3436 3439 3441 3443 3443 3444 3447 3447 3450 3452 3456 3460 3464 3465 3473 3474 3478 3478 3479 3479 3481 3485 3487 3494 3494 3497 3500 3507 3515 3520 3521 3523 3531 3532 3533 3541 3543 3544 3544 3548 3549 3550 3555 3557 3565 3565 3566 3567 3571 3573 3575 3578 3578 3579 3582 3586 3588 3591 3594 3603 3607 3609 3612 3612 3619 3620 3621 3624 3634 3634 3636 3636 3637 3638 3644 3645 3646 3657 3669 3670 3670 3672 3674 3678 3678 3679 3680 3689 3690 3693 3693 3695 3697 3702 3702 3704 3704 3710 3710 3714 3717 3722 3723 3725 3735 3737 3740 3741 3749 3749 3751 3758 3761 3771 3773 3775 3775 3775 3783 3784 3785 3791 3793 3794 3794 3797 3799 3803 3806 3811 3817 3817 3817 3817 3824 3824 3826 3826 3836 3836 3838 3840 3846 3849 3854 3862 3862 3864 3868 3876 3880 3884 3886 3886 3887 3888 3888 3888 3895 3898 3898 3898 3903 3907 3907 3912 3913 3917 3925 3930 3931 3931 3937 3939 3941 3946 3946 3949 3950 3955 3959 3963 3969 3969 3971 3973 3978 3983 3985 3989 3995 3999 4001 4001 4001 4002 4004 4010 4012 4015 4016 4016 4018 4023 4027 4032 4045 4046 4047 4049 4050 4054 4054 4059 4060 4068 4070 4072 4072 4076 4078 4081 4083 4086 4087 4093 4100 4100 4115 4116 4116 4117 4119 4120 4126 4129 4130 4141 4142 4142 4144 4147 4148 4149 4151 4155 4157 4159 4159 4160 4162 4163 4174 4175 4181 4183 4183 4187 4188 4191 4194 4201 4201 4206 4213 4228 4230 4230 4230 4232 4233 4244 4245 4245 4246 4247 4256 4256 4258 4260 4262 4262 4263 4263 4264 4272 4274 4275 4277 4277 4278 4280 4288 4289 4291 4292 4296 4296 4304 4307 4316 4318 4322 4334 4343 4348 4350 4359 4359 4359 4359 4359 4359 4364 4369 4369 4375 4377 4377 4378 4378 4391 4391 4392 4392 4393 4393 4393 4393 4393 4403 4403 4405 4405 4405 4405 4405 4409 4409 4419 4419 4435 4435 4435 4435 4506"
    # s_arg = s_arg.strip().split(' ')
    # s_arg = [int(m) for m in s_arg]
    # out = cyclopeptide_scoring(p_arg, s_arg)
    # print(out)

    #######################################################################################################
    # Leaderboard Cyclopeptide Sequencing Challenge - Given an integer and experimental spectrum, find a
    # peptide sequence that is consistent with that spectrum
    #######################################################################################################

    # n_arg = 1000
    # s_arg = "0 97 99 114 128 147 147 163 186 227 241 242 244 260 261 262 283 291 333 340 357 385 389 390 390 405 430 430 447 485 487 503 504 518 543 544 552 575 577 584 632 650 651 671 672 690 691 738 745 747 770 778 779 804 818 819 820 835 837 875 892 917 932 932 933 934 965 982 989 1030 1039 1060 1061 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1225 1322"
    # s_arg = s_arg.strip().split()
    # s_arg = [int(n) for n in s_arg]
    # out = leaderboard_cyclopeptide_sequencing(s_arg, n_arg)
    # c_out = ''
    # for i in range(len(out)):
    #     c_out = c_out + str(out[i])
    #     if i < len(out) - 1:
    #         c_out = c_out + '-'
    # print(c_out)

    #######################################################################################################
    # Spectral Convolution Problem - Given an experimental spectrum, compute the list of elements in the
    # spectrum's convolution
    #######################################################################################################

    # s_arg = "0 57 97 113 128 128 137 147 154 163 185 210 265 267 275 276 282 284 291 332 373 395 404 412 412 428 429 430 469 501 540 541 542 558 558 566 575 597 638 679 686 688 694 695 703 705 760 785 807 816 823 833 842 842 857 873 913 970"
    # s_arg = s_arg.strip().split()
    # s_arg = [int(m) for m in s_arg]

    # out = spectral_convolution(s_arg)
    # c_out = ''
    # for i in range(len(out)):
    #     c_out += str(out[i])
    #     if i < len(out) - 1:
    #         c_out += ' '

    # print(c_out)

    #######################################################################################################
    # Convolution Cyclopeptide Sequencing - Given an integer m, another integer, n, and an experimental
    # spectrum, generate the cyclic peptide whose composite amino acids are most inline with the
    # given spectrum. m indicates the number of top elements to choose from the spectral convolution and
    # n indicates the number of most consistent peptides to retain during the bound step
    #######################################################################################################

    m_arg = 20
    n_arg = 1000
    s_arg = "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"
    s_arg = s_arg.strip().split()
    s_arg = [int(m) for m in s_arg]

    out = convolution_cyclopeptide_sequencing(m_arg, n_arg, s_arg)
    s_out = ''
    for i in range(len(out)):
        s_out += str(out[i])
        if i < len(out) - 1:
            s_out += '-'

    print(s_out)




