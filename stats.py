import statistics
import re

def get_sequence_lengths(records):
    return [len(seq.sequence) for seq in records]

def calculate_basic_stats(lengths):
    if not lengths:
        return None
    stats = {
        "count": len(lengths),
        "min": min(lengths),
        "max": max(lengths),
        "median": statistics.median(lengths),
        "mean": sum(lengths)/len(lengths)
    }
    return stats

def calculate_length_distributions(lengths, bin_size):
    if not lengths:
        return None
    distribution = {}
    for length in lengths:
        bin = (length // bin_size) * bin_size
        distribution[bin] = distribution.get(bin, 0) + 1
    return dict(sorted(distribution.items()))
    
def calculate_n50_l50(lengths):
    if not lengths:
        return None
    sorted_lengths = sorted(lengths, reverse=True)
    half_sum= sum(lengths)/2
    running_total = 0
    number_contigs = 0
    for length in sorted_lengths:
        running_total += length
        number_contigs += 1
        if running_total >= half_sum:
            fifty_stats = {
                "n50" : length,
                "l50" : number_contigs
            }
            return fifty_stats

def calculate_n90_l90(lengths):
    if not lengths:
        return None
    sorted_lengths = sorted(lengths, reverse=True)
    ninety= sum(lengths)*.90
    running_total = 0
    number_contigs = 0
    for length in sorted_lengths:
        running_total += length
        number_contigs += 1
        if running_total >= ninety:
            ninety_stats = {
                "n50" : length,
                "l50" : number_contigs
            }
            return ninety_stats
    
def calculate_assembly_stats(lengths):
    if not lengths:
        return None
    
    number_bases = sum(lengths)
    large_contigs = len([l for l in lengths if l >= 10000])  # >= 10kb
    medium_contigs = len([l for l in lengths if 1000 <= l < 10000])  # 1-10kb
    small_contigs = len([l for l in lengths if l < 1000])  # < 1kb
    basic_stats = calculate_basic_stats(lengths)
    n50_l50 = calculate_n50_l50(lengths)
    n90_l90 = calculate_n90_l90 (lengths)
    assembly_stats = {
        "number_bases" : number_bases,
        **basic_stats,
        **n50_l50,
        **n90_l90,
        "large_contigs_10kb_plus" : large_contigs,
        "medium_contigs_1_10kb" : medium_contigs,
        "small_contigs_under_1kb": small_contigs,
    }
    return assembly_stats

def calculate_gc_content(records):
    GC_content_all = []
    for seq in records:
        sequence = seq.sequence.upper()
        G = sequence.count("G")
        C = sequence.count("C")
        seq_length = len(sequence)
        GC_content = round((G + C) / seq_length * 100,2)
        GC_content_all.append(GC_content)
    return GC_content_all

def get_nucleotide_composition(records):
    composition = []
    for seq in records:
        sequence = seq.sequence.upper()
        G = sequence.count("G")
        C = sequence.count("C")
        A = sequence.count("A")
        T = sequence.count("T")
        nucleotides = {
            "A" : A,
            "T" : T,
            "G" : G,
            "C" : C
        }
        composition.append(nucleotides)
    return composition

def calculate_gc_stats(gc_contents):
    if not gc_contents:
        return None
    return {
        "min_gc" : min(gc_contents),
        "max_gc" : max(gc_contents),
        "median_gc" : round(statistics.median(gc_contents),2),
        "mean_gc" : round(statistics.mean(gc_contents), 2)
    }

def count_ambiguous_bases(records):
    unambiguous_all = []
    for seq in records:
        sequence = seq.sequence.upper()
        N = sequence.count("N")
        other = len(re.findall(r'[^ATGCN]', sequence))
        unambiguous = len(sequence) - N - other
        unambiguous_all.append({
            "N" : N,
            "other_ambiguous_bases" : other,
            "total_ambiguous" : (N+other),
            "unambiguous" : unambiguous,
            "percent_ambiguous" : format((N+other)/(len(sequence))*100,".1f")
            })
    return unambiguous_all