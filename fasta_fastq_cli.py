import argparse
from parser import SequenceParser
from stats import *
from quality import *
from filtering import *
from plotting import *

def main():
    main_parser = argparse.ArgumentParser("FASTA/FASTQ Command Line Interface")
    subparsers = main_parser.add_subparsers(dest="command", required=True)

    parse_parser = subparsers.add_parser("parse", help="Parse input file and detect format.")
    parse_parser.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")

    stats_parser = subparsers.add_parser("stats", help="Calculate sequence statistics (GC content, nucleotide composition, length summary).")
    stats_parser.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    stats_group = stats_parser.add_mutually_exclusive_group(required=True)
    stats_group.add_argument("--count_sequences", action="store_true", help="Counts number of sequences in input file.")
    stats_group.add_argument("--get_sequence_lengths", action="store_true", help="Counts the length of each sequence in input file.")
    stats_group.add_argument("--calculate_basic_stats", action="store_true", help="Calculate basic length statistics (count, min, max, median, mean).")
    stats_group.add_argument("--calculate_n50_l50", action="store_true", help="Calculates N50 and L50 on input file.")
    stats_group.add_argument("--calculate_n90_l90", action="store_true", help="Calculates N90 and L90 on input file.")
    stats_group.add_argument("--calculate_assembly_stats", action="store_true", help="Calculates more extensive length statistics.")
    stats_group.add_argument("--calculate_gc_content", action="store_true", help="Calculates GC content.")
    stats_group.add_argument("--get_nucleotide_composition", action="store_true", help="Calculates nucleotide composition.")
    stats_group.add_argument("--calculate_gc_stats", action="store_true", help="Calculates GC statistics.")
    stats_group.add_argument("--count_ambiguous_bases", action="store_true", help="Calculates ambiguous base statistics.")



    quality_parser = subparsers.add_parser(
        "quality", 
        help="Analyze FASTQ quality scores and detect phred encoding."
    )
    quality_parser.add_argument(
        "--input", 
        required=True, 
        help="Input FASTQ file."
    )

    filtering_parser = subparsers.add_parser(
        "filter", 
        help="Filter sequences by length, quality, GC content, ambiguous bases, or duplicates."
    )
    filtering_parser.add_argument(
        "--input", 
        required=True, 
        help="Input FASTA or FASTQ file."
    )

    plotting_parser = subparsers.add_parser(
        "plot", 
        help="Generate plots for sequence quality and statistics."
    )
    plotting_parser.add_argument(
        "--input", 
        required=True, 
        help="Input FASTA or FASTQ file."
    )
    plotting_parser.add_argument(
        "--output", 
        required=True, 
        help="Output path for plot image."
    )

    args = main_parser.parse_args()
