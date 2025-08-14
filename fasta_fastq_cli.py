import argparse
from parser import SequenceParser
from stats import count_sequences, get_sequence_lengths, calculate_basic_stats, calculate_length_distributions, calculate_n50_l50, calculate_n90_l90, calculate_assembly_stats, calculate_gc_content, get_nucleotide_composition, calculate_gc_stats, count_ambiguous_bases
from quality import detect_phred_encoding, calculate_average_quality, analyze_quality_by_position, get_average_quality_per_sequence, get_quality_distribution, find_low_quality_positions, get_quality_stats
from filtering import filter_by_length, get_length_filter_preview, filter_by_quality, filter_by_min_quality, trim_records_by_quality, filter_by_gc_content, filter_ambiguous_sequences, remove_duplicate_sequences
from plotting import plot_length_histogram, plot_length_distribution, plot_cumulative_length, plot_quality_histogram, plot_quality_by_position, plot_gc_distribution, plot_nucleotide_composition

def load_sequences(filename):
    parser = SequenceParser(filename)
    file_format = parser.detect_format()
    if file_format == "Invalid file":
        raise ValueError("Invalid input file format.")

    with parser.open_compressed_file(filename, "rt") as fh:
        if file_format == "FASTA":
            records = list(parser.parse_fasta(fh))
        elif file_format == "FASTQ":
            records = list(parser.parse_fastq(fh))
        else:
            raise ValueError("Unknown file format.")

    return records

def main():
    main_parser = argparse.ArgumentParser(description="FASTA/FASTQ Command Line Interface")
    subparsers = main_parser.add_subparsers(dest="command", required=True)

    # STATS
    stats_count_sequences = subparsers.add_parser("count_sequences", help="Count number of sequences.")
    stats_count_sequences.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")

    stats_get_sequence_lengths = subparsers.add_parser("get_sequence_lengths", help="Counts the length of each sequence.")
    stats_get_sequence_lengths.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")

    stats_calculate_basic_stats = subparsers.add_parser("calculate_basic_stats", help="Calculate basic statistics on sequence lengths.")
    stats_calculate_basic_stats.add_argument("--input", required=True, help="Input FASTA or FASTA file.")

    stats_calculate_length_distributions = subparsers.add_parser("calculate_length_distributions", help="Groups sequence lengths into bins of size bin_size.")
    stats_calculate_length_distributions.add_argument("--input", required=True, help="Input FASTA or FASTA file.")
    stats_calculate_length_distributions.add_argument("--bin_size", type=int, required=True)

    stats_calculate_n50_l50 = subparsers.add_parser("calculate_n50_l50", help="Calculates the N50 and L50 values.")
    stats_calculate_n50_l50.add_argument("--input", required=True, help="Input FASTA or FASTA file.")

    stats_calculate_n90_l90 = subparsers.add_parser("calculate_n90_l90", help="Calculates the N90 and L90 values.")
    stats_calculate_n90_l90.add_argument("--input", required=True, help="Input FASTA or FASTA file.")

    stats_calculate_assembly_stats = subparsers.add_parser("calculate_assembly_stats", help="Calculates a summary of statistics.")
    stats_calculate_assembly_stats.add_argument("--input", required=True, help="Input FASTA or FASTA file.")

    stats_calculate_gc_content = subparsers.add_parser("calculate_gc_content", help="Calculates the GC% per sequence.")
    stats_calculate_gc_content.add_argument("--input", required=True, help="Input FASTA or FASTA file.")

    stats_get_nucleotide_composition = subparsers.add_parser("get_nucleotide_composition", help="Counts each nucleotide per sequence.")
    stats_get_nucleotide_composition.add_argument("--input", required=True, help="Input FASTA or FASTA file.")

    stats_calculate_gc_stats = subparsers.add_parser("calculate_gc_stats", help="Calculates min, max, median, and mean of GC% values.")
    stats_calculate_gc_stats.add_argument("--input", required=True, help="Input FASTA or FASTA file.")

    stats_count_ambiguous_bases = subparsers.add_parser("count_ambiguous_bases", help="Counts the number of ambiguous (N and other) and unambigous bases.")
    stats_count_ambiguous_bases.add_argument("--input", required=True, help="Input FASTA or FASTA file.")


    # FILTERING
    filter_length_parser = subparsers.add_parser("filter_by_length", help="Filter sequences by length")
    filter_length_parser.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    filter_length_parser.add_argument("--min_length", type=int, required=True)
    filter_length_parser.add_argument("--max_length", type=int, required=True)

    filter_length_preview_parser = subparsers.add_parser("get_length_filter_preview", help="Returns length statistics.")
    filter_length_preview_parser.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    filter_length_preview_parser.add_argument("--min_length", type=int, required=True)
    filter_length_preview_parser.add_argument("--max_length", type=int, required=True)

    filter_avg_quality_parser = subparsers.add_parser("filter_avg_quality", help="Filters sequences by minimum average quality.")
    filter_avg_quality_parser.add_argument("--input", required=True, help="Input FASTQ file.")
    filter_avg_quality_parser.add_argument("--min_avg_quality", type=int, required=True)
    filter_avg_quality_parser.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")

    filter_min_quality_parser = subparsers.add_parser("filter_min_quality", help="Filters sequences by minimum base quality.")
    filter_min_quality_parser.add_argument("--input", required=True, help="Input FASTQ file.")
    filter_min_quality_parser.add_argument("--min_quality", type=int, required=True)
    filter_min_quality_parser.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")

    trim_records_quality_parser = subparsers.add_parser("trim_records_by_quality", help="Trims low quality bases from both ends of all sequences.")
    trim_records_quality_parser.add_argument("--input", required=True, help="Input FASTQ file.")
    trim_records_quality_parser.add_argument("--min_quality", type=int, required=True)
    trim_records_quality_parser.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")

    filter_gc_content_parser = subparsers.add_parser("filter_gc_content", help="Filters sequences based on minimum and maximum GC %.")
    filter_gc_content_parser.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    filter_gc_content_parser.add_argument("--min_gc", type=int, required=True)
    filter_gc_content_parser.add_argument("--max_gc", type=int, required=True)

    filter_ambiguous_sequences_parser = subparsers.add_parser("filter_ambiguous_sequences", help="Filters sequences based on ambiguous base (N) %.")
    filter_ambiguous_sequences_parser.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    filter_ambiguous_sequences_parser.add_argument("--max_n_percent", type=int, required=True)

    filter_duplicate_sequences_parser = subparsers.add_parser("remove_duplicate_sequences", help="Eliminates identical sequences to retain only unique entries.")
    filter_duplicate_sequences_parser.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")


    # QUALITY
    quality_detect_phred_encoding = subparsers.add_parser("detect_phred_encoding", help="Detects encoding based on quality sequence.")
    quality_detect_phred_encoding.add_argument("--input", required=True, help="Input FASTQ file.")

    quality_calculate_average_quality = subparsers.add_parser("calculate_average_quality_all", help="Calculates the average quality of all sequences ")
    quality_calculate_average_quality.add_argument("--input", required=True, help="Input FASTQ file.")
    quality_calculate_average_quality.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")

    quality_analyze_quality_by_position = subparsers.add_parser("calculate_quality_position", help="Calculates the average quality score for each position across all sequences.")
    quality_analyze_quality_by_position.add_argument("--input", required=True, help="Input FASTQ file.")
    quality_analyze_quality_by_position.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")  

    quality_get_average_quality_per_sequence = subparsers.add_parser("calculate_average_quality_seq", help="Calculates the average quality score per sequence.")
    quality_get_average_quality_per_sequence.add_argument("--input", required=True, help="Input FASTQ file.")
    quality_get_average_quality_per_sequence.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")

    quality_get_quality_distribution = subparsers.add_parser("get_quality_distribution", help="Counts the occurrences of each quality score.")
    quality_get_quality_distribution .add_argument("--input", required=True, help="Input FASTQ file.")
    quality_get_quality_distribution .add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")

    quality_find_low_quality_positions = subparsers.add_parser("find_low_quality_positions", help="Finds the positions of each read below given threshold.")
    quality_find_low_quality_positions.add_argument("--input", required=True, help="Input FASTQ file.")
    quality_find_low_quality_positions.add_argument("--threshold", type=int, required=True, help="Input minimum quality score")
    quality_find_low_quality_positions.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")

    quality_get_quality_stats = subparsers.add_parser("get_quality_stats", help="Returns a dictionary summarizing quality stats.")
    quality_get_quality_stats.add_argument("--input", required=True, help="Input FASTQ file.")
    quality_get_quality_stats.add_argument("--threshold", type=int, required=True, help="Input minimum quality score")
    quality_get_quality_stats.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")


    # PLOTTING
    plotting_plot_length_histogram = subparsers.add_parser("plot_length_histogram", help="Generates and saves a histogram of sequence lengths.")
    plotting_plot_length_histogram.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    plotting_plot_length_histogram.add_argument("--output_path", type=str, help="Input output path.")

    plotting_plot_length_distribution = subparsers.add_parser("plot_length_distribution", help="Generates and saves a binned bar graph of sequence lengths.")
    plotting_plot_length_distribution.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    plotting_plot_length_distribution.add_argument("--bin_size", required=True, type=int, help="Input bin size.")
    plotting_plot_length_distribution.add_argument("--output_path", type=str, help="Input output path.")

    plotting_plot_cumulative_length = subparsers.add_parser("plot_cumulative_length", help="Generates and saves a line plot of cumulative sequence length.")
    plotting_plot_cumulative_length.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    plotting_plot_cumulative_length.add_argument("--output_path", type=str, help="Input output path.")

    plotting_plot_quality_histogram = subparsers.add_parser("plot_quality_histogram", help="Generates and saves a histogram of sequence quality scores.")
    plotting_plot_quality_histogram.add_argument("--input", required=True, help="Input FASTQ file.")
    plotting_plot_quality_histogram.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")
    plotting_plot_quality_histogram.add_argument("--output_path", type=str, help="Input output path.")

    plotting_plot_quality_by_position = subparsers.add_parser("plot_quality_by_position", help="Generates and saves a line plot of average quality score per position.")
    plotting_plot_quality_by_position.add_argument("--input", required=True, help="Input FASTQ file.")
    plotting_plot_quality_by_position.add_argument("--phred_offset", type=int, required=False, help="Input phred encoding, if known (33 or 64).")
    plotting_plot_quality_by_position.add_argument("--output_path", type=str, help="Input output path.")

    plotting_plot_gc_distribution = subparsers.add_parser("plot_gc_distribution", help="Generates and saves a histogram of GC content percentage.")
    plotting_plot_gc_distribution.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    plotting_plot_gc_distribution.add_argument("--bin_size", type=int, required=False, help="Specify the desired number of bins.")
    plotting_plot_gc_distribution.add_argument("--output_path", type=str, help="Input output path.")

    plotting_plot_nucleotide_composition = subparsers.add_parser("plot_nucleotide_composition", help="Generates and saves a bar graph of total counts for each nucleotide.")
    plotting_plot_nucleotide_composition.add_argument("--input", required=True, help="Input FASTA or FASTQ file.")
    plotting_plot_nucleotide_composition.add_argument("--output_path", type=str, help="Input output path.")

    args = main_parser.parse_args()

    try:
        records = load_sequences(args.input)

        # STATS
        if args.command == "count_sequences":
            result = count_sequences(records)
            print(result)
        elif args.command == "get_sequence_lengths":
            result = get_sequence_lengths(records)
            print(result)
        elif args.command == "calculate_basic_stats":
            lengths = get_sequence_lengths(records)
            result = calculate_basic_stats(lengths)
            print(result)
        elif args.command == "calculate_length_distributions":
            lengths = get_sequence_lengths(records)
            result = calculate_length_distributions(lengths, args.bin_size)
            print(result)
        elif args.command == "calculate_n50_l50":
            lengths = get_sequence_lengths(records)
            result = calculate_n50_l50(lengths)
            print(result)
        elif args.command == "calculate_n90_l90":
            lengths = get_sequence_lengths(records)
            result = calculate_n90_l90(lengths)
            print(result)
        elif args.command == "calculate_assembly_stats":
            lengths = get_sequence_lengths(records)
            result = calculate_assembly_stats(lengths)
            print(result)
        elif args.command == "calculate_gc_content":
            result = calculate_gc_content(records)
            print(result)
        elif args.command == "get_nucleotide_composition":
            result = get_nucleotide_composition(records)
            print(result)
        elif args.command == "calculate_gc_stats":
            gc_content = calculate_gc_content(records)
            result = calculate_gc_stats(gc_content)
            print(result)
        elif args.command == "count_ambiguous_bases":
            results = count_ambiguous_bases(records)
            print(results)
        # FILTERING
        elif args.command == "filter_by_length":
            results = filter_by_length(records, args.min_length, args.max_length)
            print(results)
        elif args.command == "get_length_filter_preview":
            results = get_length_filter_preview(records, args.min_length, args.max_length)
            print(results)
        elif args.command == "filter_avg_quality":
            results = filter_by_quality(records, args.min_avg_quality, args.phred_offset)
            for record in results:
                print(f"@{record.id}")
                print(record.sequence)
                print("+")
                print(record.quality)
        elif args.command == "filter_min_quality":
            results = filter_by_min_quality(records, args.min_quality, args.phred_offset)
            for record in results:
                print(f"@{record.id}")
                print(record.sequence)
                print("+")
                print(record.quality)
        elif args.command == "trim_records_by_quality":
            results = trim_records_by_quality(records, args.min_quality, args.phred_offset)
            for record in results:
                print(f"@{record.id}")
                print(record.sequence)
                print("+")
                print(record.quality)
        elif args.command == "filter_gc_content":
            results = filter_by_gc_content(records, args.min_gc, args.max_gc)
            for record in results:
                print(f"@{record.id}")
                print(record.sequence)
                if record.quality:
                    print("+")
                    print(record.quality)
        elif args.command == "filter_ambiguous_sequences":
            results = filter_ambiguous_sequences(records,args.max_n_percent)
            for record in results:
                print(f"@{record.id}")
                print(record.sequence)
                if record.quality:
                    print("+")
                    print(record.quality)
        elif args.command == "remove_duplicate_sequences":
            results = remove_duplicate_sequences(records)
            for record in results:
                print(f"@{record.id}")
                print(record.sequence)
                if record.quality:
                    print("+")
                    print(record.quality)
        # QUALITY
        elif args.command == "detect_phred_encoding":
            results = detect_phred_encoding(records)
            print(results)
        elif args.command == "calculate_average_quality_all":
            results = calculate_average_quality(records, args.phred_offset)
            print(results)
        elif args.command == "calculate_quality_position":
            results = analyze_quality_by_position(records, args.phred_offset)
            print(results)
        elif args.command == "calculate_average_quality_seq":
            results = get_average_quality_per_sequence(records, args.phred_offset)
            print(results)
        elif args.command == "get_quality_distribution":
            results = get_quality_distribution(records, args.phred_offset)
            print(results)
        elif args.command == "find_low_quality_positions":
            results = find_low_quality_positions(records, args.threshold, args.phred_offset)
            print(results)
        elif args.command == "get_quality_stats":
            results = get_quality_stats(records, args.phred_offset, args.threshold)
            print(results)
        # PLOTTING
        elif args.command == "plot_length_histogram":
            lengths = get_sequence_lengths(records)
            plot_length_histogram(lengths, args.output_path)
            print("Plotting sequence length histogram...")
        elif args.command == "plot_length_distribution":
            lengths = get_sequence_lengths(records)
            plot_length_distribution(lengths, args.output_path, args.bin_size)
            print("Plotting sequence length distribution...")
        elif args.command == "plot_cumulative_length":
            lengths = get_sequence_lengths(records)
            plot_cumulative_length(lengths, args.output_path)
            print("Plotting cumulative sequence length...")
        elif args.command == "plot_quality_histogram":
            plot_quality_histogram(records, args.phred_offset, args.output_path)
            print("Plotting sequence quality score histogram...")
        elif args.command == "plot_quality_by_position":
            plot_quality_by_position(records, args.phred_offset, args.output_path)
            print("Plotting average quality score by sequence position...")
        elif args.command ==  "plot_gc_distribution":
            plot_gc_distribution(records, args.bin_size, args.output_path)
            print("Printing histogram of GC content percentage...")
        elif args.command == "plot_nucleotide_composition":
            plot_nucleotide_composition(records, args.output_path)
            print("Plotting nucleotide count distribution...")
    
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())



