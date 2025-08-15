# fasta-fastq-parser
A command-line tool for comprehensive DNA sequence analysis from FASTA and FASTQ files, including statistics, quality assessment, filtering, and plotting. Supports compressed files (.gz, .bz2) and both single- and multi-line sequences.

## Features
- Parsing: Efficiently read FASTA and FASTQ files, including compressed `.gz` and `.bz2` files.
- Statistics: Count sequences, calculate lengths, GC content, N50/L50, N90/L90, nucleotide composition, and ambiguous bases.
- Quality Analysis: Detect Phred encoding, compute average quality per sequence and position, identify low-quality positions, and generate quality distributions.
- Filtering: Filter sequences by length, GC content, quality, ambiguous bases, or remove duplicate sequences.
- Plotting: Generate sequence length histograms, binned distributions, cumulative length plots, and quality score histograms.
- Flexible CLI: Run analyses and filtering commands with intuitive arguments and options.

## Folder Structure
```
fasta-fastq-parser/
│
├─ parser.py          # FASTA/FASTQ parsing
├─ stats.py           # Sequence statistics functions
├─ quality.py         # Quality analysis functions
├─ filtering.py       # Sequence filtering functions
├─ plotting.py        # Plotting functions
├─ fasta_fastq_cli.py # Main CLI entry point
├─ example_files/     # Sample FASTA/FASTQ files
└─ example_plots/     # Example output plots
```

## Dependencies
### Python
- Python 3.7+

### Standard Library (no installation needed)
- `os`, `gzip`, `bz2`, `statistics`, `argparse`, `collections (defaultdict)`, `typing (List, Optional, Dict, Union, Tuple)`

### Third-Party Libraries
- `numpy`
- `matplotlib`

Install third-party dependencies via pip:
```
pip install numpy matplotlib
```

## Quick Start
1. Clone the repository
```
git clone https://github.com/karo-pajdak/fasta-fastq-parser.git
cd fasta-fastq-parser
```
2. Install dependencies (see above)
3. Run the CLI with `--help` to see available commands:
```
python3 fasta_fastq_cli.py --help
```

## Example Commands
- Count sequences in a FASTA file:
```
python3 fasta_fastq_cli.py count_sequences --input example_files/sample.fasta
```

- Get sequence lengths:
```
python3 fasta_fastq_cli.py get_sequence_lengths --input example_files/sample.fastq
```

- Filter sequences by length:
```
python3 fasta_fastq_cli.py filter_by_length --input example_files/sample.fasta --min_length 100 --max_length 1000
```

- Plot sequence length histogram:
```
python3 fasta_fastq_cli.py plot_length_histogram --input example_files/sample.fasta --output example_plots/length_hist.png
```

- Filter sequences by minimum average quality (FASTQ):
```
python3 fasta_fastq_cli.py filter_avg_quality --input example_files/sample.fastq --min_avg_quality 30
```

## Supported File Types
- FASTA: `.fasta`, `.fa`
- FASTQ: `.fastq`. `.fq`
- Compressed: `.gz`, `.bz2`

## Arguments Summary
(Use `python3 fasta_fastq_cli.py --help` for the full, up-to-date list of commands.)

| Command                          | Required Arguments                        | Optional Arguments           | Description                                     |
| -------------------------------- | ----------------------------------------- | ---------------------------- | ----------------------------------------------- |
| `count_sequences`                | `--input`                                 | -                            | Count total sequences                           |
| `get_sequence_lengths`           | `--input`                                 | -                            | Print lengths of each sequence                  |
| `calculate_basic_stats`          | `--input`                                 | -                            | Compute min, max, mean, median sequence lengths |
| `calculate_length_distributions` | `--input`, `--bin_size`                   | -                            | Bin sequence lengths                            |
| `calculate_n50_l50`              | `--input`                                 | -                            | Calculate N50 and L50                           |
| `calculate_n90_l90`              | `--input`                                 | -                            | Calculate N90 and L90                           |
| `calculate_gc_content`           | `--input`                                 | -                            | GC% per sequence                                |
| `get_nucleotide_composition`     | `--input`                                 | -                            | Count A, T, C, G per sequence                   |
| `count_ambiguous_bases`          | `--input`                                 | -                            | Count ambiguous (N) bases                       |
| `filter_by_length`               | `--input`, `--min_length`, `--max_length` | -                            | Keep sequences within length range              |
| `filter_avg_quality`             | `--input`, `--min_avg_quality`            | `--phred_offset`             | Filter by average quality                       |
| `filter_min_quality`             | `--input`, `--min_quality`                | `--phred_offset`             | Filter by minimum base quality                  |
| `trim_records_by_quality`        | `--input`, `--min_quality`                | `--phred_offset`             | Trim low-quality ends of sequences              |
| `filter_gc_content`              | `--input`, `--min_gc`, `--max_gc`         | -                            | Keep sequences within GC% range                 |
| `filter_ambiguous_sequences`     | `--input`, `--max_n_percent`              | -                            | Remove sequences with high N%                   |
| `remove_duplicate_sequences`     | `--input`                                 | -                            | Remove duplicate sequences                      |
| `detect_phred_encoding`          | `--input`                                 | -                            | Detect Phred encoding of FASTQ                  |
| `calculate_average_quality_all`  | `--input`                                 | `--phred_offset`             | Average quality across all sequences            |
| `calculate_quality_position`     | `--input`                                 | `--phred_offset`             | Average quality per position                    |
| `calculate_average_quality_seq`  | `--input`                                 | `--phred_offset`             | Average quality per sequence                    |
| `get_quality_distribution`       | `--input`                                 | `--phred_offset`             | Count occurrences of each quality score         |
| `find_low_quality_positions`     | `--input`, `--threshold`                  | `--phred_offset`             | Positions below quality threshold               |
| `get_quality_stats`              | `--input`, `--threshold`                  | `--phred_offset`             | Summary statistics of quality scores            |
| `plot_length_histogram`          | `--input`                                 | `--output`                   | Histogram of sequence lengths                   |
| `plot_length_distribution`       | `--input`, `--bin_size`                   | `--output`                   | Binned length distribution                      |
| `plot_cumulative_length`         | `--input`                                 | `--output`                   | Cumulative sequence length plot                 |
| `plot_quality_histogram`         | `--input`                                 | `--phred_offset`, `--output` | Histogram of quality scores                     |



## Author
### Karolina Pajdak
Bioinformatics Portfolio Project – August 2025

github.com/karo-pajdak
