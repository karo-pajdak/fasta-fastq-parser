import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import numpy as np
import os 
from stats import calculate_length_distributions, calculate_gc_content, get_nucleotide_composition
from quality import get_quality_distribution, resolve_phred_offset, analyze_quality_by_position
from parser import SequenceRecord
from typing import List, Optional

def freedman_diaconis_bins(data):
    """Supplementary function to determine optimal number of bins for a histogram, if not provided by user."""
    data = np.asarray(data)
    q25, q75 = np.percentile(data, [25, 75])
    iqr = q75 - q25
    n = len(data)
    bin_width = 2 * iqr / n**(1/3)
    if bin_width == 0:
        return 10  
    bins = int(np.ceil((data.max() - data.min()) / bin_width))
    return bins

def plot_length_histogram(lengths: List[int], output_path: str) -> None:
    """Generates and saves a histogram of sequence lengths."""
    output_path = os.path.expanduser(output_path)
    folder = os.path.dirname(output_path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, edgecolor='black')
    plt.title("Sequence Length Distribution")
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_length_distribution(lengths: List[int], output_path: str, bin_size: int) -> None:
    """Generates and saves a binned bar graph of sequence lengths."""
    length_distribution = calculate_length_distributions(lengths, bin_size)
    if not length_distribution:
        print("No lengths provided. Nothing to graph.")
        return
    bins = list(length_distribution.keys())
    counts = list(length_distribution.values())
    output_path = os.path.expanduser(output_path)
    folder = os.path.dirname(output_path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    plt.figure(figsize=(10, 6))
    plt.bar(bins, counts, width=bin_size*0.9, edgecolor='black')
    plt.title(f"Sequence Length Distribution - Bin Size: {bin_size}")
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_cumulative_length(lengths: List[int], output_path: str) -> None:
    """Generates and saves a line plot of cumulative sequence length."""
    if not lengths:
        print("No lengths provided. Nothing to graph.")
        return
    output_path = os.path.expanduser(output_path)
    folder = os.path.dirname(output_path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    sorted_lengths = sorted(lengths)
    cumulative_counts = np.arange(1, len(sorted_lengths) + 1)
    plt.figure(figsize=(10,6))
    plt.plot(sorted_lengths, cumulative_counts, drawstyle='steps-post')
    plt.title("Cumulative Sequence Length Distribution")
    plt.xlabel("Sequence Length")
    plt.ylabel("Cumulative Count")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_quality_histogram(records: List[SequenceRecord], phred_offset: Optional[int], output_path: str) -> None:
    """Generates and saves a histogram of sequence quality scores."""
    if records[0].quality is None:
        raise ValueError("Quality filtering requires FASTQ input with quality scores.")
    phred_offset = resolve_phred_offset(records, phred_offset)
    quality_distribution = get_quality_distribution(records, phred_offset)
    if not quality_distribution:
        print("No sequences provided. Nothing to graph.")
        return
    output_path = os.path.expanduser(output_path)
    folder = os.path.dirname(output_path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    quality_score = sorted(quality_distribution.keys())
    frequency = [quality_distribution[q] for q in quality_score]
    plt.figure(figsize=(10,6))
    plt.bar(quality_score, frequency, edgecolor='black')
    plt.ticklabel_format(style='plain', axis='y')
    ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(axis='x', which='minor', length=4, color='gray') 
    plt.title("Quality Score Distribution")
    plt.xlabel("Quality Score")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_quality_by_position(records: List[SequenceRecord], phred_offset: Optional[int], output_path: str) -> None:
    """Generates and saves a line plot of average quality score per position."""
    if records[0].quality is None:
        raise ValueError("Quality filtering requires FASTQ input with quality scores.")
    position_quality = analyze_quality_by_position(records, phred_offset)
    if not position_quality:
        print("No sequences provides. Nothing to graph.")
    output_path = os.path.expanduser(output_path)
    folder = os.path.dirname(output_path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    positions = list(range(len(position_quality)))
    plt.figure(figsize=(20, 6)) 
    plt.plot(positions, position_quality, marker=".", linestyle='-')
    plt.title("Average Quality Score by Sequence Position")
    plt.xlabel("Sequence Position")
    plt.ylabel("Average Quality Score")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_gc_distribution(records: List[SequenceRecord], bins: Optional[int], output_path: str) -> None:
    """Generates and saves a histogram of GC content percentage."""
    gc_content = calculate_gc_content(records)
    if bins is None:
        bins = freedman_diaconis_bins(gc_content)
    output_path = os.path.expanduser(output_path)
    folder = os.path.dirname(output_path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    plt.figure(figsize=(10, 6))
    plt.hist(gc_content, bins, edgecolor='black')
    plt.title("GC Content Distribution")
    plt.xlabel("GC Content (%)")
    plt.ylabel("Frequency")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_nucleotide_composition(records: List[SequenceRecord], output_path: str) -> None:
    """Generates and saves a bar graph of total counts for each nucleotide."""
    composition_data = get_nucleotide_composition(records)
    output_path = os.path.expanduser(output_path)
    folder = os.path.dirname(output_path)
    if folder:
        os.makedirs(folder, exist_ok=True)
    total_counts = {"A":0,"T":0,"G":0,"C":0}
    for comp in composition_data:
        for nt in total_counts:
            total_counts[nt] += comp.get(nt,0)
    nucleotides = list(total_counts.keys())
    counts = [total_counts[nt] for nt in nucleotides]
    plt.figure(figsize=(8,6))
    plt.bar(nucleotides,counts, color=['#3CB371', '#FF6347', '#FFD700', '#6495ED'])
    plt.ticklabel_format(style='plain', axis='y')
    plt.title("Nucleotide Composition")
    plt.xlabel("Nucleotide")
    plt.ylabel("Total Count")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()