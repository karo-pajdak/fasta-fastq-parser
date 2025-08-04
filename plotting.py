import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import numpy as np
import os 
from stats import calculate_length_distributions
from quality import get_quality_distribution
from parser import SequenceRecord
from typing import List

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

def plot_quality_histogram(records: List[SequenceRecord], phred_offset: int, output_path: str) -> None:
    """Generates and saves a histogram of sequence quality scores."""
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