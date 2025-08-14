from quality import get_average_quality_per_sequence, phred_to_numeric, resolve_phred_offset, get_quality_stats
from parser import SequenceRecord
from typing import List, Optional, Dict, Union, Tuple
from stats import calculate_gc_content, count_ambiguous_bases


def filter_by_length(records: List[SequenceRecord], min_length: int, max_length: int) -> List[str]:
    """Filters sequences based on a specified length range."""
    filtered = []
    for seq in records:
       if min_length <= len(seq.sequence) <= max_length:
           filtered.append(seq.sequence)
    return filtered

def get_length_filter_preview(records: List[SequenceRecord], min_length: int, max_length: int) -> Dict[str, Union[int,float]]:
    """Returns stats on how many sequences fall within the given length range."""
    total = 0
    below_min = 0
    above_max = 0
    within_range = 0
    for seq in records:
        total += 1
        if len(seq.sequence) < min_length:
            below_min += 1
        elif len(seq.sequence) > max_length:
            above_max +=1
        else:
            within_range +=1
    length_stats =  {
        "total_sequences" : total,
        "below_min" : below_min,
        "above_max" : above_max,
        "within_range" : within_range,
        "percent_within" : round((within_range/total)*100, 1) if total else 0.0
    }
    return length_stats

def filter_by_quality(records: List[SequenceRecord], min_avg_quality: int, phred_offset: Optional[int]) -> List[SequenceRecord]:
    """Filters sequences by average quality below minimum threshold."""
    if not records:
        raise ValueError("No records provided to filter.")
    if records[0].quality is None:
        raise ValueError("Quality filtering requires FASTQ input with quality scores.")
    phred_offset = resolve_phred_offset(records, phred_offset)
    filtered_records = []
    average_quality_per_seq = get_average_quality_per_sequence(records, phred_offset)
    for seq in records:
        average_seq = average_quality_per_seq.get(seq.id,0)
        if average_seq >= min_avg_quality:
            filtered_records.append(seq)
    return filtered_records

def filter_by_min_quality(records: List[SequenceRecord], min_quality: int, phred_offset: Optional[int]) -> List[SequenceRecord]:
    """Filters out sequences containing any base with a quality score below the given threshold."""
    if not records:
        raise ValueError("No records provided to filter.")
    if records[0].quality is None:
        raise ValueError("Quality filtering requires FASTQ input with quality scores.")
    phred_offset = resolve_phred_offset(records, phred_offset)
    filtered_records = []
    for seq in records:
        quality = seq.quality
        if not quality:
            continue
        numeric = phred_to_numeric(quality, phred_offset)
        if all(score > min_quality for score in numeric):
            filtered_records.append(seq)
    return filtered_records

def trim_low_quality_ends(sequence: str, quality: str, threshold: int, phred_offset: int) -> Tuple[str, str]:
    """Trims low-quality bases from both ends of a sequence based on the threshold."""
    if not sequence or not quality:
        return "", ""
    numeric = phred_to_numeric(quality, phred_offset)
    start = 0
    end = len(numeric) - 1
    while start <= end and numeric[start] < threshold:
        start += 1
    while end >= start and numeric[end] < threshold:
        end -= 1
    if start > end:
        return "",""
    trimmed_seq = sequence[start:end+1]
    trimmed_quality = quality[start:end+1]
    return trimmed_seq, trimmed_quality

def trim_records_by_quality(records: List[SequenceRecord], threshold: int, phred_offset: Optional[int]) -> List[SequenceRecord]:
    """Trims low quality bases from both ends of all sequences in a record."""
    if not records:
        raise ValueError("No records provided to filter.")
    if records[0].quality is None:
        raise ValueError("Quality filtering requires FASTQ input with quality scores.")
    phred_offset = resolve_phred_offset(records, phred_offset)
    trimmed = []
    for seq in records:
        trimmed_seq, trimmed_quality = trim_low_quality_ends(seq.sequence, seq.quality, threshold, phred_offset)
        if trimmed_seq:
            new_record = SequenceRecord(seq.id, trimmed_seq, trimmed_quality)
            trimmed.append(new_record)
    return trimmed

def filter_by_gc_content(records: List[SequenceRecord], min_gc: float, max_gc: float) -> List[SequenceRecord]:
    """Filters sequences in a list of records based on minimum and maximum GC percentage."""
    gc_filtered = []
    gc_list = calculate_gc_content(records)
    for seq, gc in zip(records, gc_list):
        if min_gc <= gc <= max_gc:
            gc_filtered.append(seq)
    return gc_filtered

def filter_ambiguous_sequences(records: List[SequenceRecord], max_n_percent: float) -> List[SequenceRecord]:
    """Filters sequences in a list of records based on ambiguous base (N) percentage."""
    filtered = []
    ambiguous_stats = count_ambiguous_bases(records)
    for seq, stat in zip(records, ambiguous_stats):
        percent_ambiguous = float(stat.get("percent_ambiguous", 0.0))
        if percent_ambiguous <= max_n_percent:
            filtered.append(seq)
    return filtered

def remove_duplicate_sequences(records: List[SequenceRecord]) -> List[SequenceRecord]:
    """Eliminates identical sequences to retain only unique entries. """
    seen = set()
    unique_records = []
    for seq in records:
        key = (seq.sequence, seq.quality)
        if key not in seen:
            seen.add(key)
            unique_records.append(seq)
    return unique_records