from collections import defaultdict
from parser import SequenceRecord
from typing import List, Optional, Dict, Union

def resolve_phred_offset(records: List[SequenceRecord], offset: Optional[int]) -> int:
    """Returns a valid phred offset (33 or 64), detecting it from records if not provided.
Raises ValueError if detection fails or is ambiguous."""
    if offset is None:
        offset = detect_phred_encoding(records)
    if offset not in [33, 64]:
        raise ValueError("Cannot compute stats due to phred issue.")
    return offset

def phred_to_numeric(quality_string: str, phred_offset: int) -> List[int]:
    """Converts sequence quality score from ASCII to numeric based on phred_offset."""
    if phred_offset not in [33, 64]:
        raise ValueError("Invalid phred offset")
    numeric_values = []
    for char in quality_string:
        numeric_values.append(ord(char) - phred_offset)
    return numeric_values

def numeric_to_phred(numeric_scores: List[int], phred_offset: int) -> str:
    """Converts sequence quality score from numeric to ASCII based on phred_offset."""
    if phred_offset not in [33, 64]:
        raise ValueError("Invalid phred offset")
    phred = ""
    for num in numeric_scores:
        phred += (chr(num+phred_offset))
    return phred

def detect_phred_encoding(records: List[SequenceRecord]) -> Union[str,int]:
    """Detects encoding based on quality sequence in records list. Returns 33, 64, or unknown/ambiguous."""
    min_value = 130 #arbitrary high number
    max_value = 0 
    for seq in records:
        for char in seq.quality:
            val = ord(char)
            min_value = min(min_value, val)
            max_value = max(max_value, val)

    if 33 <= min_value <= 73 and max_value <= 73:
        return 33
    elif 64 <= min_value and max_value <= 104:
        return 64
    elif 64 <= min_value <= 73 and max_value <= 73:
        return "Ambiguous (Phred+33 or Phred+64)"
    else:
        return "Unknown"

def calculate_average_quality(
        records: List[SequenceRecord], 
        phred_offset: Optional[int]) -> float:
    """Calculates the average quality of all sequences in a list of records."""
    phred_offset = resolve_phred_offset(records, phred_offset)
    total_score = 0
    total_bases = 0 
    for seq in records:
        quality = seq.quality
        quality_numeric = phred_to_numeric(quality, phred_offset)
        total_score += sum(quality_numeric)
        total_bases += len(quality_numeric)
    if total_bases == 0:
        raise ValueError("Cannot compute on empty file.")
    average = round((total_score / total_bases),1)
    return average

def analyze_quality_by_position(
        records: List[SequenceRecord], 
        phred_offset: Optional[int]) -> List[float]:
    """Calculates the average quality score for each position across all sequences, returning a list of averages."""
    phred_offset = resolve_phred_offset(records, phred_offset)
    position_sums = {}
    position_counts = {}
    for seq in records:
        quality = seq.quality
        numeric = phred_to_numeric(quality, phred_offset)
        for index, score in enumerate(numeric):
            position_sums[index] = position_sums.get(index, 0) + score
            position_counts[index] = position_counts.get(index, 0) + 1
    average_per_position = []
    for i in sorted(position_sums.keys()):
        average_per_position.append(round(position_sums[i] / position_counts[i],1))
    return average_per_position

def get_average_quality_per_sequence(
        records: List[SequenceRecord], 
        phred_offset: Optional[int]) -> Dict[str, float]:
    """Calculates the average quality score per sequence in a list of sequence records, returning a dictionary of averages."""
    phred_offset = resolve_phred_offset(records, phred_offset)
    average_quality = {}
    for seq in records:
        quality = seq.quality
        if quality:
            numeric = phred_to_numeric(quality, phred_offset)
            average_seq = round(sum(numeric) / len(numeric),2)
        else:
            average_seq = 0
        average_quality[seq.id] = average_seq
    return average_quality

def get_quality_distribution(
        records: List[SequenceRecord], 
        phred_offset: Optional[int]) -> Dict[int, int]:
    """Counts the occurrences of each quality score in a list of sequence records."""
    phred_offset = resolve_phred_offset(records, phred_offset)
    distribution = defaultdict(int)
    for seq in records:
        quality = seq.quality
        numeric = phred_to_numeric(quality, phred_offset)
        for score in numeric:
            distribution[score] += 1
    return dict(sorted(distribution.items())) 

def find_low_quality_positions(
        records: List[SequenceRecord], 
        threshold: int, phred_offset: Optional[int]) -> Dict[int, List[int]]:
    """Finds the positions of each read below given threshold in a list of sequence records."""
    phred_offset = resolve_phred_offset(records, phred_offset)
    low_quality = {}
    for seq in records:
        quality = seq.quality
        seq_id = seq.id
        numeric = phred_to_numeric(quality, phred_offset)
        low_quality_positions = []
        for index, score in enumerate(numeric):
            if score < threshold:
                low_quality_positions.append(index)
        low_quality[seq_id] = low_quality_positions
    return low_quality

def get_quality_stats(
    records: List[SequenceRecord],
    phred_offset: Optional[int],
    threshold: int
) -> Dict[str, Union[float, int, str, Dict[int, int], List[float], int]]:
    """Returns a dictionary summarizing quality stats."""
    phred_offset = resolve_phred_offset(records, phred_offset)
    overall_avg = calculate_average_quality(records, phred_offset)
    avg_by_position = analyze_quality_by_position(records, phred_offset)
    distribution = get_quality_distribution(records, phred_offset)
    low_quality = find_low_quality_positions(records, threshold, phred_offset)
    return {
        "phred_encoding": f"Phred+{phred_offset}",
        "average_quality_overall": round(overall_avg, 2),
        "average_quality_by_position": avg_by_position,
        "quality_score_distribution": distribution,
        "low_quality_reads": low_quality
    }