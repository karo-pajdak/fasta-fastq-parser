from quality import get_average_quality_per_sequence
from quality import phred_to_numeric
from stats import calculate_gc_content, count_ambiguous_bases


def filter_by_length(records, min_length, max_length):
    filtered = []
    for seq in records:
       if len(seq.sequence) >= min_length and len(seq.sequence) <= max_length:
           filtered.append(seq.sequence)
    return filtered

def get_length_filter_preview(records, min_length, max_length):
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
        "total_sequence" : total,
        "below_min" : below_min,
        "above_max" : above_max,
        "within_range" : within_range,
        "percent_within" : round((within_range/total)*100, 1)
    }
    return length_stats

def filter_by_quality(records, min_avg_quality, phred_offset):
    filtered_records = []
    average_quality_per_seq = get_average_quality_per_sequence(records, phred_offset)
    for seq in records:
        average_seq = average_quality_per_seq.get(seq.id,0)
        if average_seq >= min_avg_quality:
            filtered_records.append(seq)
    return filtered_records

def filter_by_position_quality(records, position_threshold, phred_offset):
    filtered_records = []
    for seq in records:
        quality = seq.quality
        if not quality:
            continue
        numeric = phred_to_numeric(quality, phred_offset)
        if all(score > position_threshold for score in numeric):
            filtered_records.append(seq)
    return filtered_records

def trim_low_quality_ends(sequence, quality, threshold, phred_offset):
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

def trim_records_by_quality(records, threshold, phred_offset):
    trimmed = []
    for seq in records:
        trimmed_seq, trimmed_quality = trim_low_quality_ends(seq.sequence, seq.quality, threshold, phred_offset)
        if trimmed_seq:
            new_record = SequenceRecord(seq.id, trimmed_seq, trimmed_quality)
            trimmed.append(new_record)
    return trimmed

def filter_by_gc_content(records, min_gc, max_gc):
    gc_filtered = []
    gc_list = calculate_gc_content(records)
    for seq, gc in zip(records, gc_list):
        if min_gc <= gc <= max_gc:
            gc_filtered.append(seq)
    return gc_filtered

def filter_ambiguous_sequences(records, max_n_percent):
    ambig_filtered = []
    ambig_stats = count_ambiguous_bases(records)
    for seq, stat in zip(records,ambig_stats):
        if float(stat["percent_ambiguous"]) <= max_n_percent:
            ambig_filtered.append(seq)
    return ambig_filtered

def remove_duplicate_sequences(records):
    seen = set()
    unique_records = []
    for seq in records:
        key = (seq.sequence, seq.quality)
        if key not in seen:
            seen.add(key)
            unique_records.append(seq)
    return unique_records