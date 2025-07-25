from collections import defaultdict

def phred_to_numeric(quality_string, phred_offset):
    if phred_offset not in [33, 64]:
        raise ValueError("Invalid phred offset")
    numeric_values = []
    for char in quality_string:
        numeric_values.append(ord(char) - phred_offset)
    return numeric_values

def numeric_to_phred(numeric_scores, phred_offset):
    if phred_offset not in [33, 64]:
        raise ValueError("Invalid phred offset")
    quality_string = []
    for num in numeric_scores:
        quality_string.append(chr(num+phred_offset))
    return quality_string

## def detect_phred_encoding

def calculate_average_quality(quality_string, phred_offset):
    quality_numeric = phred_to_numeric(quality_string, phred_offset)
    average = round((sum(quality_numeric) / len(quality_numeric)),1)
    return average

def analyze_quality_by_position(records, phred_offset):
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

def get_quality_distribution(records, phred_offset):
    distribution = defaultdict(int)
    for seq in records:
        quality = seq.quality
        numeric = phred_to_numeric(quality, phred_offset)
        print(numeric)
        for score in numeric:
            distribution[score] += 1
    return dict(sorted(distribution.items())) 




#print(phred_to_numeric("IHHGF", 33))
#print(numeric_to_phred([40, 39, 39, 38, 37], 33))
#print(calculate_average_quality("IHHGF", 33))