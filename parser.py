import gzip
import bz2
import os

class SequenceRecord:
    def __init__(self, id, sequence, quality=None):
        self.id = id
        self.sequence = sequence
        self.quality = quality

class SequenceParser:
    def __init__(self, filename):
        self.filename = filename

    def open_compressed_file(self, filename, mode):
        # open compressed or uncompressed files depending on extension
        try:
            extension = os.path.splitext(filename)[1]
            if extension == ".gz":
                return gzip.open(filename, mode)
            elif extension == ".bz2":
                return bz2.open(filename, mode)
            else:
                return open(filename, mode)
        except FileNotFoundError:
            print(f"Error: File {filename} not found!")
            return None
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return None

    def detect_format(self):
        # open file in read format and check first character
        f = self.open_compressed_file(self.filename, "r")
        if f is None:
            return "Invalid file"
        try:
            with f:
                fast_file = f.readline()
                if not fast_file:
                    return "Invalid file"
                if fast_file[0] == ">":
                    return "FASTA"
                elif fast_file[0] == "@":
                    return "FASTQ"
                else:
                    return "Invalid file"
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return "Invalid file"
    
    def parse_fasta(self, file_handle):
        sequence_id = None
        sequence = []
        # loop through lines and extract sequence id and sequence
        for line in file_handle:
            line = line.strip()
            if line.startswith(">"):
                # already have a sequence id; yield previous sequence id and sequence 
                if sequence_id is not None: 
                    yield SequenceRecord(sequence_id, ''.join(sequence))
                # reset for next sequence
                sequence_id = line[1:] 
                sequence = []
            else:
                # if sequence spans multiple lines, add to previous
                sequence.append(line)
        # yield last sequence id and sequence in file
        if sequence_id is not None:
            yield SequenceRecord(sequence_id, ''.join(sequence))
        
    def parse_fastq(self, file_handle):
        sequence_id = None
        sequence = []
        quality = []
        reading_quality = False
        # loop through lines and extract sequence id, sequence, and quality score
        for line in file_handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("@"):
                # already have a sequence id; yield previous sequence id, sequence, and quality score
                if sequence_id is not None:
                    yield SequenceRecord(sequence_id, ''.join(sequence), ''.join(quality))
                #reset for next sequence
                sequence_id = line[1:] 
                sequence = []   
                quality = []
                reading_quality = False
            elif line.startswith("+"):
                # read quality lines
                reading_quality = True
            elif reading_quality:
                quality.append(line)
            else:
                # if sequence spans multiple lines, add to previous
                sequence.append(line)
        # yield last sequence id and sequence in file
        if sequence_id is not None:
            yield SequenceRecord(sequence_id, ''.join(sequence), ''.join(quality))