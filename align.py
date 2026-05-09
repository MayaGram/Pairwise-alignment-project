#!/usr/bin/env python3
import sys

#Defining all IUPAC letters for the diffrent sequence types
DNA_CHARS = set("ACGTRWSWKMBDHVN")
RNA_CHARS = set("ACGURWSWKMBDHVN")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWY")
ALLOWED_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fsa", ".fas")

def usage(error_message):
    """Printing the error in the code and ending the script"""
    print("Fejl:", error_message)
    print("Input should looke like: python3 align.py <fastafil> <type> <local/global> <customize yes or no>", file=sys.stderr)
    sys.exit(1)

def check_command_line(argv):
    """Check that the command line input is correct."""

    # Checking input length
    if len(argv) != 5:
        usage("Input length is not correct")

    filename = argv[1].lower()

    if not filename.endswith(ALLOWED_FASTA_SUFFIXES):
        usage("Input file must be a FASTA file (.fasta, .fa, .fna, .fsa, .fas).")

    # Getting the sequence type
    type_seq = argv[2].upper()

    if type_seq not in ("DNA", "RNA", "AA"):
        usage("The sequence type is not correct. It should be DNA, RNA, or AA.")
    
    #Type of alignment local or global
    alignment_type = argv[3].lower()
    if alignment_type not in ("global", "local"):
        usage("You must input global or local for alignment type")

    # Customize section
    customs = argv[4].lower()

    if customs not in ("no", "yes"):
        usage("Input must be either 'no' or 'yes' for custom options.")

    if customs == "yes":
        settings = customize(type_seq)
    else:
        settings = None

    return filename, type_seq, customs, settings, alignment_type

def customize(type_seq):
    type_seq = type_seq.upper()

    if type_seq in ("RNA", "DNA"):

        #Checking if 
        print("You must now provide the input for the 3 values. They must be integers.")
        try:                    
            match = int(input("match (mus be >0): "))
            if match <= 0:
                usage("match must be above 0")
            mismatch = int(input("mismatch: "))
            gap_penalty = int(input("gap penalty: "))
            filename_matrix = None
        except ValueError:
            usage("The values provided must be integers.")

        return {
            "match": match,
            "mismatch": mismatch,
            "gap_penalty": gap_penalty,
        }

    elif type_seq == "AA":
        print("You must now provide the 1 integer value and a filename for a substitution matrix.")
        try:
            gap_penalty = int(input("gap penalty: "))
            filename_matrix = input("Input filename for substitution matrix: ").strip()
        except ValueError:
            usage("Gap and gap penalty must be integers.")

        if filename_matrix == "":
            usage("You must provide a filename for the substitution matrix.")

        return {
            "match": None,
            "mismatch": None,
            "gap_penalty": gap_penalty,
            "filename_matrix": filename_matrix
        }

    else:
        usage("Sequence type must be DNA, RNA, or AA.")


def fastaread(filename):
    """Reads a fasta file by given filename and returns a list with headers and a list with sequences"""
    headers = []
    sequences = []
    # Open the file
    with open(filename, "r") as infile:
        for line in infile:
            if line.startswith('>'):
                headers.append(line.rstrip())
                sequences.append('')
            elif len(sequences) == 0:
                continue        # ignore leading non header in file
            else:
                sequences[-1] += ''.join(line.split())
            if len(headers) > 2:
                usage("The amount of sequences are greater than 2")

    return headers, sequences

def sequence_type(sequence:str):
    """This function tjeks if the sequence is DNA, RNA og aminoacid sequence"""

    if set(sequence).issubset(DNA_CHARS):
        return "DNA"
    elif set(sequence).issubset(RNA_CHARS):
        return "RNA"
    elif set(sequence).issubset(PROTEIN_CHARS):
        return "AA"
    
    #Giving the reason why the sequence does not pass
    else:
        for i, char in enumerate(sequence):
            if char not in DNA_CHARS and char not in RNA_CHARS and char not in PROTEIN_CHARS:
                usage(f"Invalid character '{char}' at position {i}")
        usage("The sequence provided is not DNA, RNA or amino acid.")

def score(a, b, Matrix=None, match=1, mismatch=-1):
    if Matrix:
        return Matrix[a][b]
    return match if a == b else mismatch


def nw_score(seq1, seq2, gap_penalty=-1, Matrix=None, match=1, mismatch=-1):
    prev = [j * gap_penalty for j in range(len(seq2)+1)] # Initialize the first row of the DP table

    for i in range(1, len(seq1)+1): # Iterate through each character in seq1
        curr = [i * gap_penalty] + [0]*len(seq2) # Initialize the current row of the DP table

        for j in range(1, len(seq2)+1): # Iterate through each character in seq2
            s = score(seq1[i-1], seq2[j-1], Matrix, match, mismatch) # getting score from substitution matrix

            match_score = prev[j-1] + s 
            delete = prev[j] + gap_penalty
            insert = curr[j-1] + gap_penalty

            curr[j] = max(match_score, delete, insert) # Update the current cell with the maximum score

        prev = curr # Move to the next row by setting curr to prev

    return prev

def hirschberg(seq1, seq2, Matrix=None, gap_penalty=-1, match=1, mismatch=-1):

    #base cases - if one of the sequences is empty, return the other sequence with gaps
    if len(seq1) == 0:
        return "-"*len(seq2), seq2
    if len(seq2) == 0:
        return seq1, "-"*len(seq1)
    
    #base case - if one of the sequences has length 1, use Needleman-Wunsch to align them

    if len(seq1) == 1 or len(seq2) == 1:
        return needleman_wunsch(seq1, seq2, Matrix, gap_penalty, match, mismatch)

    #divide the first sequence into two halves
    mid = len(seq1)//2

    #compute the score of aligning the first half of seq1 with all of seq2 and the second half of seq1 with all of seq2 in reverse
    scoreL = nw_score(seq1[:mid], seq2, gap_penalty=gap_penalty, Matrix=Matrix, match=match, mismatch=mismatch)
    scoreR = nw_score(seq1[mid:][::-1], seq2[::-1], gap_penalty=gap_penalty, Matrix=Matrix, match=match, mismatch=mismatch)

    # Find the optimal split point
    m = len(seq2)
    split = max(range(m+1), key=lambda j: scoreL[j] + scoreR[m-j])

    # Recursively align the two halves of seq1 with the corresponding halves of seq2
    left1, left2 = hirschberg(seq1[:mid], seq2[:split], Matrix, gap_penalty, match, mismatch)
    right1, right2 = hirschberg(seq1[mid:], seq2[split:], Matrix, gap_penalty, match, mismatch)

    return left1 + right1, left2 + right2

def needleman_wunsch(seq1, seq2, Matrix=None, gap_penalty=-1, match=1, mismatch=-1):

    # Initialize the scoring values and the dimensions of the DP table
    gap = gap_penalty
    n, m = len(seq1), len(seq2)
    dp = [[0]*(m+1) for _ in range(n+1)]

    # Initialize the first row and column of the DP table with gap penalties
    for i in range(n+1):
        dp[i][0] = i * gap
    for j in range(m+1):
        dp[0][j] = j * gap

    # Fill the DP table using the scoring scheme and the values from the substitution matrix
    for i in range(1, n+1):
        for j in range(1, m+1):
            s = score(seq1[i-1], seq2[j-1], Matrix, match, mismatch)

            match_score = dp[i-1][j-1] + s
            delete = dp[i-1][j] + gap
            insert = dp[i][j-1] + gap

            dp[i][j] = max(match_score, delete, insert)

    # Backtrack through the DP table to construct the optimal alignment
    i, j = n, m
    a1, a2 = "", ""

    # Backtrack through the DP table to construct the optimal alignment
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + score(seq1[i-1], seq2[j-1], Matrix, match, mismatch):
            a1 += seq1[i-1]
            a2 += seq2[j-1]
            i -= 1
            j -= 1
        elif i > 0:
            a1 += seq1[i-1]
            a2 += "-"
            i -= 1
        else:
            a1 += "-"
            a2 += seq2[j-1]
            j -= 1

    return a1[::-1], a2[::-1]

def read_matrix(filename):
    # Reads a substitution matrix from a file and returns it as a dictionary of dictionaries.
    try:
        with open(filename) as f:
            line = f.readline()
            f.seek(0)
            if '\t' in line:
                lines = [line1.strip().split('\t') for line1 in f if line1.strip()] 
            elif ' ' in line:
                lines = [line1.strip().split() for line1 in f if line1.strip()]
            else:
                usage('can only take matrixes splittet by " " or tap')
    except FileNotFoundError:
        usage(f"Substitution matrix file '{filename}' not found.")


    headers = lines[0]
    Matrix = {}

    # Validate characters in the matrix against the expected character sets for the specified sequence type
    for row in lines[1:]:
        char = row[0]
        if sys.argv[2].upper() == 'AA' and char not in PROTEIN_CHARS:
            usage(f"Invalid character '{char}' in substitution matrix for amino acids.")
        elif sys.argv[2].upper() == 'DNA' and char not in DNA_CHARS:
            usage(f"Invalid character '{char}' in substitution matrix for DNA.")
        elif sys.argv[2].upper() == 'RNA' and char not in RNA_CHARS:
            usage(f"Invalid character '{char}' in substitution matrix for RNA.")

        Matrix[char] = {headers[i]: int(row[i+1]) for i in range(len(headers))}

    return Matrix

def global_alignment(seq1, seq2, Matrix=None, gap_penalty=-1, match=1, mismatch=-1):

    # Use Hirschberg's algorithm to compute the optimal global alignment of seq1 and seq2
    align1, align2 = hirschberg(seq1, seq2, Matrix, gap_penalty, match, mismatch)

    # Calculate the total score of the alignment by summing the scores of each aligned pair of characters, using the provided scoring scheme and substitution matrix if available
    total_score = sum(
        score(a, b, Matrix, match, mismatch)
        if a != "-" and b != "-" else gap_penalty
        for a, b in zip(align1, align2))

    return total_score, align1, align2

def smith_waterman(seq1, seq2, Matrix=None, gap_penalty=-1, match=1, mismatch=-1):
    gap = gap_penalty
    M, N = len(seq1), len(seq2)

    #Here my alignment matrix is getting scores and path matrix is for keeping track of the scores
    alignment_matrix = [[0] * (N + 1) for i in range(M + 1)]
    path_matrix = [[0] * (N + 1) for i in range(M + 1)]

    max_score = 0
    max_position = (0, 0)

    for i in range(1, M + 1):
        for j in range(1, N + 1):

            #Species keep track of the nukleotides in the position
            species1 = seq1[i - 1]
            species2 = seq2[j - 1]

            #Defining the different ways to move in the matix
            diag_score = alignment_matrix[i - 1][j - 1] + score(species1, species2, Matrix, match, mismatch)
            down_score = alignment_matrix[i - 1][j] + gap
            left_score = alignment_matrix[i][j - 1] + gap

            #Kepping track of the max score it can't ever be negative
            current_score = max(0, diag_score, down_score, left_score)
            alignment_matrix[i][j] = current_score

            if current_score == 0:
                path_matrix[i][j] = -1
            elif current_score == diag_score:
                path_matrix[i][j] = 0
            elif current_score == down_score:
                path_matrix[i][j] = 1
            else:
                path_matrix[i][j] = 2

            if current_score > max_score:
                max_score = current_score
                max_position = (i, j)

    i, j = max_position
    a1, a2 = "", ""

    while i > 0 and j > 0 and path_matrix[i][j] != -1:
        if path_matrix[i][j] == 0:
            a1 += seq1[i - 1]
            a2 += seq2[j - 1]
            i -= 1
            j -= 1

        elif path_matrix[i][j] == 1:
            a1 += seq1[i - 1]
            a2 += "-"
            i -= 1

        else:
            a1 += "-"
            a2 += seq2[j - 1]
            j -= 1

    return max_score, a1[::-1], a2[::-1]

def local_alignment(seq1, seq2, Matrix=None, gap_penalty=-1, match=1, mismatch=-1):
    # Use Smith-Waterman to compute the best local alignment of seq1 and seq2
    total_score, align1, align2 = smith_waterman(seq1, seq2, Matrix, gap_penalty, match, mismatch)
    return total_score, align1, align2

def alignment_printer(aligned):
    seq1 = aligned[0]
    seq2 = aligned[1]
    block_size = 60

    for start in range(0, len(seq1), block_size):
        end = start + block_size
        part1 = seq1[start:end]
        part2 = seq2[start:end]

        print(f"{part1}\t{start}-{len(part1)+start}")

        Connecting_string = ""
        for i, char in enumerate(part1):
            if char == part2[i]:
                Connecting_string += "|"
            elif char == "-" or part2[i] == "-":
                Connecting_string += " "
            else:
                Connecting_string += ":"

        print(Connecting_string + "\t")
        print(part2 + "\t")
        print()

    return None

def runner():
    if len(sys.argv) != 5:
        usage("Wrong input length")

    filename, type_seq, customs, settings, alignment_type = check_command_line(sys.argv)
    matrix = None
    gap_penalty = -1
    match = 1
    mismatch = -1

    if customs == "yes":

        if type_seq == "AA":
            gap_penalty = settings["gap_penalty"]
            file_name_matrix = settings["filename_matrix"]
            matrix = read_matrix(file_name_matrix)

        elif type_seq in ("RNA", "DNA"):
            match = settings["match"]
            mismatch = settings["mismatch"]
            gap_penalty = settings["gap_penalty"]
        
    if customs == "no" and type_seq == "AA":
        matrix = read_matrix("blosum62.txt")

    headers, sequences = fastaread(filename)

    if len(sequences) != 2 or len(headers) != 2:
        usage("The file contains more than two headers and/or sequences")

    if sequence_type(sequences[0]) != sequence_type(sequences[1]):
        usage("The two sequences are not the same type")

    if alignment_type == "global":
        total_score, align1, align2 = global_alignment(
            sequences[0],
            sequences[1],
            matrix,
            gap_penalty,
            match,
            mismatch,
        )
    elif alignment_type == "local":
        total_score, align1, align2 = local_alignment(
            sequences[0],
            sequences[1],
            matrix,
            gap_penalty,
            match,
            mismatch,
        )

    print("=" * 60)
    print(f"Result for pairwise alignment of {filename}")
    print(f"The type of sequences alignt was {type_seq}")
    print("=" * 60)
    print(f"Sequence 1:{headers[0]}")
    print("Was aligned with")
    print(f"Sequence 2:{headers[1]}")
    print("=" * 60)
    print(f"The totalt alignment score was: {total_score}")
    print("Alignment:")
    alignment_printer([align1, align2])

if __name__ == "__main__":
    runner()