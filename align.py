#!/usr/bin/env python3
import sys

#Defining all IUPAC letters for the diffrent sequence types
DNA_CHARS = set("ANCGTRYSWKMBDHV")
RNA_CHARS = set("ANCGURYSWKMBDHV")
PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
ALLOWED_FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fsa", ".fas")


def usage(error_message):
    """Printing the error in the code and ending the script"""
    print("Fejl:", error_message, file=sys.stderr)
    print("Input should looke like: python3 align.py <fastafil>", file=sys.stderr)
    sys.exit(1)

def check_command_line(argv):
    """Tjek at der kun er 1 argument og at det er en FASTA-fil."""

    #Cheking input length
    if len(argv) != 2:
        raise UsageError("Der skal gives præcis én input-fil.")

    filename = argv[1]
    lower_name = filename.lower()

    if not lower_name.endswith(ALLOWED_FASTA_SUFFIXES):
        raise UsageError("Input-filen skal være en FASTA-fil (.fasta, .fa, .fna, .fsa, .fas).")

    try:
        with open(filename, "r"):
            pass
    except OSError:
        raise UsageError("The file can't be read" + filename)

    return filename

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


def nw_score(seq1, seq2):
    gap_penalty = -1
    match_score = 1
    mismatch_penalty = -1

    prev = [j*gap_penalty for j in range(len(seq2)+1)]

    for i in range(1,len(seq1)+1):
        curr = [i*gap_penalty] + [0]*len(seq2)

        for j in range(1,len(seq2)+1):
            match = prev[j-1] + (match_score if seq1[i-1]==seq2[j-1] else mismatch_penalty)
            delete = prev[j] + gap_penalty
            insert = curr[j-1] + gap_penalty

            curr[j] = max(match,delete,insert)

        prev = curr

    return prev

def hirschberg(seq1, seq2):
    if len(seq1) == 0:
        return "-"*len(seq2), seq2
    if len(seq2) == 0:
        return seq1, "-"*len(seq1)

    if len(seq1) == 1 or len(seq2) == 1:
        return needleman_wunsch(seq1, seq2)

    mid = len(seq1)//2

    scoreL = nw_score(seq1[:mid], seq2)
    scoreR = nw_score(seq1[mid:][::-1], seq2[::-1])

    m = len(seq2)
    split = max(range(m+1), key=lambda j: scoreL[j] + scoreR[m-j])

    left1, left2 = hirschberg(seq1[:mid], seq2[:split])
    right1, right2 = hirschberg(seq1[mid:], seq2[split:])

    return left1 + right1, left2 + right2

def needleman_wunsch(seq1, seq2):
    match, mismatch, gap = 1, -1, -1
    n, m = len(seq1), len(seq2)

    dp = [[0]*(m+1) for _ in range(n+1)]

    for i in range(n+1):
        dp[i][0] = i * gap
    for j in range(m+1):
        dp[0][j] = j * gap

    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = dp[i-1][j] + gap
            insert = dp[i][j-1] + gap
            dp[i][j] = max(match_score, delete, insert)

    i, j = n, m
    a1, a2 = "", ""

    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + (
            match if seq1[i-1] == seq2[j-1] else mismatch
        ):
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


def global_alignment(seq1, seq2):
    align1, align2 = hirschberg(seq1, seq2)
    score = sum(
        1 if a == b else -1 if a != "-" and b != "-" else -1
        for a, b in zip(align1, align2)
    )
    return score, align1, align2

def alignment_printer(aligned):
    print(aligned[0])

    Connecting_string = ""
    for i, char in enumerate(aligned[0]):
        if char == aligned[1][i]:
            Connecting_string += "|"
        elif char == "-" or aligned[1][i] == "-":
            Connecting_string += " "
        else:
            Connecting_string += ":"
    
    print(Connecting_string)
    print(aligned[1])
    return None

def runner():
    filename = check_command_line(sys.argv)
    headers, sequences = fastaread(filename)
    print(headers)
    print(sequences)

    if len(sequences) != 2 or len(headers) != 2:
        usage("The file contains more than two headers and/or sequences")
    
    if sequence_type(sequences[0]) != sequence_type(sequences[1]):
        usage("The two sequences are not the same type")
    
    score, align1, align2 = global_alignment(sequences[0], sequences[1])
    alignment_printer([align1, align2])
    print(score)


    print(sequence_type(sequences[0]))


    return None

runner()