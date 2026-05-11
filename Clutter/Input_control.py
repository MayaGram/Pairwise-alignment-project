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

def runner():
    filename = check_command_line(sys.argv)
    headers, sequences = fastaread(filename)
    print(headers)
    print(sequences)

    if len(sequences) != 2 or len(headers) != 2:
        usage("The file contains more than two headers and/or sequences")
    
    if sequence_type(sequences[0]) != sequence_type(sequences[1]):
        usage("The two sequences are not the same type")

    print(sequence_type(sequences[0]))
    return None

runner()




