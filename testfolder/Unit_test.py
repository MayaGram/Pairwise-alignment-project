

import pytest
import sys
from align import sequence_type


#test for sequence type function
def test_sequence_type_dna():
    assert sequence_type("ACGT") == "DNA"

def test_sequence_type_rna():
    assert sequence_type("ACGU") == "RNA"

def test_sequence_type_protein():
    assert sequence_type("ACDEFGHIK") == "AA"

#invalid sequence type
def test_sequence_type_invalid():
    with pytest.raises(SystemExit):
        sequence_type("ACGT@")


from align import score


#test for score function with a simple match and mismatch
def test_score_match_no_matrix():
    assert score("A", "A", Matrix=None, match=2, mismatch=-1) == 2


def test_score_mismatch_no_matrix():
    assert score("A", "T", Matrix=None, match=2, mismatch=-1) == -1

from align import nw_score

#test for nw_score function with a simple match and mismatch
def test_nw_score_simple():
    seq1 = "A"
    seq2 = "A"

    result = nw_score(seq1, seq2, gap_penalty=-1, Matrix=None, match=1, mismatch=-1)

    assert result[-1] == 1

from align import needleman_wunsch, global_alignment

#test for needleman_wunsch function with a simple case
def test_needleman_wunsch_simple():
    align1, align2 = needleman_wunsch(
        "GATTACA",
        "GCATGCU",
        Matrix=None,
        gap_penalty=-1,
        match=1,
        mismatch=-1
    )

    assert len(align1) == len(align2)
    assert align1.replace("-", "") == "GATTACA"
    assert align2.replace("-", "") == "GCATGCU"

#test for global_alignment function with a simple case
def test_global_alignment_score():
    total_score, align1, align2 = global_alignment(
        "AAA",
        "AAA",
        Matrix=None,
        gap_penalty=-1,
        match=1,
        mismatch=-1
    )

    assert total_score == 3
    assert align1 == "AAA"
    assert align2 == "AAA"

from align import fastaread

#test for fastaread function with an invalid file (3 files in one txtx file")
def fastaread_test():
    header, seq = fastaread('test_filer/3filer.fasta')

    with pytest.raises(SystemExit):
        fastaread('test_filer/3filer.fasta')

from align import check_command_line

#test for check_command_line function with valid and invalid arguments
def test_check_command_line_valid():
    args=['python3','test_files/DNA.Fasta', 'DNA','global', 'no']

    filename, type_seq, customs, settings, alignment_type = check_command_line(args)
    assert filename == 'test_files/DNA.Fasta'
    assert type_seq == 'DNA'
    assert customs == 'no'
    assert settings == None
    assert alignment_type == 'global'


def test_check_command_line_invalid():
    args=['python3','test_files/DNA.Fasta','global', 'no']

    with pytest.raises(SystemExit):
        check_command_line(args)

from align import read_matrix

#test for read_matrix function with a valid matrix file
def test_read_matrix_valid(tmp_path, monkeypatch):
    matrix_file = tmp_path / "dna_matrix.txt"
    matrix_file.write_text(
        """A C G T
A 2 -1 -1 -1
C -1 2 -1 -1
G -1 -1 2 -1
T -1 -1 -1 2
"""
    )

    monkeypatch.setattr(sys, "argv", ["align.py", "DNA.fasta", "DNA",'global', "yes"])



    result = read_matrix(str(matrix_file))

    expected = {
        "A": {"A": 2, "C": -1, "G": -1, "T": -1},
        "C": {"A": -1, "C": 2, "G": -1, "T": -1},
        "G": {"A": -1, "C": -1, "G": 2, "T": -1},
        "T": {"A": -1, "C": -1, "G": -1, "T": 2},
    }

    assert result == expected


from align import smith_waterman

#test for smith_waterman function with a simple case
def test_smith_waterman_simple():
    max_score, align1, align2 = smith_waterman(
        "AAAA",
        "AAAT",
        Matrix=None,
        gap_penalty=-1,
        match=1,
        mismatch=-1
    )

    assert len(align1) == len(align2)
    assert align1.replace("-", "") in "AAAA"
    assert align2.replace("-", "") in "AAAT"
    assert max_score == 3

#test for smith_waterman function with no match
def test_smith_waterman_no_match():
    max_score, align1, align2 = smith_waterman(
        "AAAA",
        "TTTT",
        Matrix=None,
        gap_penalty=-1,
        match=1,
        mismatch=-1
    )

    assert align1 == ""
    assert align2 == ""
    assert max_score == 0

from align import local_alignment

#test for local_alignment function with a simple case
def test_local_alignment_simple():
    max_score, align1, align2 = local_alignment(
        "AAAA",
        "AAAT",
        Matrix=None,
        gap_penalty=-1,
        match=1,
        mismatch=-1
    )

    assert len(align1) == len(align2)
    assert align1.replace("-", "") in "AAAA"
    assert align2.replace("-", "") in "AAAT"
    assert max_score == 3

#test for local_alignment function with no match
def test_local_alignment_no_match():
    max_score, align1, align2 = local_alignment(
        "AAAA",
        "TTTT",
        Matrix=None,
        gap_penalty=-1,
        match=1,
        mismatch=-1
    )

    assert align1 == ""
    assert align2 == ""
    assert max_score == 0










