

import pytest
from align import sequence_type

def test_sequence_type_dna():
    assert sequence_type("ACGT") == "DNA"

def test_sequence_type_rna():
    assert sequence_type("ACGU") == "RNA"

def test_sequence_type_protein():
    assert sequence_type("ACDEFGHIK") == "AA"

def test_sequence_type_invalid():
    with pytest.raises(SystemExit):
        sequence_type("ACGT@")


from align import score

def test_score_match_no_matrix():
    assert score("A", "A", Matrix=None, match=2, mismatch=-1) == 2

def test_score_mismatch_no_matrix():
    assert score("A", "T", Matrix=None, match=2, mismatch=-1) == -1

from align import nw_score

def test_nw_score_simple():
    seq1 = "A"
    seq2 = "A"

    result = nw_score(seq1, seq2, gap_penalty=-1, Matrix=None, match=1, mismatch=-1)

    assert result[-1] == 1

from align import needleman_wunsch, global_alignment

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

def fastaread_test():
    header, seq = fastaread('test_filer/3filer.fasta')

    with pytest.raises(SystemExit):
        fastaread('test_filer/3filer.fasta')









