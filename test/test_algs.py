import numpy as np
from smith_waterman import algs

def test_roc():
    return None

def test_smithwaterman():
    score, align1, align2, align = algs.sw("TTTT", "TTTT", "BLOSUM50")
    print('a1')
    print(align1)
    assert(align1 == align2 == align == "TTTT")

    score, align1, align2, align = algs.sw("GCTA", "GCTTTA", "BLOSUM50")
    assert(align1 == "GC--TA" or align1 == "GCT--A")
    assert(align2 =="GCTTTA")

def test_scoring():
    score, align1, align2, align = algs.sw("GCTA", "GCTTTA", "BLOSUM50")
    print(score)
    print(algs.get_score(align1, align2, align, "BLOSUM50"))
    assert(algs.get_score(align1, align2, align, "BLOSUM50") == score)
