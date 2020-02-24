import numpy as np
from smith_waterman import algs

def test_roc():
    return None

def test_smithwaterman():
    score, align1, align2, align = algs.sw("TTTT", "TTTT", "BLOSUM50")
    print(align1)
    assert(align1 == align2 == align == "TTTT")
    
    score, align1, align2, align = algs.sw("GCTA", "GCTTTA", "BLOSUM50")
    assert(align1 == "GC--TA" or assign1 == "GCT--A")
    assert(align2 =="GCTTTA")

def test_scoring():
    assert(algs.score("AAAAAAA", "AAAAAAA", "BLOSUM50") == 35)
    assert(algs.score("AAAAAAA", "GGGGGGG", "BLOSUM50") == 0)
