import pandas as pd
import numpy as np
from Bio import SeqIO
import itertools

def sw(seq1, seq2, method, gap_start, gap_extend):
    M, X, Y = initialize_scoring_matrices(seq1, seq2)
    scoring = read_scoring_matrix(method)
    M, X, Y = fill_matrix(M, X, Y, scoring, seq1, seq2, gap_start, gap_extend)
    ali, a_score = traceback(M, X, Y, seq2)
    return ali[0], a_score


def score():
    return None

def roc():
    return None


def initialize_scoring_matrices(a, b):
    M = np.zeros((len(a) + 1, len(b) + 1), np.int)
    X = np.zeros((len(a) + 1, len(b) + 1), np.int)
    Y = np.zeros((len(a) + 1, len(b) + 1), np.int)
    return M, X, Y

def fill_matrix(M, X, Y, scoring, a, b, gap_start, gap_extend):
    for i, j in itertools.product(range(1, M.shape[0]-1), range(1, M.shape[1]-1)):
        match_score = scoring.loc[a[i], b[j]]
        M[i,j] = match_score + max(M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1], 0)

        x1 = gap_start + gap_extend + M[i,j-1]
        x2 = gap_extend + X[i, j-1]
        x3 = gap_start + gap_extend + Y[i, j-1]
        X[i,j] = max(x1, x2, x3, 0)

        y1 = gap_start + gap_extend + M[i-1,j]
        y2 = gap_start + gap_extend + X[i-1,j]
        y3 = gap_extend + Y[i-1,j]
        Y[i,j] = max(y1, y2, y3, 0)
    return M, X, Y

def traceback(M, X, Y, b, b_='', old_i=0):
    H, score = choose_trace_mat(M, X, Y)
    i, j = get_max_loc(H)
    return trace_recursion(H[0:i, 0:j], b, b_, i), score


def read_seq_file(filename):
    seq = list(SeqIO.read(filename, "fasta").seq)
    return seq

def read_scoring_matrix(method):
    scoring = pd.read_table(method, comment='#', sep=' ', skipinitialspace=True, usecols=range(24))
    return scoring.set_index(scoring.columns)


def get_max_loc(H):
    H_flip = np.flip(np.flip(H, 0), 1)
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))
    return i, j


def choose_trace_mat(M, X, Y):
    mi, mj = get_max_loc(M)
    M_max = M[mi,mj]
    xi, xj = get_max_loc(X)
    X_max = X[xi,xj]
    yi, yj = get_max_loc(Y)
    Y_max = Y[yi,yj]

    m=[[M,M_max], [X,X_max], [Y,Y_max]]

    return max(m, key=lambda item: item[1])


def trace_recursion(H, b, b_='', old_i=0):
    i, j = get_max_loc(H)
    if H[i, j] == 0:
        return b_, j
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
    return trace_recursion(H[0:i, 0:j], b, b_, i)
