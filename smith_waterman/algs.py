import pandas as pd
import numpy as np
from Bio import SeqIO
import itertools

def sw(seq1, seq2, method, gap_start = 10, gap_extend = 2):
    M, X, Y = initialize_scoring_matrices(seq1, seq2)
    scoring = read_scoring_matrix(method)
    M, X, Y , PM, PX, PY= fill_matrix(M, X, Y, scoring, seq1, seq2, gap_start, gap_extend)
    SMat, T, score = choose_trace_start(M, X, Y, PM, PX, PY)
    print('hi')
    print(seq1)
    align1, sym, align2, identity = traceback(T, seq1, seq2)
    return score, align1, align2, sym


def score(seq1, seq2, method, gap_start = 10, gap_extend = 2):
    M, X, Y = initialize_scoring_matrices(seq1, seq2)
    scoring = read_scoring_matrix(method)
    M, X, Y , PM, PX, PY= fill_matrix(M, X, Y, scoring, seq1, seq2, gap_start, gap_extend)
    SMat, T, score = choose_trace_start(M, X, Y, PM, PX, PY)
    return score

def roc():
    return None


def initialize_scoring_matrices(a, b):
    M = np.zeros((len(a) + 1, len(b) + 1), np.int)
    X = np.zeros((len(a) + 1, len(b) + 1), np.int)
    Y = np.zeros((len(a) + 1, len(b) + 1), np.int)
    return M, X, Y

def fill_matrix(M, X, Y, scoring, a, b, gap_start, gap_extend):
    PM = M.copy()
    PX = X.copy()
    PY = Y.copy()
    for i, j in itertools.product(range(1, M.shape[0]), range(1, M.shape[1])):
        match_score = scoring.loc[a[i-1], b[j-1]]
        m1 = [3, M[i-1, j-1]]
        m2 = [2, X[i-1, j-1]]
        m3 = [1, Y[i-1, j-1]]
        m4 = [0,0]
        ms= [m1, m2, m3, m4]
        mat, score = max(ms, key=lambda item: item[1])
        M[i,j] = match_score + score
        PM[i,j] = mat

        x1 = [3, gap_start + gap_extend + M[i,j-1]]
        x2 = [3, gap_extend + X[i, j-1]]
        x3 = [1, gap_start + gap_extend + Y[i, j-1]]
        x4 = [0,0]
        xs = [x1, x2, x3, x4]
        mat, score = max(xs, key=lambda item: item[1])
        X[i,j] = score
        PX[i,j] = mat


        y1 = [3, gap_start + gap_extend + M[i-1,j]]
        y2 = [2, gap_start + gap_extend + X[i-1,j]]
        y3 = [1, gap_extend + Y[i-1,j]]
        y4 = [0,0]
        ys = [y1, y2, y3, y4]
        mat, score = max(ys, key=lambda item: item[1])
        Y[i,j] = score
        PY[i,j] = mat
    return M, X, Y, PM, PX, PY



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


def choose_trace_start(M, X, Y, PM, PX, PY):
    mi, mj = get_max_loc(M)
    M_max = M[mi,mj]
    xi, xj = get_max_loc(X)
    X_max = X[xi,xj]
    yi, yj = get_max_loc(Y)
    Y_max = Y[yi,yj]

    m=[[M, PM, M_max], [X, PX, X_max], [Y, PX, Y_max]]

    return max(m, key=lambda item: item[2])



def traceback(T, s1, s2):
    i, j = get_max_loc(T)
    align1, align2 = '', ''
    while T[i][j] != 0:
        if T[i][j] == 3:
            a1 = s1[i-1]
            a2 = s2[j-1]
            i -= 1
            j -= 1
        elif T[i][j] == 2:
            a1 = '-'
            a2 = s2[j-1]
            j -= 1
        elif T[i][j] == 1:
            a1 = s1[i-1]
            a2 = '-'
            i -= 1
        align1 += a1
        align2 += a2

    align1 = align1[::-1]
    align2 = align2[::-1]
    sym = ''
    iden = 0
    for i in range(len(align1)):
        a1 = align1[i]
        a2 = align2[i]
        if a1 == a2:
            sym += a1
            iden += 1
        elif a1 != a2 and a1 != '-' and a2 != '-':
            sym += '.'
        elif a1 == '-' or a2 == '-':
            sym += '.'

    identity = iden / len(align1) * 100
    print('Identity = %f percent' % identity)
    print(align1)
    print(sym)
    print(align2)

    return align1, sym, align2, identity
