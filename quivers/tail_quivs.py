from quivers.knot_quiver import colored_homfly_vectors_and_quiver
import sympy as sp

def flip_vectors_and_quiver(S, A, Q):
    flip_Q = sp.Matrix.eye(Q.rows) - sp.Matrix.ones(Q.rows) - Q
    flip_A = [-a for a in A]
    flip_S = [-s for s in S]
    return flip_S, flip_A, flip_Q

def normalize_vectors_and_quiver(S, A, Q):
    diag = [Q[i, i] for i in range(Q.rows)]
    idx = diag.index(min(diag))
    new_Q = Q - Q[idx, idx] * sp.Matrix.ones(Q.rows, Q.cols)
    new_A = [a - A[idx] for a in A]
    new_S = [s - S[idx] for s in S]
    return new_S, new_A, new_Q

def permute_key(M, i):
    n = M.rows
    right_vals = [M[i, i + k] for k in range(1, n - i)]
    padding = [-sp.oo] * (n - 1 - len(right_vals))
    return (M[i, i], *right_vals, *padding)

def permute_vectors_and_quiver(S, A, Q):
    diag = [Q[i, i] for i in range(Q.rows)]
    perm = sorted(range(len(diag)), key=lambda i: permute_key(Q, i))
    Q_sorted = Q[perm, perm] 
    A_sorted = [A[p] for p in perm]
    S_sorted = [S[p] for p in perm]
    return S_sorted, A_sorted, Q_sorted

def sub_vectors_and_quiver(S, A, Q):
    "must already be sorted/permuted"
    idx = [i for i in range(1, Q.cols) if Q[i, 0] == 0]
    sub_Q = Q[idx, idx]
    sub_A = [A[i] for i in idx]
    sub_S = [S[i] for i in idx]
    return sub_S, sub_A, sub_Q

def colored_homfly_tail_vectors_and_quiver(u, v):
    S, A, Q = colored_homfly_vectors_and_quiver(u, v)
    norm_S, norm_A, norm_Q = normalize_vectors_and_quiver(S, A, Q)
    sort_S, sort_A, sort_Q = permute_vectors_and_quiver(norm_S, norm_A, norm_Q)
    sub_S, sub_A, sub_Q = sub_vectors_and_quiver(sort_S, sort_A, sort_Q)
    return sub_S, sub_A, sub_Q

def colored_homfly_head_vectors_and_quiver(u, v):
    S, A, Q = colored_homfly_vectors_and_quiver(u, v)
    flip_S, flip_A, flip_Q = flip_vectors_and_quiver(S, A, Q)
    norm_S, norm_A, norm_Q = normalize_vectors_and_quiver(flip_S, flip_A, flip_Q)
    sort_S, sort_A, sort_Q = permute_vectors_and_quiver(norm_S, norm_A, norm_Q)
    sub_S, sub_A, sub_Q = sub_vectors_and_quiver(sort_S, sort_A, sort_Q)
    return sub_S, sub_A, sub_Q