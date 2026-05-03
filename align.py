import numpy as np
"""Needleman-Wunsch Score Matrix"""
def create_score_matrix(x, y, match=3, mismatch=-3, gap=-2):
    
    s = np.zeros((len(x) + 1, len(y) + 1))

    for i in range(1, len(x) + 1):
        s[i, 0] = s[i-1, 0] + gap
    for j in range(1, len(y) + 1):
        s[0, j] = s[0, j-1] + gap

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            if x[i-1] == y[j-1]:
                a = s[i-1, j-1] + match
            else:
                a = s[i-1, j-1] + mismatch
                
            b = s[i-1, j] + gap
            c = s[i, j-1] + gap

            s[i, j] = max(a, b, c)
            
    return s

# Traceback
def align_nw(m, s1, s2, match=3, mismatch=-3, gap=-2):
    se1 = ""
    se2 = ""
    
    i = len(s1)
    j = len(s2)

    while i > 0 or j > 0:
        if i == 0:
            se1 = '-' + se1
            se2 = s2[j-1] + se2
            j -= 1
        elif j == 0:
            se1 = s1[i-1] + se1
            se2 = '-' + se2
            i -= 1
        else:
            # Çaprazdan gelme (Eşleşme durumu)
            if s1[i-1] == s2[j-1]:
                se1 = s1[i-1] + se1
                se2 = s2[j-1] + se2
                i -= 1
                j -= 1
            # Çaprazdan gelme (Eşleşmeme durumu)
            elif m[i, j] == m[i-1, j-1] + mismatch:
                se1 = s1[i-1] + se1
                se2 = s2[j-1] + se2
                i -= 1
                j -= 1
            # Soldan gelme (1. dizide boşluk)
            elif m[i, j] == m[i, j-1] + gap:
                se1 = '-' + se1
                se2 = s2[j-1] + se2
                j -= 1
            # Yukarıdan gelme (2. dizide boşluk)
            else:
                se1 = s1[i-1] + se1
                se2 = '-' + se2
                i -= 1

    return se1, se2

if __name__ == "__main__":
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    
    matrix = create_score_matrix(seq1, seq2)
    res1, res2 = align_nw(matrix, seq1, seq2)
    
    print(res1)
    print(res2)