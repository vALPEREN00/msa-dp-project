import numpy as np

def calc_sp_score(c1, c2, c3, match=3, mismatch=-3, gap=-2):
    """Sum of Pairs (SP) mantığı ile 3 karakterin toplam skorunu hesaplar."""
    score = 0
    for p1, p2 in [(c1, c2), (c1, c3), (c2, c3)]:
        if p1 == '-' and p2 == '-':
            continue 
        elif p1 == '-' or p2 == '-':
            score += gap
        elif p1 == p2:
            score += match
        else:
            score += mismatch
    return score

def create_score_matrix(x, y, match=3, mismatch=-3, gap=-2):
    """2 dizi için Needleman-Wunsch skor matrisini oluşturur."""
    s = np.zeros((len(x) + 1, len(y) + 1), dtype=int)

    for i in range(1, len(x) + 1):
        s[i, 0] = s[i-1, 0] + gap
    for j in range(1, len(y) + 1):
        s[0, j] = s[0, j-1] + gap

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            cost = match if x[i-1] == y[j-1] else mismatch
            s[i, j] = max(
                s[i-1, j-1] + cost, 
                s[i-1, j] + gap,    
                s[i, j-1] + gap     
            )
    return s

def align_nw(m, s1, s2, match=3, mismatch=-3, gap=-2):
    """Skor matrisi üzerinden Traceback"""
    se1, se2 = "", ""
    i, j = len(s1), len(s2)

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            cost = match if s1[i-1] == s2[j-1] else mismatch
            if m[i, j] == m[i-1, j-1] + cost:
                se1 = s1[i-1] + se1
                se2 = s2[j-1] + se2
                i -= 1
                j -= 1
            elif m[i, j] == m[i, j-1] + gap:
                se1 = '-' + se1
                se2 = s2[j-1] + se2
                j -= 1
            else:
                se1 = s1[i-1] + se1
                se2 = '-' + se2
                i -= 1
        elif i > 0:
            se1 = s1[i-1] + se1
            se2 = '-' + se2
            i -= 1
        elif j > 0:
            se1 = '-' + se1
            se2 = s2[j-1] + se2
            j -= 1

    return se1, se2

def align_three_sequences(seq1, seq2, seq3):
    L1, L2, L3 = len(seq1), len(seq2), len(seq3)
    s = np.zeros((L1 + 1, L2 + 1, L3 + 1), dtype=int)
    
    for i in range(1, L1 + 1):
        s[i, 0, 0] = s[i-1, 0, 0] + calc_sp_score(seq1[i-1], '-', '-')
    for j in range(1, L2 + 1):
        s[0, j, 0] = s[0, j-1, 0] + calc_sp_score('-', seq2[j-1], '-')
    for k in range(1, L3 + 1):
        s[0, 0, k] = s[0, 0, k-1] + calc_sp_score('-', '-', seq3[k-1])
        
    for i in range(1, L1 + 1):
        for j in range(1, L2 + 1):
            s[i, j, 0] = max(
                s[i-1, j-1, 0] + calc_sp_score(seq1[i-1], seq2[j-1], '-'),
                s[i-1, j, 0] + calc_sp_score(seq1[i-1], '-', '-'),
                s[i, j-1, 0] + calc_sp_score('-', seq2[j-1], '-')
            )
    for i in range(1, L1 + 1):
        for k in range(1, L3 + 1):
            s[i, 0, k] = max(
                s[i-1, 0, k-1] + calc_sp_score(seq1[i-1], '-', seq3[k-1]),
                s[i-1, 0, k] + calc_sp_score(seq1[i-1], '-', '-'),
                s[i, 0, k-1] + calc_sp_score('-', '-', seq3[k-1])
            )
    for j in range(1, L2 + 1):
        for k in range(1, L3 + 1):
            s[0, j, k] = max(
                s[0, j-1, k-1] + calc_sp_score('-', seq2[j-1], seq3[k-1]),
                s[0, j-1, k] + calc_sp_score('-', seq2[j-1], '-'),
                s[0, j, k-1] + calc_sp_score('-', '-', seq3[k-1])
            )

    # 3D Matrisi Doldurma
    for i in range(1, L1 + 1):
        for j in range(1, L2 + 1):
            for k in range(1, L3 + 1):
                s[i, j, k] = max(
                    s[i-1, j-1, k-1] + calc_sp_score(seq1[i-1], seq2[j-1], seq3[k-1]),
                    s[i-1, j-1, k] + calc_sp_score(seq1[i-1], seq2[j-1], '-'),
                    s[i-1, j, k-1] + calc_sp_score(seq1[i-1], '-', seq3[k-1]),
                    s[i, j-1, k-1] + calc_sp_score('-', seq2[j-1], seq3[k-1]),
                    s[i-1, j, k] + calc_sp_score(seq1[i-1], '-', '-'),
                    s[i, j-1, k] + calc_sp_score('-', seq2[j-1], '-'),
                    s[i, j, k-1] + calc_sp_score('-', '-', seq3[k-1])
                )
                
    #Traceback
    align1, align2, align3 = "", "", ""
    i, j, k = L1, L2, L3
    
    while i > 0 or j > 0 or k > 0:
        current = s[i, j, k]
        
        # Öncelik 1: 3'lü çapraz hareket (Kimse boşluk almaz)
        if i > 0 and j > 0 and k > 0 and current == s[i-1, j-1, k-1] + calc_sp_score(seq1[i-1], seq2[j-1], seq3[k-1]):
            align1, align2, align3 = seq1[i-1] + align1, seq2[j-1] + align2, seq3[k-1] + align3
            i -= 1; j -= 1; k -= 1
            
        # Öncelik 2: 2'li hareketler (Biri boşluk alır)
        elif i > 0 and j > 0 and current == s[i-1, j-1, k] + calc_sp_score(seq1[i-1], seq2[j-1], '-'):
            align1, align2, align3 = seq1[i-1] + align1, seq2[j-1] + align2, '-' + align3
            i -= 1; j -= 1
        elif i > 0 and k > 0 and current == s[i-1, j, k-1] + calc_sp_score(seq1[i-1], '-', seq3[k-1]):
            align1, align2, align3 = seq1[i-1] + align1, '-' + align2, seq3[k-1] + align3
            i -= 1; k -= 1
        elif j > 0 and k > 0 and current == s[i, j-1, k-1] + calc_sp_score('-', seq2[j-1], seq3[k-1]):
            align1, align2, align3 = '-' + align1, seq2[j-1] + align2, seq3[k-1] + align3
            j -= 1; k -= 1
            
        # Öncelik 3: 1'li hareketler (İkisi boşluk alır)
        elif i > 0 and current == s[i-1, j, k] + calc_sp_score(seq1[i-1], '-', '-'):
            align1, align2, align3 = seq1[i-1] + align1, '-' + align2, '-' + align3
            i -= 1
        elif j > 0 and current == s[i, j-1, k] + calc_sp_score('-', seq2[j-1], '-'):
            align1, align2, align3 = '-' + align1, seq2[j-1] + align2, '-' + align3
            j -= 1
        elif k > 0 and current == s[i, j, k-1] + calc_sp_score('-', '-', seq3[k-1]):
            align1, align2, align3 = '-' + align1, '-' + align2, seq3[k-1] + align3
            k -= 1

    return align1, align2, align3

if __name__ == "__main__":
    sq1 = "GATTACA"
    sq2 = "GCATGCU"
    sq3 = "GATCA"
    
    matrix = create_score_matrix(sq1, sq2)
    res1, res2 = align_nw(matrix, sq1, sq2)
    print("--- 2'li MSA Sonucu ---")
    print("Dizi 1:", res1)
    print("Dizi 2:", res2)
    
    r1, r2, r3 = align_three_sequences(sq1, sq2, sq3)
    print("\n--- 3'lü MSA Sonucu ---")
    print("Dizi 1:", r1)
    print("Dizi 2:", r2)
    print("Dizi 3:", r3)