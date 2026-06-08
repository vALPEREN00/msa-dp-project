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

def align_three_sequences(seq1, seq2, seq3):
    """3D Dinamik Programlama matrisi ile 3 diziyi hizalar."""
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
    
    optimal_score = s[L1, L2, L3]
                
    align1, align2, align3 = "", "", ""
    i, j, k = L1, L2, L3
    
    while i > 0 or j > 0 or k > 0:
        current = s[i, j, k]
        
        # 1: 3'lü çapraz hareket (Kimse boşluk almaz)
        if i > 0 and j > 0 and k > 0 and current == s[i-1, j-1, k-1] + calc_sp_score(seq1[i-1], seq2[j-1], seq3[k-1]):
            align1, align2, align3 = seq1[i-1] + align1, seq2[j-1] + align2, seq3[k-1] + align3
            i -= 1; j -= 1; k -= 1
            
        # 2: 2'li hareketler (Biri boşluk alır)
        elif i > 0 and j > 0 and current == s[i-1, j-1, k] + calc_sp_score(seq1[i-1], seq2[j-1], '-'):
            align1, align2, align3 = seq1[i-1] + align1, seq2[j-1] + align2, '-' + align3
            i -= 1; j -= 1
        elif i > 0 and k > 0 and current == s[i-1, j, k-1] + calc_sp_score(seq1[i-1], '-', seq3[k-1]):
            align1, align2, align3 = seq1[i-1] + align1, '-' + align2, seq3[k-1] + align3
            i -= 1; k -= 1
        elif j > 0 and k > 0 and current == s[i, j-1, k-1] + calc_sp_score('-', seq2[j-1], seq3[k-1]):
            align1, align2, align3 = '-' + align1, seq2[j-1] + align2, seq3[k-1] + align3
            j -= 1; k -= 1
            
        # 3: 1'li hareketler (İkisi boşluk alır)
        elif i > 0 and current == s[i-1, j, k] + calc_sp_score(seq1[i-1], '-', '-'):
            align1, align2, align3 = seq1[i-1] + align1, '-' + align2, '-' + align3
            i -= 1
        elif j > 0 and current == s[i, j-1, k] + calc_sp_score('-', seq2[j-1], '-'):
            align1, align2, align3 = '-' + align1, seq2[j-1] + align2, '-' + align3
            j -= 1
        elif k > 0 and current == s[i, j, k-1] + calc_sp_score('-', '-', seq3[k-1]):
            align1, align2, align3 = '-' + align1, '-' + align2, seq3[k-1] + align3
            k -= 1

    return align1, align2, align3, optimal_score