from msa import align_nw, align_three_sequences, create_score_matrix

print("MSA Kütüphanesi Testi")
seq1 = "GATTACA"
seq2 = "GCATGCU"
seq3 = "GATCA"

print("\n--- 2'li Hizalama Testi ---")
matris = create_score_matrix(seq1, seq2)
res1, res2 = align_nw(matris, seq1, seq2)
print("Dizi 1:", res1)
print("Dizi 2:", res2)

print("\n--- 3'lü Hizalama Testi ---")
r1, r2, r3 = align_three_sequences(seq1, seq2, seq3)
print("Dizi 1:", r1)
print("Dizi 2:", r2)
print("Dizi 3:", r3)