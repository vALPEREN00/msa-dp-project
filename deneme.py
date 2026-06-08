from msa import align_three_sequences
seq1 = "GATTACA"
seq2 = "GCATGCU"
seq3 = "GATCA"

print("Hizalanacak Diziler:")
print(f"1: {seq1}")
print(f"2: {seq2}")
print(f"3: {seq3}\n")

r1, r2, r3, skor = align_three_sequences(seq1, seq2, seq3)

print("--- Hizalama Sonuçları ---")
print(f"Dizi 1: {r1}")
print(f"Dizi 2: {r2}")
print(f"Dizi 3: {r3}")
print(f"\nOptimal Hizalama Skoru: {skor}")