import sys

# Funkcja do wczytywania sekwencji z pliku FASTA
def read_fasta(filename):
    sequences = []
    with open(filename, 'r') as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ""
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

# Funkcja Needleman-Wunsch Alignment
def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    m = len(seq1)
    n = len(seq2)

    # Inicjalizacja macierzy scoringowej
    init_mat = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    # Wypełnianie pierwszego wiersza i kolumny karami za przerwy
    for i in range(m + 1):
        init_mat[i][0] = i * gap
    for j in range(n + 1):
        init_mat[0][j] = j * gap

    # Wypełnianie macierzy scoringowej
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diag = init_mat[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            score_up = init_mat[i - 1][j] + gap
            score_left = init_mat[i][j - 1] + gap
            init_mat[i][j] = max(score_diag, score_up, score_left)

    # Backtracking - znajdowanie najlepszego dopasowania
    seq1_align = ""
    seq2_align = ""
    i, j = m, n
    while i > 0 or j > 0:
        score_current = init_mat[i][j]
        score_diag = init_mat[i - 1][j - 1] if i > 0 and j > 0 else None
        score_up = init_mat[i - 1][j] if i > 0 else None
        score_left = init_mat[i][j - 1] if j > 0 else None

        if i > 0 and j > 0 and score_current == score_diag + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
            seq1_align += seq1[i - 1]
            seq2_align += seq2[j - 1]
            i -= 1
            j -= 1
        elif i > 0 and score_current == score_up + gap:
            seq1_align += seq1[i - 1]
            seq2_align += "-"
            i -= 1
        elif j > 0 and score_current == score_left + gap:
            seq1_align += "-"
            seq2_align += seq2[j - 1]
            j -= 1

    # Odwracanie wyrównanych sekwencji
    seq1_align = seq1_align[::-1]
    seq2_align = seq2_align[::-1]

    # Stworzenie ciągu dopasowań
    match_string = ""
    for i in range(len(seq1_align)):
        if seq1_align[i] == seq2_align[i]:
            match_string += "|"
        elif seq1_align[i] == "-" or seq2_align[i] == "-":
            match_string += " "
        else:
            match_string += "*"

    # Wynik dopasowania
    alignment_score = init_mat[m][n]

    # Zwrócenie wyników do zapisu
    result = f"Aligned Sequences:\n{seq1_align}\n{match_string}\n{seq2_align}\nAlignment score: {alignment_score}\n"
    return result

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python nw.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    sequences = read_fasta(fasta_file)

    if len(sequences) != 2:
        print("Error: The FASTA file must contain exactly two sequences.")
        sys.exit(1)

    seq1, seq2 = sequences
    result = needleman_wunsch(seq1, seq2)

    # Zapisanie wyników do pliku wynikowego
    with open("alignment_result.txt", "w") as output_file:
        output_file.write(result)

    print("Alignment completed. Results saved to 'alignment_result.txt'.")
