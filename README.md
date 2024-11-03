# needleman_wunsch
Ten plik uruchamiany jest poniższą komendą wpisywaną w terminalu w VS Code:
    python nw.py sekwencje.fa
gdzie nw.py jest nazwą pliku z kodem w języku python, a plik sekwencje.fa jest plikiem zawierającym 2 sekwencje w formacie fasta. Jeżeli plik nie zawiera 2 sekwencji program zwróci informację o błedzie - plik musi zawierać 2 sekwencje.
Utworzony zostaje pik wynikowy o nazwie alignment_result.txt który przechowuje informacje o wyniku dopasowania.
Przykładowa zawartość pliku wejściowego:
    > informacje o sekwencji
    ADADGDGDADADCCCGAGGGGADADACCCAD
    > informacje o sekwencji
    ADDAGGGCCACACDCDCAGAGDCDCAAAADCAGDDA
Oczekiwany wynik:
    A D D A G G G C C A C A - C D C D
    | | * * | * | *   | * |   | * | *
    A D A D G D G D - A D A D C C C G
    Score = -3
gdzie "|" oznacza "match", a "*" oznacza "mismatch".
