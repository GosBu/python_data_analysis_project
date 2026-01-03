# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 15:30:26 2025

@author: magor
"""
# Projekt końcowy

# Otwórz plik z sekwencjami DNA
from Bio import SeqIO

for seq_records in SeqIO.parse("sekwencje.txt", "fasta"):
    print("\nID:", seq_records.id)
    print("Sekwencja:", seq_records.seq)
    print("Długosć sekwencji:", len(seq_records))

try:
    with open("sekwencje.txt", "r") as plik:
        zawartosc = plik.read()
        print(zawartosc)
except FileNotFoundError:
    print("Nie znaleziono pliku z sekwencją DNA")
    
# Zapisanie sekwencji z pliku jako lista
lista_sekwencji = []

for seq_records in SeqIO.parse("sekwencje.txt", "fasta"):
    lista_sekwencji.append(seq_records.seq)


# Przypisanie sekwencji do oddzielnych zmiennych
dna_seq1 = lista_sekwencji[0]
dna_seq_ostatnia = lista_sekwencji[-1]


# Pobranie sekwencji użytkownika po ID z NCBI
from Bio import Entrez

mail_uzytkownika = input("Podaj swojego maila: ").strip()

if mail_uzytkownika == "":
    print("Nie podałes maila. Podanie maila jest wymagane przez NCBI. Proszę podaj mail: ")

else:
    print(f"\nDziękuję za podanie maila: {mail_uzytkownika}")
    ID_sekwencji = input("Podaj ID z NCBI sekwencji DNA, którą chcesz dodać do pliku: ").strip()

# Pobranie sekwencji DNA podanej przez użytkownika
try:
    handle = Entrez.efetch(db="nucleotide", id=ID_sekwencji, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")

except Exception as error:
    print(f"Nie udało się pobrać sekwencji {error}")
    exit()

print(f"\nID: {record.id}" )
print(f"\nSekwencja od użytkownika:\n{record.seq[:200]}")


# Zapisanie podanej sekwencji do pliku
with open("sekwencje.txt", "a") as plik:
    SeqIO.write(record, plik, "fasta")


# Funkcja do zliczenia nukleotydów
def liczenie_nukleotydow(sekwencja):
    liczenie = {
        "A": sekwencja.count("A"),
        "C": sekwencja.count("C"),
        "G": sekwencja.count("G"),
        "T": sekwencja.count("T")
    }
    return liczenie


# Wywołanie funkcji dla każdej sekwencji z pliku
for i, seq in enumerate(lista_sekwencji):
    print(f"\nPrzetwarzam sekwencję: {i+1}: {seq}")
    
    policzone_nukleotydy = liczenie_nukleotydow(seq)
    print("\nZliczone nukleotydy: ")
    
    for nukleotyd, liczba in policzone_nukleotydy.items():
        print(f"{nukleotyd}: {liczba}")


# Zliczanie par GC     
    ilosc_par_gc = seq.count("GC") 
    print("\nIlosć par GC: ", ilosc_par_gc)
    
# Pierwsze wystąpienie pary GC
lokalizacja_sekwencji_gc = seq.find("GC")   
print("\nIndeks pierwszego wystąpienia pary GC: ", lokalizacja_sekwencji_gc)


# Wyswietlenie informacji o sekwencjach z pliku
for i, seq_records in enumerate(SeqIO.parse("sekwencje.txt", "fasta")):
    print("")
    print(f"Numer sekwencji: {i+1}")
    print("\nID:", seq_records.id)
    print("\nSekwencja:", seq_records.seq)
    print("\nDługosć sekwencji:", len(seq_records))

# Tworzenie słownika z sekwencjami z pliku, aby uzyskać ID sekwencji z pliku
slownik_sekwencji = {}

for seq_records in SeqIO.parse("sekwencje.txt", "fasta"):
    slownik_sekwencji[seq_records.id] = seq_records.seq
    print(slownik_sekwencji.keys())
    print(len(slownik_sekwencji))


# Pobranie sekwencji referencyjnych z NCBI
Entrez.email = input("Podaj swojego maila: ").strip()

ids = slownik_sekwencji.keys()

handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
records = SeqIO.parse(handle, "fasta")

# Zapisanie sekwencji do pliku FASTA
with open("sekwencje_referencyjne_ncbi.txt", "w") as plik:
    SeqIO.write(records, plik, "fasta")


# Usuwanie duplikatów
plik_wyjsciowy_sekwencje = SeqIO.parse("sekwencje.txt", "fasta")
zbior_sekwencji = set()
lista_sekwencji_bez_duplikatow = []

for record in plik_wyjsciowy_sekwencje:
    if str(record.seq) not in zbior_sekwencji:
        lista_sekwencji_bez_duplikatow.append(record)
        zbior_sekwencji.add(str(record.seq))


# Porównanie sekwencji z pliku z sekwencjami referencyjnymi
plik_sekwencje_referencyjne = SeqIO.parse("sekwencje_referencyjne_ncbi.txt", "fasta")
zbior_sekwencje_referencyjne = set(str(record.seq) for record in plik_sekwencje_referencyjne)

lista_poprawne_sekwencje = []
lista_zmodyfikowane_sekwencje = []

for record in lista_sekwencji_bez_duplikatow:
    if str(record.seq) in zbior_sekwencje_referencyjne:
        lista_poprawne_sekwencje.append(record)
    else:
        lista_zmodyfikowane_sekwencje.append(record)

# Zapisanie wyników do pliku
with open("sekwencje_czyste.fasta", "w") as plik:
    SeqIO.write(lista_poprawne_sekwencje, plik, "fasta")


# Funkcja do policzenia par GC
def liczenie_par_gc(sekwencja):
    return sekwencja.count("GC")  


# Utworzenie słownika i DataFrame
import pandas as pd

# Utworzenie danych do DataFrame
sekwencje = []
ids = []
dlugosci = []
motyw_gc = []
nukleotydy = []

for seq in SeqIO.parse("sekwencje_czyste.fasta", "fasta"):
    seq_str = str(seq.seq)
    
    sekwencje.append(seq_str)
    ids.append(seq.id)
    dlugosci.append(len(seq_str))
    motyw_gc.append(liczenie_par_gc(seq_str))
    nukleotydy.append(liczenie_nukleotydow(seq_str))
    
df = pd.DataFrame({
    "sekwencja": sekwencje,
    "ID": ids,
    "dlugosc": dlugosci,
    "motyw_gc": motyw_gc,
    "zliczone_nt": nukleotydy
})

print(df)

# Histogram: wykres przedstawiający długosć sekwencji
import matplotlib.pyplot as plt
import seaborn as sns

# Wykres słupkowy: wizualizacja wyników zliczania wysp GC
plt.figure(figsize=(10,6))
sns.boxplot(x=df["ID"], y=df["motyw_gc"])
plt.xlabel("ID sekwencji")
plt.ylabel("Liczba wystąpień GC")
plt.title("Wykres ilosci wystąpień GC")
plt.xticks(rotation=90)
plt.show()
