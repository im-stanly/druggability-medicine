# Fpocket

[Manual](https://github.com/Discngine/fpocket/blob/master/doc/MANUAL.md)

[Getting started](https://github.com/Discngine/fpocket/tree/master?tab=readme-ov-file#getting-started)

## Zasada działania
Metoda opiera się na teselacji Voronoia oraz koncepcji sfer alpha.

**Sfera alpha** to geometryczna kula stykająca się z czterema atomami na swojej granicy ale nie zawiera wewnątrz żadnego innego atomu. Promienie sfery określają nam krzywiznę - bardzo małe sfery znajdują się wewnątrz struktury białka, bardzo duże zaś - na zewnątrz. Szukamy sfer o relatywnie średnim promieniu, które odpowiadają kieszeniom na powierzchni.

Działanie algorytmu można podzielić na trzy kroki
1. Wyznaczanie i filtirowanie sfer alpha - za pomocą pakietu QHull przeprowadzana jest tesalacja Voronoia na strukturze białka, aby wtznaczyć współrzędne wierzchołków Voronoia. Te właśnie wierzchołki środkami sfer alpha. Następnie odrzucamy sfery które są za małe i za duże (zdefiniowany parametr minimalny i maksymalny). Otrzymane sfery alfa są dodatkowo etykietowane jako polarne lub apolarne, na podstawie elektroujemności atomów, z którymi się stykają
2. Grupowanie sfer alfa (Clustering): Odfiltrowane sfery są łączone w klastry za pomocą trzystopniowego procesu. Najpierw wykonywana jest zgrubna segmentacja oceniająca odległości między połączonymi wierzchołkami. Następnie odrzuca się pojedyncze sfery na powierzchni białka i łączy małe klastry w większe grupy za pomocą obliczonych środków ciężkości. Na koniec przeprowadzane jest precyzyjne grupowanie metodą wielokrotnych powiązań (multiple linkage clustering), które weryfikuje wzajemne odległości sfer. Klastry niespełniające kryteriów, np. zbyt małe, są z reguły odrzucane
3. Charakteryzacja i ocena (Scoring): Wykryte w ten sposób kieszenie są charakteryzowane w celu przewidzenia ich potencjalnej zdolności do wiązania małych cząsteczek. Program stosuje prostą funkcję punktacyjną bazującą na regresji PLS (Partial Least Squares). Ocena tej zdolności opiera się na pięciu parametrach wyliczanych dla danej kieszeni: znormalizowanej liczbie sfer alfa, znormalizowanej średniej lokalnej gęstości hydrofobowej, odsetku apolarnych sfer alfa, całkowitym wyniku polarności zaangażowanych aminokwasów oraz gęstości sfer alfa (średniej odległości między ich parami)

[Le Guilloux, V., Schmidtke, P. & Tuffery, P. Fpocket: An open source platform for ligand pocket detection. BMC Bioinformatics 10, 168 (2009).](https://doi.org/10.1186/1471-2105-10-168)
## Tech
Wspiera pliki o rozszerzeniu PDB/PDBx/mmcif.
Mamy możliwość oglądania cząsteczek w VMD albo w PyMOL.
Są również wyświetlane sfery alpha, razem z ich środkami. Sfery są w tym samym kolorze jeśli są w tym samym klastrze.
Podobno można również jakoś ustawić żeby wyświetlało czy sfera jest polarna, czy apolarna.


## Ekstrakcja deskryptorów - dpocket
Jest tu moduł do ekstrakcji deskryptorów fizykochemicznych dla kieszeni. Podajemy mu jakieś .pdb i ligand, na tej podstawie wybiera kieszenie z którymi ligand mógłby się związać i wypisuje. Również liczy destryptory dla wszystkich kieszeni, a także dla zestawu atomów, w razie gdybyśmy mieli *explicite* zdefiniowany ligand z białkiem.
[lista deskryptorów](https://github.com/Discngine/fpocket/blob/master/doc/GETTINGSTARTED.md#pocket-descriptors)


## Pliki wyjściowe (z dokumentacji)

###  Output files description

fpocket yields output directly in the directory of the data file, creating a directory using the name of the PDB file followed by the _out extension. Here, the command ll sample/3LKF_out of the current sample run would look something like this:

        total 332
        -rw-r--r-- 1 peter users    769 Nov 29 00:14 3LKF.pml
        -rw-r--r-- 1 peter users    698 Nov 29 00:14 3LKF.tcl
        -rwxr-xr-x 1 peter users     30 Nov 29 00:14 3LKF_PYMOL.sh
        -rwxr-xr-x 1 peter users     41 Nov 29 00:14 3LKF_VMD.sh
        -rw-r--r-- 1 peter users 245835 Nov 29 00:14 3LKF_out.pdb
        -rw-r--r-- 1 peter users   6725 Nov 29 00:14 3LKF_pockets.info
        -rw-r--r-- 1 peter users  49355 Nov 29 00:14 3LKF_pockets.pqr
        -rw-r--r-- 1 peter users   4073 Nov 29 00:14 3LKF_info.txt
        drwxr-xr-x 2 peter users   4096 Nov 29 00:14 pockets

As you can see, fpocket provides a lot of files and another subdirectory. However, majority of these files are necessary for easy visualization of binding pockets. Lets explain the content and utility of each file:

* `3LKF_info.txt`: this file contains human readable information (descriptors) about the pockets found on the protein. Notably this file contains a pocket score (likeliness this is a small molecule binding site) and a druggability score (how druggable the binding site is) Here an extract:
    Pocket 1 :
            Score :         0.490
            Druggability Score :    0.019
            Number of Alpha Spheres :       21
            Total SASA :    19.687
            Polar SASA :    7.611
            Apolar SASA :   12.076
            Volume :        270.934
            Mean local hydrophobic density :        3.000
            Mean alpha sphere radius :      3.816
            Mean alp. sph. solvent access :         0.519
            Apolar alpha sphere proportion :        0.190
            Hydrophobicity score:   23.889
            ...

* `3LKF.pml`: this is a PyMOL script for visualization of binding pockets using PyMOL
* `3LKF.tcl`: this a tcl script for visualization of binding pockets using VMD
* `3LKF_PYMOL.sh`: this is the executable script to launch fast visualization using PYMOL
* `3LKF_VMD.sh`: this is the executable script to launch fast visualization using VMD
* `3LKF_out.pdb`: this is the most important file, it contains the initial PDB structure given as input. Non cofactor HETATM occurrences will be stripped off in this file compared to the original PDB input file. The PDB file contains centers of alpha spheres using the HETATM definition as dummy atoms. These alpha sphere centers are attached in the end of the PDB file, using the STP residue name (for site point). Apolar alpha spheres carry the atom name APOL, polar alpha spheres the atom name POL. Pockets are sets of alpha spheres. They can be distinguished by residue number. Thus residue STP 1 would be the first binding pocket according to fpocket. To show this more clearly here is an extract of the `3LKF_out.pdb`:
        
        ATOM   2349   CD LYS A 299       9.679  16.827 105.636  0.00  0.00           C 0
        ATOM   2350   CE LYS A 299      10.371  16.314 104.370  0.00  0.00           C 0
        ATOM   2351   NZ LYS A 299      11.749  15.794 104.597  0.00  0.00           N 0
        ATOM   2352  OXT LYS A 299       5.240  20.009 107.670  0.58  9.64           O 0
        HETATM    1 APOL STP C   1      27.849  33.435 123.906  0.00  0.00          Ve  
        HETATM    2 APOL STP C   1      29.108  33.195 122.206  0.00  0.00          Ve  
        HETATM    3 APOL STP C   1      28.611  33.141 119.797  0.00  0.00          Ve  
        HETATM    4 APOL STP C   1      26.830  32.143 118.779  0.00  0.00          Ve  

* `3LKF_pockets.pqr`: This file contains all alpha sphere centers, as the 3LKF_out.pdb file, but contains no information about the protein structure. Furthermore using the pqr format enables writing of the van der Waals radius of atoms explicitely in this file. Here this possibility was used to output the radii of alpha spheres of a pocket. Charging this pqr file, one can analyze more precisely the volume recognized by fpocket. Note that, currently only VMD supports reading this format correctly. PyMOL is able to read pqr file, but does not interpret van der Waals radii.    

* `pockets/`: Well, again a subdirectory. But I promise, it's the last one. For development purposes or easy analysis, fpocket proposes this directory which contains according to the current example:

        pocket0_atm.pdb   pocket2_vert.pqr  pocket5_atm.pdb   pocket7_vert.pqr
        pocket0_vert.pqr  pocket3_atm.pdb   pocket5_vert.pqr  pocket8_atm.pdb
        pocket1_atm.pdb   pocket3_vert.pqr  pocket6_atm.pdb   pocket8_vert.pqr
        pocket1_vert.pqr  pocket4_atm.pdb   pocket6_vert.pqr  pocket9_atm.pdb
        pocket2_atm.pdb   pocket4_vert.pqr  pocket7_atm.pdb   pocket9_vert.pqr

* `*_atm.pdb`: These files contain only the atoms contacted by alpha spheres in the given pocket. Complementary to this information, `*_vert.pqr` files contain only the centers and radii of alpha spheres within the respective pocket. As extensions mention, atoms are output in the PDB file format and alpha sphere centers in the PQR file format.


## Opis formatów plików wyjściowych

### 1. Format PDB (Protein Data Bank)
To standardowy format zapisu współrzędnych atomowych. W pliku `_out.pdb` generowanym przez fpocket, atomy białka są zapisane jako `ATOM`, a centra alfa-sfer jako `HETATM`.

#### Struktura kolumnowa PDB (Rekord ATOM/HETATM)

| Kolumny | Pole | Opis | Znaczenie we fpocket |
| :--- | :--- | :--- | :--- |
| **1-6** | `Record name` | `ATOM  ` lub `HETATM` | `HETATM` dla alfa-sfer. |
| **7-11** | `Serial` | Numer porządkowy atomu. | ID punktu. |
| **13-16** | `Atom name` | Nazwa atomu (np. `CA`). | `APOL` (apolarna) lub `POL` (polarna). |
| **18-20** | `ResName` | Nazwa aminokwasu (np. `LYS`). | **`STP`** (Site Point). |
| **22** | `Chain ID` | Identyfikator łańcucha. | Zazwyczaj puste lub `A`. |
| **23-26** | `ResSeq` | Numer aminokwasu. | **Numer kieszeni** (np. 1, 2, 3). |
| **31-38** | `X` | Współrzędna X (Å). | Położenie centrum sfery. |
| **39-46** | `Y` | Współrzędna Y (Å). | Położenie centrum sfery. |
| **47-54** | `Z` | Współrzędna Z (Å). | Położenie centrum sfery. |
| **55-60** | `Occupancy` | Współczynnik zajętości. | Zazwyczaj `0.00` lub `1.00`. |
| **61-66** | `B-factor` | Czynnik temperaturowy. | Zazwyczaj `0.00` (brak info o promieniu). |

> **Przykład w fpocket:**
> `HETATM    1 APOL STP C   1      27.849  33.435 123.906  0.00  0.00`

---

### 2. Format PQR
Format PQR jest rozszerzeniem PDB. Został stworzony, aby obok współrzędnych przechowywać dane o **ładunku** ($q$) i **promieniu** ($r$), co jest niezbędne do obliczeń fizycznych. Fpocket używa go do zapisania wielkości alfa-sfer.

#### Struktura kolumnowa PQR

| Kolumny | Pole | Opis | Znaczenie we fpocket |
| :--- | :--- | :--- | :--- |
| **1-54** | *Jak wyżej* | Tak samo jak w PDB. | Współrzędne i nazewnictwo sfer. |
| **55-62** | **`Charge`** | Ładunek elektryczny ($q$). | We fpocket często powielony **numer kieszeni**. |
| **63-70** | **`Radius`** | Promień atomowy ($r$). | **Realny promień alfa-sfery (w Å)**. |

> **Przykład w fpocket:**
> `ATOM      1 APOL STP     1      27.849  33.435 123.906  1.0000  2.8500`
> *(Tutaj `1.0000` to ID kieszeni, a `2.8500` to promień sfery).*

---

## 3. Bezpośrednie porównanie (Kluczowe różnice)

Najłatwiej zrozumieć różnicę, patrząc na końcówkę linii w obu plikach dla tego samego punktu:

| Format | Końcówka rekordu (od 55. kolumny) | Co to oznacza? |
| :--- | :--- | :--- |
| **PDB** | ` 0.00  0.00` | Brak danych o promieniu. Wizualizacja to tylko mały punkt. |
| **PQR** | ` 1.0000  2.8500` | Wiemy, że to kieszeń nr 1 i ma promień $2.85\ \text{Å}$. |

### Dlaczego to jest ważne przy lokalizacji?
* Jeśli otworzysz **PDB**, zobaczysz "chmurę kropek". Trudno ocenić, gdzie kieszeń się kończy, a gdzie zaczyna.
* Jeśli otworzysz **PQR** w programie obsługującym ten format (np. VMD), program narysuje kule o podanych promieniach. Zobaczysz wtedy **pełny kształt (odlew) kieszeni**, co pozwoli Ci precyzyjnie ocenić jej objętość i granice.

---

### Szybka ściąga dla lokalizacji kieszeni:
1.  **Chcesz listę aminokwasów?** Sprawdź pliki `pocketX_atm.pdb` (format PDB).
2.  **Chcesz środek geometryczny i objętość?** Sprawdź pliki `pocketX_vert.pqr` (format PQR).
3.  **Chcesz parametry liczbowe (score, volume)?** Sprawdź plik `_info.txt`.

Czy potrzebujesz jakiegoś konkretnego skryptu (np. w Pythonie), który wyciągnie te dane z kolumn PQR do tabeli Excela?