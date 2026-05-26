# Bakalauro darbas

Monte Karlo simuliacijomis pagrįsta juridinių asmenų nemokumo tikimybės (PD) pasikliovimo intervalų analizė esant nepriklausomiems ir koreliuotiems įsipareigojimų nevykdymo įvykiams.

## Darbo tema

Šiame darbe nagrinėjamos juridinių asmenų nemokumo tikimybės (angl. *probability of default*, PD) įverčių pasikliovimo intervalų savybės.

Analizuojami klasikiniai binominės proporcijos pasikliovimo intervalai:

- Wald;
- Clopper-Pearson;
- Wilson;
- Agresti-Coull.

Tyrime lyginami du atvejai:

- nepriklausomi įsipareigojimų nevykdymo įvykiai;
- koreliuoti įsipareigojimų nevykdymo įvykiai, generuojami taikant vienfaktorinį Vašiček modelį.

Papildomai nagrinėjama efektyvaus imties dydžio korekcija, skirta įvertinti priklausomybės tarp skolininkų poveikį pasikliovimo intervalų savybėms.

## Repozitorijos struktūra
```text
.
├── Bakalauras.R
├── data/
├── plots/
├── README.md
└── LICENSE
```text


## Failų aprašymas

### `Bakalauras.R`

Pagrindinis R programos failas, kuriame pateikiamas visas tyrimo kodas:

- nepriklausomų įsipareigojimų nevykdymo įvykių generavimas;
- koreliuotų įsipareigojimų nevykdymo įvykių generavimas taikant vienfaktorinį Vašiček modelį;
- Wald, Clopper--Pearson, Wilson ir Agresti--Coull pasikliovimo intervalų skaičiavimas;
- Monte Karlo simuliacijų vykdymas;
- empirinės padengimo tikimybės skaičiavimas;
- vidutinio intervalo pločio skaičiavimas;
- stebimo nemokumo dažnio standartinio nuokrypio skaičiavimas;
- apatinės intervalo ribos, lygios nuliui, dažnio skaičiavimas;
- efektyvaus imties dydžio korekcijos taikymas;
- rezultatų lentelių generavimas;
- grafikų generavimas.

### `data/`

Aplanke pateikiami CSV formato rezultatų failai, sugeneruoti vykdant Monte Karlo simuliacijas.

Šie failai naudojami:

- pagrindinėms darbo lentelėms sudaryti;
- prieduose pateikiamoms pilnoms rezultatų lentelėms sudaryti;
- nepriklausomo ir koreliuoto atvejų palyginimui;
- efektyvaus imties dydžio korekcijos rezultatų analizei.

### `plots/`

Aplanke pateikiami tyrime naudoti grafikai.

Grafikai naudojami praktinėje darbo dalyje vizualiai palyginti:

- pasikliovimo intervalų empirinę padengimo tikimybę;
- vidutinį pasikliovimo intervalų plotį;
- stebimo nemokumo dažnio standartinį nuokrypį;
- nepriklausomo ir koreliuoto modelių skirtumus;
- efektyvaus imties dydžio korekcijos poveikį.


## Naudojama programinė įranga

Analizė atlikta naudojant R programavimo kalbą.

Pagrindiniai naudoti R paketai:

- `dplyr`;
- `tidyr`;
- `ggplot2`;
- `patchwork`;
- `purrr`.

## Pastaba

Darbe naudojami sintetiniai duomenys, sugeneruoti Monte Karlo simuliacijomis. Realūs bankų kredito portfelių duomenys darbe nėra naudojami.
