# Model-til-forudsigelse-af-svovbrintekorrison-af-m-rtel

Denne model er delt op i fire årstider, hvor vandføring og temperatur varierer. I "Korrosionsrate samlet" kombineres de fire scripts.
For at bruge dem skal hver enkelt script for årstiderne køres uden at workspace slettes, og til sidst køres "korrosionsrate samlet".
Modellen er tilpasset en gravtationsledning af beton placeret lige før Ikast Renseanlæg.
Modellen beregner sulfid- og svovlbrintekoncentrationer i gas- og vandfasen, og herudfra beregnes en korrosionsrate.

Excel-arket "Opholdstider" udregner opholdstider af spildevand i en trykledning, der er tilkoblet gravitationsledningen, givet ud fra en vandmængde og en antaget døgnvariation. 
Ønskes der at ændre i mængden af spildevand der tilledes trykledningen, skal det derfor gøres i excel-arket. Det samme gælder oplysninger om f.eks. pumpekapacitet, pumpesump osv.

Der er ikke inkluderet kode til figurer.
