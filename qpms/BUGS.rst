gaunt.c
=======
abort při určitých vstupech, např. -6 10 8 10
(ačkoliv fortran originál paří bez problémů)


vswf.c
======
Rozklad a znovusložení rovinné vlny vrací nesprávné výsledky pro radiální část podélných
vln, viz test_planewave_decomposition.c. Při reálné amplitudě (E) rovnoběžné s vlnovým
vektorem (k) to vypadá, že reálná část „radiální“ (ve vztahu k počátku souřadnic) složky 
je správně, ale imaginární část neodpovídá. Pakliže amplituda E je komplexní, chyba je 
i v reálné části (pro čistě imaginární amplitudu je chyba *právě* v reálné části).

Chybu jsem hledal především ve výpočtu kulových vln v bodě, tj. ve funkci
qpms_vswf_fill(), ale tam se mi ji zatím najít nepodařilo. Pochopitelně
také může být blbě přímo rozklad, tj. qpms_planewave2vswf_fill_sph().
Výše uvedené (chyba jen v imaginární části pro reál. amplitudu) snad časem napoví.

Neradiální část (podélných vln) je snad v pořádku. Rozklad čistě příčných vln (k.E = 0) 
funguje dobře, což pro elektromagnetismus přirozeně stačí, ale stejně....

Alternativní funkce qpms_vswf_fill_alternative() *nevrací* stejné výsledky
jako qpms_vswf_fill, oproti qpms_vswf_fill() je nejspíš ještě chybovější.
