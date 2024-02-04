# Semicond-MIR THz-generálást szimuláló szoftver
## Telepítés
A program futtatásához szükség van a ```julia``` interpreterre, lehetőség szerint a legfrissebb verziójára valamint a ```git``` verziókövető szoftverre. A repo-t cloneozva nyissunk meg egy terminál ablakot. Indítsuk el a ```julia```-t a következő parancs segítségével:
```
julia --threads=2 --heap-size-hint=4G --project=.
```
Az első alkalommal szükség van a függőségek telepítésére, a következő parancs segítségével:
```
]instantiate
```
Amennyiben lefutott a függőségek telepítése a szoftver használhatóvá válik.
## Használat
### Szimulációk futtatása
Szimuláció futtatásához nyissuk meg az interpretert:
```
julia --threads=2 --heap-size-hint=4G --project=.
```
Ezt követően inicializáljuk a szimulációt
```
include("main.jl")
```
A parancs inicializálja a megfelelő függvényeket valamint paramétereket. Ekkor még lehetőség van ezek megváltoztatására az ```input.jl``` fájl tartalmának módosításával. Amennyiben megváltoztattuk az ```input.jl``` fájlt inicializáljunk újra.
Amennyiben mindent rendben találtunk hívjuk meg a ```runcalc``` függvényt üres argumentummal. Ekkor elindul a szimuláció. A program a szimuláció végeztével promptot fog kijelezni. Amennyiben szeretnénk kiíratni a futás idejét a függvényhívás előtt használjuk az ```@elapsed``` makrót.

### Szimulációk kiértékelése
Első lépésként inicializáljuk a kiértékelő függvényeket, illetve szimulációs objektumot.
```
include("reader.jl")
```
A szimulációk kiértékelése során lehetőség van:
1. Adatok kiíratására a terminálra illetve változókba mentése: ```return...``` illetve ```print...``` kezdetű függvények
2. Adatok ábrázolása interaktív grafikonokon: ```plot...``` kezdetű függvények segítségével
3. Adatok exportálása külső programban való használatra: ```export...``` kezdetű függvények.

A függvények használatához szükség van a szimuláció adatait tartalmazó objektum példányosítására:
```
sim = readData("DB-name")
```
Ez után a ```sim``` objektummal interaktálva elvégezhetjük az adatok kiértékelését.
Példa:
