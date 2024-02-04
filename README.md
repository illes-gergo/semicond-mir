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
```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.10.0 (2023-12-25)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> include("reader.jl");

julia> sim = readData("DB-500fs-2um");

julia> sim.exportEffic()

julia> sim.plotEffic()
[ Info: Listening on: 127.0.0.1:9321, thread id: 2

julia> sim.plotTHzField(1.3)
1.3 mm kristályhossz kiválasztva

julia> sim.plotTHzSpect(1.3)
1.3 mm kristályhossz kiválasztva
```
![Képernyőkép 2024-02-04 121847](https://github.com/illes-gergo/semicond-mir/assets/79720047/5dae1d19-7faf-48aa-99fa-5f9117fc8d57)
