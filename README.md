# Semicond-MIR THz-generálást szimuláló szoftver | Semicond-MIR THz generation simulating toolbox
## Telepítés
A program futtatásához szükség van a ```julia``` interpreterre, lehetőség szerint a legfrissebb verziójára valamint a ```git``` verziókövető szoftverre. A repo-t cloneozva nyissunk meg egy terminál ablakot. Indítsuk el a ```julia```-t a következő parancs segítségével:

For running the code it is necessary to install the latest version of ```julia interpreter``` and ```git```version tracker. After cloning the repository start ```julia``` in terminal using the command below:
```
julia --threads=2 --heap-size-hint=4G --project=.
```
Az első alkalommal szükség van a függőségek telepítésére, a következő parancs segítségével:

At first start we have to install the dependencies with the following command:
```
]instantiate
```
Amennyiben lefutott a függőségek telepítése a szoftver használhatóvá válik.

After the process finished everything is ready to go.
## Használat | Usage
### Szimulációk futtatása | Running a simulation
Szimuláció futtatásához nyissuk meg az interpretert:
In order to run a simulation open the terminal with the following command:
```
julia --threads=2 --heap-size-hint=4G --project=.
```
Ezt követően inicializáljuk a szimulációt:
Next the simulation have to be initialized:
```
include("main.jl")
```
A parancs inicializálja a megfelelő függvényeket valamint paramétereket. Ekkor még lehetőség van ezek megváltoztatására az ```input.jl``` fájl tartalmának módosításával. Amennyiben megváltoztattuk az ```input.jl``` fájlt inicializáljunk újra.
Amennyiben mindent rendben találtunk hívjuk meg a ```runcalc``` függvényt üres argumentummal. Ekkor elindul a szimuláció. A program a szimuláció végeztével promptot fog kijelezni. Amennyiben szeretnénk kiíratni a futás idejét a függvényhívás előtt használjuk az ```@elapsed``` makrót.

The last command initializes the run, loads all the functions dependencies and input data for the simulation. At this time it is still possible to change user input data by changing the content of ```input.jl```. If the input changed one has to reinitialize the simulation. If everything is set than one can call the ```runcalc``` function with empty arguments to run the simulation. If you want to see the time the simulation took, use the ```@elapsed``` macro.

### Szimulációk kiértékelése | Evaluating simulations
Első lépésként inicializáljuk a kiértékelő függvényeket, illetve szimulációs objektumot.

The first step is to initialize the evaluating functions and the simulation object.
```
include("reader.jl")
```
A szimulációk kiértékelése során lehetőség van:
1. Adatok kiíratására a terminálra illetve változókba mentése: ```return...``` illetve ```print...``` kezdetű függvények
2. Adatok ábrázolása interaktív grafikonokon: ```plot...``` kezdetű függvények segítségével
3. Adatok exportálása külső programban való használatra: ```export...``` kezdetű függvények.

When evaluating simulations you have several wrapper functions to:
1. Save the raw data to a variable or print it to terminal by using functions starting with ```return...``` and ```print...``` respectfully
2. You can plot the results interactivly by using functions starting with ```plot...```
3. Finally you can export data to ASCII files to process them externaly using the ```export...``` functions.

A függvények használatához szükség van a szimuláció adatait tartalmazó objektum példányosítására:

In order to use the evaluation functions you have to instatiate a simulation object first by doing:
```
sim = readData("DB-name")
```
Ez után a ```sim``` objektummal interaktálva elvégezhetjük az adatok kiértékelését.
Példa:

Then by interacting with the ```sim``` object you can evaluate the run.
Example:
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
