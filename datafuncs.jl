include("typedefs.jl")

function printInputs2Console(inp::userInputs)
  outstring = "Bemeneti paraméterek:\n"
  outstring *= "Központi hullámhossz: $(inp.lambda0*1e6) μm\n"
  outstring *= "Intenzitás félértékszélesség: $(inp.tau*1e12) ps\n"
  outstring *= "Csúcsintenzitás: $(inp.I0*1e-13) GW/cm^2\n"
  outstring *= "Sebességillesztési frekvencia: $(inp.nu0*1e-12) THz\n"
  outstring *= "Kristály teljes hossza: $(inp.z_end*1e3) mm\n"
  if inp.cry == 3
    material = "GaP"
  elseif inp.cry == 4
    material = "GaAs"
  else
    error("Ismeretlen kristály anyag")
  end
  outstring *= "Kristály anyaga: $(material)\n"
  outstring *= "Hőmérséklet: $(inp.T) K\n\n"
  
  outstring *= "Futtatási paraméterek:\n"
  outstring *= "Adatbázis neve: \"$(inp.DB_Name)\"\n"
  outstring *= "Térbeli lépésköz: $(inp.dz)\n"
  outstring *= "Felvett adatpontok száma: $(inp.N)\n\n"

  outstring*= "Számolás indítása a runcalc() függvényhívással"
  print(outstring)
  return nothing
end
