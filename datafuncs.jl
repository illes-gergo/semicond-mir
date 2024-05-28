include("typedefs.jl")

function printInputs2Console(inp::userInputs)
  outstring = "Input Parameters:\n"
  outstring *= "Pump Wavelength: $(inp.lambda0*1e6) μm\n"
  outstring *= "Half-Width Half-Maximum: $(inp.tau*1e12) ps\n"
  outstring *= "Peak Intensity: $(inp.I0*1e-13) GW/cm^2\n"
  outstring *= "THz Matched Frequency: $(inp.nu0*1e-12) THz\n"
  outstring *= "Crystal Length: $(inp.z_end*1e3) mm\n"
  if inp.cry == 3
    material = "GaP"
  elseif inp.cry == 4
    material = "GaAs"
  elseif inp.cry == 0
    material = "LN"
  else
    error("Ismeretlen kristály anyag")
  end
  outstring *= "Crystal: $(material)\n"
  outstring *= "Temperature: $(inp.T) K\n\n"
  
  outstring *= "Run Parameters:\n"
  outstring *= "Többfotonos abszorpció rendje: $(inp.MPAorder)\n"
include("typedefs.jl")

function printInputs2Console(inp::userInputs)
  outstring = "Input Parameters:\n"
  outstring *= "Pump Wavelength: $(inp.lambda0*1e6) μm\n"
  outstring *= "Half-Width Half-Maximum: $(inp.tau*1e12) ps\n"
  outstring *= "Peak Intensity: $(inp.I0*1e-13) GW/cm^2\n"
  outstring *= "THz Matched Frequency: $(inp.nu0*1e-12) THz\n"
  outstring *= "Crystal Length: $(inp.z_end*1e3) mm\n"
  if inp.cry == 3
    material = "GaP"
  elseif inp.cry == 4
    material = "GaAs"
  elseif inp.cry == 0
    material = "LN"
  else
    error("Ismeretlen kristály anyag")
  end
  outstring *= "Crystal: $(material)\n"
  outstring *= "Temperature: $(inp.T) K\n\n"
  
  outstring *= "Run Parameters:\n"
  outstring *= "Multi-photon Absorption Order: $(inp.MPAorder)\n"
  outstring *= "Database Filename: \"$(inp.DB_Name)\"\n"
  outstring *= "Spatial Step Length (dz): $(inp.dz)\n"
  outstring *= "Number of Points Recorded: $(inp.N)\n\n"

  outstring*= """Start Calculation typing "runcalc()"\n"""
  print(outstring)
  return nothing
end
  outstring *= "Adatbázis neve: \"$(inp.DB_Name)\"\n"
  outstring *= "Térbeli lépésköz: $(inp.dz)\n"
  outstring *= "Felvett adatpontok száma: $(inp.N)\n\n"

  outstring*= "Számolás indítása a runcalc() függvényhívással"
  print(outstring)
  return nothing
end
