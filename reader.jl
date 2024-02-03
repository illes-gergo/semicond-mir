using HDF5, PlotlyJS, DelimitedFiles

include("plotLayouts.jl")

function readData(s::String)
  FID = h5open(s)
  return SimData(file=FID)
end

function printInputs2Console(file)
  outstring = "Bemeneti paraméterek:\n"
  outstring *= "Központi hullámhossz: $(read(file["inp/lambda0"])*1e6) μm\n"
  outstring *= "Intenzitás félértékszélesség: $(read(file["inp/tau"])*1e12) ps\n"
  outstring *= "Csúcsintenzitás: $(read(file["inp/I0"])*1e-13) GW/cm^2\n"
  outstring *= "Sebességillesztési frekvencia: $(read(file["inp/nu0"])*1e-12) THz\n"
  outstring *= "Sebességillesztés szöge: $(rad2deg(read(file["gamma"])))\n"
  outstring *= "Kristály teljes hossza: $(read(file["inp/z_end"])*1e3) mm\n"
  matchoice = read(file["/inp/cry"])
  if matchoice == 3
    material = "GaP"
  elseif matchoice == 4
    material = "GaAs"
  else
    error("Ismeretlen kristály anyag")
  end
  outstring *= "Kristály anyaga: $(material)\n"
  outstring *= "Hőmérséklet: $(read(file["inp/T"])) K\n\n"

  outstring *= "Futtatási paraméterek:\n"
  outstring *= "Adatbázis eredeti neve: \"$(read(file["inp/DB_Name"]))\"\n"
  outstring *= "Térbeli lépésköz: $(read(file["inp/dz"])*1e6) μm\n"
  outstring *= "Felvett adatpontok száma: $(read(file["inp/N"])) darab\n\n"
  print(outstring)
  return nothing
end

function plotEfficiency(file)
  z = read(file["z"]) * 1e3
  effic = read(file["effic"]) * 100
  return plot(scatter(x=z, y=effic), merge(lout_general, lout_effic))
end

function plotFreeCarriers(file)
  z = read(file["z"]) * 1e3
  Nc = read(file["Nc"])
  return plot(scatter(x=z, y=Nc), merge(lout_general, lout_freecarr))
end

function findnearest(zsave::Vector, value::Number)
  workArray = abs.(zsave .- value)
  _, index = findmin(workArray)
  println("$(zsave[index]*1e3) mm kristályhossz kiválasztva")
  return index
end

function pltPmpSpct(zval, file)
  idx = findnearest(read(file["zsave"]), zval)
  nu = read(file["nu"]) * 1e-12
  spct = read(file["$(idx)/Aop"])
  return plot(scatter(x=nu, y=spct), merge(lout_general, lout_spectr))
end

function pltTHzSpct(zval, file)
  idx = findnearest(read(file["zsave"]), zval)
  nu = read(file["nu"]) * 1e-12
  spct = read(file["$(idx)/ATHz"])
  return plot(scatter(x=nu, y=spct), merge(lout_general, lout_spectr))
end
function pltTHzField(zval, file)
  idx = findnearest(read(file["zsave"]), zval)
  t = read(file["t"]) * 1e12
  field = read(file["$(idx)/ETHz"]) * 1e-5
  return plot(scatter(x=t, y=field), merge(lout_general, lout_thz))
end

function pltPmpInt(zval, file)
  idx = findnearest(read(file["zsave"]), zval)
  t = read(file["t"]) * 1e12
  field = read(file["$(idx)/Eop"]) * 1e-13
  return plot(scatter(x=t, y=field), merge(lout_general, lout_pmp))
end

function checkpath(foldName::String)
  if !ispath(foldName)
    mkdir(foldName)
  end
  return nothing
end
function printInputs2File(file)
  DB_Name = read(file["inp/DB_Name"]) * "-exports"
  checkpath(DB_Name)

  outstring = "Bemeneti paraméterek:\n"
  outstring *= "Központi hullámhossz: $(read(file["inp/lambda0"])*1e6) μm\n"
  outstring *= "Intenzitás félértékszélesség: $(read(file["inp/tau"])*1e12) ps\n"
  outstring *= "Csúcsintenzitás: $(read(file["inp/I0"])*1e-13) GW/cm^2\n"
  outstring *= "Sebességillesztési frekvencia: $(read(file["inp/nu0"])*1e-12) THz\n"
  outstring *= "Sebességillesztés szöge: $(rad2deg(read(file["gamma"])))\n"
  outstring *= "Kristály teljes hossza: $(read(file["inp/z_end"])*1e3) mm\n"
  matchoice = read(file["/inp/cry"])
  if matchoice == 3
    material = "GaP"
  elseif matchoice == 4
    material = "GaAs"
  else
    error("Ismeretlen kristály anyag")
  end
  outstring *= "Kristály anyaga: $(material)\n"
  outstring *= "Hőmérséklet: $(read(file["inp/T"])) K\n\n"

  outstring *= "Futtatási paraméterek:\n"
  outstring *= "Adatbázis eredeti neve: \"$(read(file["inp/DB_Name"]))\"\n"
  outstring *= "Térbeli lépésköz: $(read(file["inp/dz"])*1e6) μm\n"
  outstring *= "Felvett adatpontok száma: $(read(file["inp/N"])) darab\n\n"
  write(DB_Name * "/inputs.txt", outstring)
  return nothing
end

function exportEfficiency(file)
  z = read(file["z"]) * 1e3
  effic = read(file["effic"]) * 100
  DB_Name = read(file["inp/DB_Name"]) * "-exports"
  checkpath(DB_Name)
  writedlm(DB_Name * "/effic.txt", [z;; effic])
  return nothing
end

function exportCarriers(file)
  z = read(file["z"]) * 1e3
  Nc = read(file["Nc"])
  DB_Name = read(file["inp/DB_Name"]) * "-exports"
  checkpath(DB_Name)
  writedlm(DB_Name * "/Nc.txt", [z;; Nc])
  return nothing
end

function Dummy(file)
  println(read(file["inp/lambda0"]))
  println("Dummy function ran")
end

function exportPmpSpct(zval, file)
  idx = findnearest(read(file["zsave"]), zval)
  nu = read(file["nu"]) * 1e-12
  spct = read(file["$(idx)/Aop"])
  DB_Name = read(file["DB_Name"])
  checkpath(DB_Name)
  zselect = read(file["zsave"])[idx] * 1e3
  writedlm(DB_Name * "/pump_spectrum-$(zselect)mm.txt".[nu;; spct])
  return nothing
end

function exportTHzSpct(zval, file)
  idx = findnearest(read(file["zsave"]), zval)
  nu = read(file["nu"]) * 1e-12
  spct = read(file["$(idx)/ATHz"])
  DB_Name = read(file["DB_Name"])
  checkpath(DB_Name)
  zselect = read(file["zsave"])[idx] * 1e3
  writedlm(DB_Name * "/thz_spectrum-$(zselect)mm.txt".[nu;; spct])
  return nothing
end
function exportTHzField(zval, file)
  idx = findnearest(read(file["zsave"]), zval)
  t = read(file["t"]) * 1e12
  field = read(file["$(idx)/ETHz"]) * 1e-5
  DB_Name = read(file["DB_Name"])
  checkpath(DB_Name)
  zselect = read(file["zsave"])[idx] * 1e3
  writedlm(DB_Name * "/thz_field-$(zselect)mm.txt".[t;; field])
  return nothing
end

function exportPmpInt(zval, file)
  idx = findnearest(read(file["zsave"]), zval)
  t = read(file["t"]) * 1e12
  field = read(file["$(idx)/Eop"]) * 1e-13
  DB_Name = read(file["DB_Name"])
  checkpath(DB_Name)
  zselect = read(file["zsave"])[idx] * 1e3
  writedlm(DB_Name * "/pump_intensity-$(zselect)mm.txt".[t;; field])
  return nothing
end

@kwdef struct SimData
  file
  printUserInputs::Function = () -> printInputs2Console(file)
  plotEffic::Function = () -> plotEfficiency(file)
  plotNc::Function = () -> plotFreeCarriers(file)
  plotPumpSpect::Function = (z) -> pltPmpSpct(z * 1e-3, file)
  plotPumpInt::Function = (z) -> pltPmpInt(z * 1e-3, file)
  plotTHzSpect::Function = (z) -> pltTHzSpct(z * 1e-3, file)
  plotTHzField::Function = (z) -> pltTHzField(z * 1e-3, file)

  exportEffic::Function = () -> exportEfficiency(file)
  exportNc::Function = () -> exportCarriers(file)
  exportUserInputs::Function = () -> printInputs2File(file)
  exportPumpSpect::Function = (z) -> exportPmpSpct(z * 1e-3, file)
  exportPumpInt::Function = (z) -> exportPmpInt(z * 1e-3, file)
  exportTHzSpect::Function = (z) -> exportTHzSpct(z * 1e-3, file)
  exportTHzField::Function = (z) -> exportTHzField(z * 1e-3, file)
end
