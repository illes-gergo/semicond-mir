# Sweep for τ variable

include("main.jl")

# taus = 25e-15:25e-15:200e-15
taus = [100e-15, 500e-15]

@threads for tau in taus 
  runcalc(userInputs(tau= tau))
end
