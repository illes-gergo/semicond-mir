# Sweep for τ variable

include("main.jl")

taus = 25e-15:25e-15:200e-15

@threads for tau in taus 
  runcalc(tau)
end
