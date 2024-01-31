# Setting Parameters
include("typedefs.jl")

function setInput()::userInputs
  lambda0 = 1.0e-6
  tau = 0.5e-12
  I0 = 80e13
  dz = 5e-6
  nu0 = 0.5e12
  z_end = 10e-3
  DB_name = "DB"
  differential_equation = 3 # 4 = multiphoton, 2 = shg, 1 = n2
  cry = 3
  T = 300

  return userInputs(lambda0=lambda0, tau=tau, I0=I0, dz=dz, nu0=nu0, z_end=z_end, DB_Name=DB_name, DifferentialEquationSum=differential_equation, cry=cry, T=T)
end
