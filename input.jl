# Setting Parameters

#= function setInput()::userInputs
  lambda0 = 2.35e-6
  tau = 0.1e-12
  I0 = 20e13
  dz = 5e-6
  nu0 = 0.5e12
  z_end = 2e-3
  DB_name = "DB-20GW-HJ-PARAMS-4"
  differential_equation = 3 # 4 = multiphoton, 2 = shg, 1 = n2
  cry = 3
  T = 300
  N = 4e4
  MPAorder = 4
  return userInputs(lambda0=lambda0, tau=tau, I0=I0, dz=dz, nu0=nu0, z_end=z_end, DB_Name=DB_name, DifferentialEquationSum=differential_equation, cry=cry, T=T, N=N, MPAorder=MPAorder)
end =#

@kwdef struct userInputs
  lambda0::Float64 = 1.75e-6
  tau::Float64 = 25e-15
  I0::Float64 = 17.5e13
  dz::Float64 = 5e-6
  nu0::Float64 = 2.0e12
  z_end::Float64 = 2e-3
  DB_Name::String = "test_dir/luis-proba-2"
  DifferentialEquationSum::Int = 0
  cry::Int = 3
  T::Float64 = 300
  N::Int=2e4
  MPAorder::Int = 4
end
