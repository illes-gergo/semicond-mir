import Base.+, Base.*, Base./

@kwdef struct userInputs
  lambda0::Float64
  tau::Float64
  I0::Float64
  dz::Float64
  nu0::Float64
  z_end::Float64
  DB_Name::String
  DifferentialEquationSum::Int
  cry::Int
  T::Float64
end

@kwdef struct differentialEqInputs
  ATHz::Vector{ComplexF64}
  Aop::Vector{ComplexF64}
  ASH::Vector{ComplexF64}
end

@kwdef struct naturalConstants
  c0::Float64 = 3e8
  e0::Float64 = 8.8541878128e-12
end

@kwdef struct runTimeConstants
  gamma::Float64
  khi_eff::Float64
  deff::Float64
  n2::Float64
  omega::Vector{Float64}
  omega0::Float64
  domega::Float64
  k_omega::Vector{Float64}
  k_OMEGA::Vector{Float64}
  k_omegaSH::Vector{Float64}
  dnu::Float64
  absorption::Vector{Float64}
end

@kwdef struct miscInput
  NC::naturalConstants
  IN::userInputs
  RTC::runTimeConstants
end

# ERROR: MethodError: no method matching *(::differentialEqInputs, ::Float64)

function *(a::differentialEqInputs, b::Float64)
  return (differentialEqInputs(ATHz=a.ATHz .* b,
    Aop=a.Aop .* b, ASH=a.ASH .* b))
end

function *(b::Float64, a::differentialEqInputs)
  return (differentialEqInputs(ATHz=a.ATHz .* b,
    Aop=a.Aop .* b, ASH=a.ASH .* b))
end
# ERROR: MethodError: no method matching /(::differentialEqInputs, ::Int64)

function /(a::differentialEqInputs, b::Int64)
  return (differentialEqInputs(ATHz=a.ATHz ./ b,
    Aop=a.Aop ./ b, ASH=a.ASH ./ b))
end

#  MethodError: no method matching +(::differentialEqInputs, ::differentialEqInputs)

function +(a::differentialEqInputs, b::differentialEqInputs)
  return (differentialEqInputs(ATHz=a.ATHz .+ b.ATHz,
    Aop=a.Aop .+ b.Aop, ASH=a.ASH .+ b.ASH))
end

# ERROR: MethodError: no method matching *(::Int64, ::differentialEqInputs)

function *(b::Int64, a::differentialEqInputs)
  return (differentialEqInputs(ATHz=a.ATHz .* b,
    Aop=a.Aop .* b, ASH=a.ASH .* b))
end
