import Base.+, Base.*, Base./
include("input.jl")

@kwdef struct differentialEqInputs
  ATHz::Vector{ComplexF64}
  Aop::Vector{ComplexF64}
  # ASH::Vector{ComplexF64}
  cumulativePhase::Vector{Float64}
  Nc::Float64
end

@kwdef struct COdifferentialEqInputs
  ATHz::Vector{ComplexF64}
  Aop::Vector{ComplexF64}
  ASH::Vector{ComplexF64}
end

@kwdef struct LNeqInputs
  ATHz::Vector{ComplexF64}
  Aop::Vector{ComplexF64}
end

@kwdef struct naturalConstants
  c0::Float64 = 3e8
  e0::Float64 = 8.8541878128e-12
  hv::Float64 = 6.626e-34 / 2 / pi
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
  pumpRefInd::Float64
  betaN::Float64
  NN::Int
  dt::Float64
  period::Float64 # Domain period. Each Wafer has length period/2, producing half cycle.
end

@kwdef struct miscInput
  NC::naturalConstants
  IN::userInputs
  RTC::runTimeConstants
end

# ERROR: MethodError: no method matching *(::differentialEqInputs, ::Float64)

function *(a::differentialEqInputs, b::Float64)
  return (differentialEqInputs(ATHz=a.ATHz .* b,
    Aop=a.Aop .* b, Nc=a.Nc .* b, cumulativePhase=a.cumulativePhase .* b))
end

function *(b::Float64, a::differentialEqInputs)
  return (differentialEqInputs(ATHz=a.ATHz .* b,
    Aop=a.Aop .* b, Nc=a.Nc .* b, cumulativePhase=a.cumulativePhase .* b))
end
# ERROR: MethodError: no method matching /(::differentialEqInputs, ::Int64)

function /(a::differentialEqInputs, b::Int64)
  return (differentialEqInputs(ATHz=a.ATHz ./ b,
    Aop=a.Aop ./ b, Nc=a.Nc ./ b, cumulativePhase=a.cumulativePhase ./ b))
end

#  MethodError: no method matching +(::differentialEqInputs, ::differentialEqInputs)

function +(a::differentialEqInputs, b::differentialEqInputs)
  return (differentialEqInputs(ATHz=a.ATHz .+ b.ATHz,
    Aop=a.Aop .+ b.Aop, Nc=a.Nc .+ b.Nc, cumulativePhase=a.cumulativePhase .+ b.cumulativePhase))
end

# ERROR: MethodError: no method matching *(::Int64, ::differentialEqInputs)

function *(b::Int64, a::differentialEqInputs)
  return (differentialEqInputs(ATHz=a.ATHz .* b,
    Aop=a.Aop .* b, Nc=a.Nc .* b, cumulativePhase=a.cumulativePhase .* b))
end

function *(a::COdifferentialEqInputs, b::Number)
  return COdifferentialEqInputs(Aop=a.Aop .* b, ATHz=a.ATHz .* b, ASH=a.ASH .* b)
end

function *(b::Number, a::COdifferentialEqInputs)
  return COdifferentialEqInputs(Aop=a.Aop .* b, ATHz=a.ATHz .* b, ASH=a.ASH .* b)
end

function /(a::COdifferentialEqInputs, b::Number)
  return COdifferentialEqInputs(Aop=a.Aop ./ b, ATHz=a.ATHz ./ b, ASH=a.ASH ./ b)
end

function +(a::COdifferentialEqInputs, b::COdifferentialEqInputs)
  return COdifferentialEqInputs(Aop=a.Aop .+ b.Aop, ATHz = a.ATHz .+ b.ATHz, ASH = a.ASH .+ b.ASH)
end

function *(a::LNeqInputs, b::Number)
  return LNeqInputs(Aop=a.Aop .* b, ATHz=a.ATHz .* b)
end

function *(b::Number, a::LNeqInputs)
  return LNeqInputs(Aop=a.Aop .* b, ATHz=a.ATHz .* b)
end

function /(a::LNeqInputs, b::Number)
  return LNeqInputs(Aop=a.Aop ./ b, ATHz=a.ATHz ./ b)
end

function +(a::LNeqInputs, b::LNeqInputs)
  return LNeqInputs(Aop=a.Aop .+ b.Aop, ATHz = a.ATHz .+ b.ATHz)
end
