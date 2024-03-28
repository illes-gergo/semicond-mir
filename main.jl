using Symbolics
using FFTW
using DSP
using Base.Threads
using HDF5

include("fuggvenyek.jl")
include("typedefs.jl")
include("input.jl")
include("datafuncs.jl")

function runcalc()
  inputs = userInputs()
  natConsts = naturalConstants()
  cry = inputs.cry
  T = inputs.T
  c = natConsts.c0
  lambda0 = inputs.lambda0
  N::Int = inputs.N # 4e4 by default
  tau = inputs.tau
  I0 = inputs.I0
  khi_eff = 2 * deffTHz(cry)
  e0 = natConsts.e0
  nu0 = inputs.nu0
  dz = inputs.dz
  z_vegso = inputs.z_end
  z = 0:dz:z_vegso
  saveInterval = 1.0e-4
  zsave = 0:saveInterval:z_vegso
  omega0 = 2 * pi * c / lambda0
  MPAorder = inputs.MPAorder
  #beta = betaN(cry,MPAorder)

  elochirp = 0 * 1 * z_vegso / 2

  omegaMAX = 2.5 * omega0

  dt = 2 * pi / omegaMAX
  t = (0:N-1) * dt
  t = t .- t[end] / 2

  domega = omegaMAX / N
  dnu = domega / 2 / pi

  omega = (0:N-1) * domega
  nu = omega / 2 / pi


  lambda = 2 * pi * c ./ omega
  lambda[1] = lambda[2]
  ngp0 = ngp(lambda0, T, cry)
  np0 = neo(lambda0, T, cry)
  ngpSH = ngp(lambda0 / 2, T, cry)
  npSH = neo(lambda0 / 2, T, cry)

  nTHz = nTHzo(2 * pi * nu0, T, cry)

  gamma = acos.(ngp0 / nTHz)

  A0t = sqrt(2 * I0 / neo(lambda0, T, cry) / e0 / c)

  n_omega = neo(lambda, T, cry)
  k_OMEGA = real(omega .* nTHzo(omega, T, cry) / c)
  k_OMEGA0 = real(omega .* nTHzo(2 * pi * nu0, T, cry) / c)

  ddk_omega = -ngp0 .^ 2 / omega0 / c / np0 * tan(gamma)^2
  k_omega = real(1 / cos(gamma) .* (omega .* n_omega / c + 1 * (omega .- omega0) .^ 2 / 2 .* ddk_omega))
  ddk_omegaSH = -ngpSH .^ 2 / omega0 / 2 / c / npSH * tan(gamma)^2
  k_omegaSH = real(1 / cos(gamma) .* (omega .* n_omega / c + 1 * (omega .- 2 * omega0) .^ 2 / 2 .* ddk_omegaSH))
  k_omega0 = real.(1 ./ cos(gamma) .* (omega .* ngp0 / c))


  Aot = A0t * exp.(-2 * log(2) * t .^ 2 / tau^2) .* exp.(1im * omega0 * t)
  Aop = fft(Aot) .* exp.(1im * (k_omega - k_omega0) * elochirp) * dt / 2 / pi

  ATHz = zeros(size(Aop))

  pF = sum(abs.(Aop) .^ 2) * np0

  FI = 4 * nTHz .^ 2 ./ (1 + nTHz) .^ 2
  FA = 2 * nTHz ./ (1 + nTHz)

  #A_komp::differentialEqInputs = differentialEqInputs(Aop=Aop, ATHz=ATHz, cumulativePhase=zeros(N), Nc=0)
  A_komp::COdifferentialEqInputs = COdifferentialEqInputs(Aop=Aop, ATHz=ATHz, ASH=zeros(size(Aop)))

  effic = zeros(size(z))
  efficSH = zeros(size(effic))
  n2 = n2value(cry)
  RTC = runTimeConstants(gamma=gamma, khi_eff=khi_eff, deff=deff(cry), n2=n2,
    omega=omega,
    omega0=omega0, domega=domega, k_omega=k_omega,
    k_OMEGA=k_OMEGA, k_omegaSH=k_omegaSH, dnu=dnu, pumpRefInd=np0, NN=N, dt=dt, betaN=0)

  misc = miscInput(NC=natConsts, IN=inputs, RTC=RTC)
  NcArray = zeros(length(z))
  FID = h5open(inputs.DB_Name, "w")
  saveCounter::Int = 1
  for ii in eachindex(z)[2:end]
    A_komp = RK4_M(diffegy_conv, dz, z[ii], A_komp, misc)
    effic[ii] = sum(abs.(A_komp.ATHz) .^ 2 .* FI) / pF
    NcArray[ii] = 0
    if zsave[saveCounter] == z[ii-1] || ii == length(z)
      ATHz = A_komp.ATHz
      Aop = A_komp.Aop
      kdz = RTC.k_OMEGA .* z[ii]
      Iop = real.(np0 * e0 * c / 2 * abs.((ifft(Aop .* exp.(-1im * (k_omega - k_omega0) * z[ii]))) * omegaMAX) .^ 2)
      ETHz = real.((ifft(FA .* ATHz .* exp.(-1im * ((-k_OMEGA0) * z[ii] .+ kdz))))) * omegaMAX
      DataBaseWriter(FID, saveCounter, Aop, Iop, ATHz, ETHz)
      saveCounter += 1
    end
    print("$(ii) / $(length(z)) \r")
    flush(stdout)
  end
  DataBaseEnder(FID, z, t, nu, effic, NcArray, inputs, zsave, gamma)
  close(FID)

  #display(plot(z, effic))
  println("Végeztem!")
  return nothing
end

printInputs2Console(userInputs())
