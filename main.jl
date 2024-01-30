using Symbolics
using FFTW
using DSP
using Base.Threads
using HDF5

include("fuggvenyek.jl")
include("typedefs.jl")
include("input.jl")

function runcalc()
  inputs = setInput()
  natConsts = naturalConstants()
  cry = inputs.cry
  T = inputs.T
  c = natConsts.c0
  lambda0 = inputs.lambda0
  N = 4e4
  tau = inputs.tau
  I0 = inputs.I0
  khi_eff = 2 * deffTHz(cry)
  e0 = natConsts.e0
  nu0 = 0.5e12

  dz = inputs.dz
  z_vegso = inputs.z_end
  z = 0:dz:z_vegso
  omega0 = 2 * pi * c / lambda0

  elochirp = 0 * 1 * z_vegso / 2

  omegaMAX = 5 * 2 * omega0

  dt = 2 * pi / omegaMAX
  t = (0:N-1) * dt
  t = t .- t[end] / 2

  domega = omegaMAX / N
  dnu = domega / 2 / pi

  omega = (0:N-1) * domega
  nu = omega / 2 / pi

  deltaOmega = 2 * sqrt(2 * log(2)) / tau

  lambda = 2 * pi * c ./ omega
  lambda[1] = lambda[2]
  ngp0 = ngp(lambda0, T, cry)
  np0 = neo(lambda0, T, cry)
  ngpSH = ngp(lambda0 / 2, T, cry)
  npSH = neo(lambda0 / 2, T, cry)

  nTHz = nTHzo(2 * pi * nu0, T, cry)
  vfTHz = c ./ nTHz

  gamma = acos.(ngp0 / nTHz)

  A0t = sqrt(2 * I0 / neo(lambda0, T, cry) / e0 / c)

  n_omega = neo(lambda, T, cry)
  k_OMEGA = real(omega .* nTHzo(omega, T, cry) / c)
  k_OMEGA0 = real(omega .* nTHzo(2 * pi * nu0, T, cry) / c)

  ddk_omega = -ngp0 .^ 2 / omega0 / c / np0 * tan(gamma)^2
  k_omega = real(1 / cos(gamma) .* (omega .* n_omega / c + 1 * (omega .- omega0) .^ 2 / 2 .* ddk_omega))
  # + 1e5;
  ddk_omegaSH = -ngpSH .^ 2 / omega0 / 2 / c / npSH * tan(gamma)^2
  k_omegaSH = real(1 / cos(gamma) .* (omega .* n_omega / c + 1 * (omega .- 2 * omega0) .^ 2 / 2 .* ddk_omegaSH))
  # + 1e5;
  k_omega0 = real.(1 ./ cos(gamma) .* (omega .* ngp0 / c))
  k_omegaSH0 = real.(1 ./ cos(gamma) .* (omega .* ngpSH / c))

  A0 = sqrt(2 * I0 / neo(lambda0, T, cry) / e0 / c) * tau / (2 * sqrt(2 * pi * log(2)))

  Aot = A0t * exp.(-2 * log(2) * t .^ 2 / tau^2) .* exp.(1im * omega0 * t)
  Aop = fft(Aot) .* exp.(1im * (k_omega - k_omega0) * elochirp) * dt / 2 / pi
  Aop0 = copy(Aop)

  ATHz = zeros(size(Aop))
  ASH = zeros(size(Aop))

  pF = sum(abs.(Aop) .^ 2) * np0

  FI = 4 * nTHz .^ 2 ./ (1 + nTHz) .^ 2
  FA = 2 * nTHz ./ (1 + nTHz)

  A_komp = differentialEqInputs(ATHz, Aop0, ASH)

  effic = zeros(size(z))
  efficSH = zeros(size(effic))
  n2 = n2value(cry)
  RTC = runTimeConstants(gamma=gamma, khi_eff=khi_eff, deff=deff(cry), n2=n2,
    omega=omega,
    omega0=omega0, domega=domega, k_omega=k_omega,
    k_OMEGA=k_OMEGA, k_omegaSH=k_omegaSH, dnu=dnu,
    absorption=aTHzo(omega, T, cry))

  misc = miscInput(NC=natConsts, IN=inputs, RTC=RTC)
  #FID = h5open(inputs.DB_Name, "w")
  for ii in eachindex(z)[2:end]
    A_komp = RK4_M(diffegy_conv, dz, z[ii-1], A_komp, misc)
    ATHz = A_komp.ATHz
    Aop = A_komp.Aop
    ASH = A_komp.ASH
    effic[ii] = sum(abs.(ATHz) .^ 2 .* FI) / pF
    efficSH[ii] = sum(abs.(ASH) .^ 2 .* npSH) / pF

    Iop = np0 * e0 * c / 2 * abs.((ifft(Aop .* exp.(-1im * (k_omega - k_omega0) * z[ii]))) * omegaMAX) .^ 2
    ETHz = real.((ifft(FA .* ATHz .* exp.(-1im * (k_OMEGA - k_OMEGA0) * z[ii])))) * omegaMAX
    ISH = abs.((ifft(ASH .* exp.(-1im * (k_omegaSH - k_omegaSH0) * z[ii]))) * omegaMAX) .^ 2
    #DataBaseWriter(FID, z[ii], Aop, Iop, ATHz, ETHz, ASH, ISH)
    println(ii)
  end
  #DataBaseEnder(FID, z, t, nu, effic, efficSH)
  #close(FID)

  #display(plot(z, effic))
  println("Végeztem!")
  return nothing
end
