include("typedefs.jl")

function deffTHz(cry)
  if cry == 0 # LN
    deff_ = 168e-12
  elseif cry == 2 # ZnTe
    deff_ = 0
  elseif cry == 3 # GaP
    deff_ = 0
  elseif cry == 4 # GaAs
    deff_ = 2 / sqrt(3) * 86.5e-12 #2023-08-04
  elseif cry == 7 # ZnSe
    deff_ = 0
  end
  return deff_
end

function neo(lambda, T, cry)
  if cry == 4 #GaAs Skauli et al. 2003 0.97-17 um
    l = lambda * 1e6
    a0 = 4.372514
    a = [5.466742 0.0242996 1.957522]
    b = [0.4431307 0.8746453 36.9166]

    n = real(sqrt.(Complex.(a0 .+ 1 .+ a[1] * l .^ 2 ./ (l .^ 2 .- b[1]^2) + a[2] * l .^ 2 ./ (l .^ 2 .- b[2]^2) + a[3] * l .^ 2 ./ (l .^ 2 .- b[3]^2))))
  elseif cry == 3
    l = lambda * 1e6
    a1 = 1.39
    a2 = 0.172
    b1 = 4.131
    b2 = 0.234
    c1 = 2.57
    c2 = 0.345
    d1 = 2.056
    d2 = 27.52

    n = @. sqrt(complex(1 + a1 * l .^ 2 ./ (l .^ 2 - a2^2) + b1 * l .^ 2 ./ (l .^ 2 - b2^2) + c1 * l .^ 2 ./ (l .^ 2 - c2^2) + d1 * l .^ 2 ./ (l .^ 2 - d2^2)))
  end
  return n
end

function ngp(lambda, T, cry)
  lambda1 = lambda * 1e6
  if cry == 4
    @variables l

    a0 = 4.372514
    a = [5.466742 0.0242996 1.957522]
    b = [0.4431307 0.8746453 36.9166]

    n0 = real.(sqrt.(a0 + 1 + a[(1)] * l .^ 2 ./ (l .^ 2 - b[(1)]^2) + a[(2)] * l .^ 2 ./ (l .^ 2 - b[(2)]^2) + a[(3)] * l .^ 2 ./ (l .^ 2 - b[(3)]^2)))

    a = n0 - l * Symbolics.derivative(n0, l)
    #l = lambda1;
    ng = Symbolics.value(substitute(a, l => lambda1))
    return ng
  elseif cry == 3
    @variables l
    lambda1 = lambda * 1e6
    a1 = 1.39
    a2 = 0.172
    b1 = 4.131
    b2 = 0.234
    c1 = 2.57
    c2 = 0.345
    d1 = 2.056
    d2 = 27.52

    n0 = sqrt(1 + a1 * l .^ 2 ./ (l .^ 2 - a2^2) + b1 * l .^ 2 ./ (l .^ 2 - b2^2) + c1 * l .^ 2 ./ (l .^ 2 - c2^2) + d1 * l .^ 2 ./ (l .^ 2 - d2^2))


    a = n0 - l * Symbolics.derivative(n0, l)
    #l = lambda1
    ng = Symbolics.value(substitute(a, l => lambda1))
  end
  return ng
end

function nTHzo(omega, T, cry, Nc=0)
  if cry == 4
    nTHz = real.(sqrt.(er(omega, T, cry)))
  elseif cry == 3
    nTHz = real.(sqrt.(complex.(er(omega, T, cry, Nc))))
  end
  return nTHz
end

function diffegy_conv(z, A_kompozit::differentialEqInputs, misc::miscInput)
  ATHz = A_kompozit.ATHz
  Aop = A_kompozit.Aop
  #  ASH = A_kompozit.ASH
  kdz = A_kompozit.cumulativePhase
  NN = misc.RTC.NN


  At = ifft(Aop .* 2 * pi * misc.RTC.dnu * NN)
  It = misc.NC.e0 / 2 * misc.NC.c0 * misc.RTC.pumpRefInd * abs.(ifftshift(ifft(Aop .* exp.(-1im * (misc.RTC.k_omega) * z) * 2 * pi * misc.RTC.dnu * NN))) .^ 2

  Nt = misc.RTC.beta4 .* cumsum(It .^ 4) .* misc.RTC.dt / 4 / misc.NC.hv / misc.RTC.omega0

  ITHzt = abs.(ifftshift(ifft(ATHz .* exp.(-1im * kdz)))) .^ 2
  ITHzt ./= maximum(ITHzt)
  THzint = sum(ITHzt)

  THzint > 0 ? Neff = sum(Nt .* ITHzt) ./ THzint : Neff = 0

  dkdz = real.(misc.RTC.omega .* nTHzo(misc.RTC.omega, misc.IN.T, misc.IN.cry, Neff)) / misc.NC.c0

  n2pm = fft(1im * misc.NC.e0 * misc.RTC.omega0 * misc.RTC.pumpRefInd * misc.RTC.n2 / 2 * abs.(At) .^ 2 .* At) / misc.RTC.dnu / 2 / pi / length(misc.RTC.omega)

  thzAbsorption = 2 * misc.RTC.omega / misc.NC.c0 .* imag.(sqrt.(complex.(er(misc.RTC.omega, misc.IN.T, misc.IN.cry, Neff))))

  t1 = @spawn begin
    temp11 = conv(reverse(conj(Aop) .* exp.(1im .* misc.RTC.k_omega .* z)), (Aop .* exp.(-1im * misc.RTC.k_omega .* z)))
    temp11 = temp11[NN:end] .* exp.(1im .* kdz) .* (-1 .* 1im .* misc.RTC.khi_eff .* misc.RTC.omega .^ 2 / 2 / misc.NC.c0^2 ./ dkdz) .* misc.RTC.domega - 1 .* thzAbsorption / 2 .* ATHz
    temp11[1] = 0
    return temp11
  end

  t2 = @spawn begin
    temp21 = conv(reverse(conj(ATHz) .* exp.(1im .* kdz)), Aop .* exp.(-1im .* misc.RTC.k_omega .* z))
    temp21 = temp21[NN:end] .* exp.(1im .* misc.RTC.k_omega .* z)
    temp22 = conv(Aop .* exp.(-1im .* misc.RTC.k_omega .* z), ATHz .* exp.(-1im .* kdz))
    temp22 = temp22[1:NN] .* exp.(1im .* misc.RTC.k_omega .* z)
    temp20 = -1 * n2pm - 1 * 1im * misc.RTC.khi_eff .* misc.RTC.omega .^ 2 / 2 / misc.NC.c0^2 ./ misc.RTC.k_omega .* (temp21 + temp22) .* misc.RTC.domega
    temp20[1] = 0
    #= temp23 = conv(reverse(conj(Aop) .* exp.(1im .* misc.RTC.k_omega .* z)), ASH .* exp.(-1im .* misc.RTC.k_omegaSH .* z)) .* misc.RTC.domega
    temp23 = -1 ./ cos(misc.RTC.gamma) .* 1im .* misc.RTC.deff .* misc.RTC.omega .^ 2 / misc.NC.c0^2 ./ misc.RTC.k_omega .* temp23[NN:end] .* exp.(1im .* misc.RTC.k_omega .* z)
    temp24 = temp20 + temp23
    temp24[1] = 0=#
    return temp20
  end

  #=t3 = @spawn begin
    temp31 = conv(Aop .* exp.(-1im .* misc.RTC.k_omega .* z), Aop .* exp.(-1im .* misc.RTC.k_omega .* z)) * misc.RTC.domega
    temp31 = -1 ./ cos(misc.RTC.gamma) * 1im * misc.RTC.deff .* misc.RTC.omega .^ 2 / 2 / misc.NC.c0^2 ./ misc.RTC.k_omega .* temp31[1:NN] .* exp.(1im .* misc.RTC.k_omegaSH .* z)
    temp31[1] = 0
    return temp31
  end =#
  wait.([t1, t2])
  return differentialEqInputs(ATHz=t1.result, Aop=t2.result, Nc=Neff, cumulativePhase=dkdz)
end

function RK4_M(f::Function, step::Float64, T::Float64, Y::differentialEqInputs, misc::miscInput)::differentialEqInputs
  k1 = f(T, Y, misc)
  k2 = f(T + step / 2, Y + k1 * step / 2, misc)
  k3 = f(T + step / 2, Y + k2 * step / 2, misc)
  k4 = f(T + step, Y + k3 * step, misc)
  Y += 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4) * step
  return Y
end

function n2value(cry)
  if cry == 4 # GaAs
    n2_ = 5.9e-18 #2023-08-04
  elseif cry == 3
    n2_ = 11.5e-18 #https://doi.org/10.1038/s41566-019-0537-9
  end
  return n2_
end

function deff(cry)
  if cry == 4 #% GaAs
    deff_ = 2 / sqrt(3) * 80e-12 #2023-08-04
  end
  return deff_
end

function aTHzo(omega, T, cry)
  if cry == 4
    alpha = -2 .* omega / 3e8 .* imag(sqrt.(er(omega, T, cry)))
  end
  return alpha
end

function er(omega, T, cry, Nc=0)
  nu = omega / 2 / pi / 3e8 * 0.01
  if cry == 4 #GaAs
    if T == 300 #ord
      e_inf = 11
      nu_L = 292.1
      nu_T = 268.7
      G = 2.4
      nu = omega / 2 / pi / 3e8 * 1e-2

      er_ = e_inf * (1 .+ (nu_L^2 .- nu_T^2) ./ (nu_T^2 .- nu .^ 2 .+ 1im * G * nu))
    end
  end

  if cry == 3 # GaP
    tsc = 180e-15
    meff = 0.25 * 9.109e-31
    q = 1.602e-19
    e0 = 8.8541878e-12

    op2 = q^2 * Nc / e0 / meff

    er_ =@. 9.09 + 2.06 * 363.4^2 ./ (363.4^2 - omega .^ 2 * 0.01^2 / (2 * pi * 3e8)^2 - 2 * 1im * 0.55 * omega * 0.01 / (2 * pi * 3e8))
    -op2 ./ (omega .^ 2 + 1im * omega / tsc)
  end
  return er_
end

function DataBaseWriter(FID, z, Aop, Iop, ATHz, ETHz, ASH, ISH)
  FID[string(Int(round(z * 1e6)))*"/Aop"] = collect(abs.(Aop))
  FID[string(Int(round(z * 1e6)))*"/Eop"] = collect(Iop)
  FID[string(Int(round(z * 1e6)))*"/ATHz"] = collect(abs.(ATHz))
  FID[string(Int(round(z * 1e6)))*"/ETHz"] = collect(ETHz)
  FID[string(Int(round(z * 1e6)))*"/ASH"] = collect(abs.(ASH))
  FID[string(Int(round(z * 1e6)))*"/ESH"] = collect(ISH)
end

function DataBaseEnder(FID, z, t, nu, effic, efficSH)
  FID["z"] = collect(transpose(z))
  FID["effic"] = collect(transpose(effic))
  FID["efficSH"] = collect(transpose(efficSH))
  FID["t"] = collect(transpose(t))
  FID["nu"] = collect(transpose(nu))
end

