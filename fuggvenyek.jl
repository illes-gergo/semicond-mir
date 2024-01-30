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

    n = real(sqrt.(Complex.(a0 .+ 1 .+ a[1] * l .^ 2 ./ (l .^ 2 .- b[1]^2) + a[2] * l .^ 2 ./ (l .^ 2 .- b[2]^2) + a[(3)] * l .^ 2 ./ (l .^ 2 .- b[(3)]^2))))
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

  end
  return ng
end

function nTHzo(omega, T, cry)
  if cry == 4

    nTHz = real.(sqrt.(er(omega, T, cry)))
  end
  return nTHz
end

function diffegy_conv(z, A_kompozit::differentialEqInputs, misc::miscInput)#, omega, T, k_omega, k_OMEGA, k_omegaSH, khi_eff, dnu, domega, k_omega0, omega0, gamma, cry)
  #n2 = _differential_equation == 1 || _differential_equation == 3 ? n2value(misc.cry) : 0
  #deff_ = _differential_equation == 2 || _differential_equation == 3 ? 1 * deff(cry) : 0
  #c = 3e8    #%m/s
  #e0 = 8.854e-12
  #abszorpcio = aTHzo(omega, T, cry)
  #abszorpcio[abszorpcio.>1e5] .= 1e5
  ATHz = A_kompozit.ATHz
  Aop = A_kompozit.Aop
  ASH = A_kompozit.ASH
  NN = length(Aop)
  At = ifft(Aop .* exp.(-0im * (misc.RTC.k_omega) * z) * 2 * pi * misc.RTC.dnu * length(misc.RTC.omega)) # ???? Szerintem -1im kell
  n2pm = fft(1im * misc.NC.e0 * misc.RTC.omega0 * neo(2 * pi * 3e8 / misc.RTC.omega0, misc.IN.T, misc.IN.cry) * misc.RTC.n2 / 2 * abs.(At) .^ 2 .* At) / misc.RTC.dnu / 2 / pi / length(misc.RTC.omega) .* exp.(0im .* misc.RTC.k_omega .* z)
  t1 = @spawn begin
    temp11 = conv(reverse(conj(Aop) .* exp.(1im .* misc.RTC.k_omega .* z)), (Aop .* exp.(-1im * misc.RTC.k_omega .* z)))
    temp11 = temp11[NN:end] .* exp.(1im .* misc.RTC.k_OMEGA .* z) .* (-1 .* 1im .* misc.RTC.khi_eff .* misc.RTC.omega .^ 2 / 2 / misc.NC.c0^2 ./ misc.RTC.k_OMEGA) .* misc.RTC.domega - 1 .* misc.RTC.absorption / 2 .* ATHz
    temp11[1] = 0
    return temp11
  end

  t2 = @spawn begin
    #       println(threadid())
    temp21 = conv(reverse(conj(ATHz) .* exp.(1im .* misc.RTC.k_OMEGA .* z)), Aop .* exp.(-1im .* misc.RTC.k_omega .* z))
    temp21 = temp21[NN:end] .* exp.(1im .* misc.RTC.k_omega .* z)
    temp22 = conv(Aop .* exp.(-1im .* misc.RTC.k_omega .* z), ATHz .* exp.(-1im .* misc.RTC.k_OMEGA .* z))
    temp22 = temp22[1:NN] .* exp.(1im .* misc.RTC.k_omega .* z)
    temp20 = -1 * n2pm - 1 * 1im * misc.RTC.khi_eff .* misc.RTC.omega .^ 2 / 2 / misc.NC.c0^2 ./ misc.RTC.k_omega .* (temp21 + temp22) .* misc.RTC.domega
    temp20[1] = 0
    temp23 = conv(reverse(conj(Aop) .* exp.(1im .* misc.RTC.k_omega .* z)), ASH .* exp.(-1im .* misc.RTC.k_omegaSH .* z)) .* misc.RTC.domega
    temp23 = -1 ./ cos(misc.RTC.gamma) .* 1im .* misc.RTC.deff .* misc.RTC.omega .^ 2 / misc.NC.c0^2 ./ misc.RTC.k_omega .* temp23[NN:end] .* exp.(1im .* misc.RTC.k_omega .* z)
    temp24 = temp20 + temp23
    temp24[1] = 0
    return temp24
  end
  t3 = @spawn begin
    temp31 = conv(Aop .* exp.(-1im .* misc.RTC.k_omega .* z), Aop .* exp.(-1im .* misc.RTC.k_omega .* z)) * misc.RTC.domega
    temp31 = -1 ./ cos(misc.RTC.gamma) * 1im * misc.RTC.deff .* misc.RTC.omega .^ 2 / 2 / misc.NC.c0^2 ./ misc.RTC.k_omega .* temp31[1:NN] .* exp.(1im .* misc.RTC.k_omegaSH .* z)
    temp31[1] = 0
    return temp31
  end
  wait.([t1, t2, t3])
  return differentialEqInputs(ATHz=t1.result, Aop=t2.result, ASH=t3.result)
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
function er(omega, T, cry)
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

