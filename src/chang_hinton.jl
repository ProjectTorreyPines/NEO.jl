"""
    neoclassical_changhinton(
        eq1d::IMAS.equilibrium__time_slice___profiles_1d, 
        cp1d::IMAS.core_profiles__profiles_1d,
        rho_fluxmatch::Real,
        iion::Integer)

Calculates the neoclassical flux using Chang-Hinton model which has has been modified assuming Zi = 1, and ni=ne
"""
function changhinton(
    eqt::IMAS.equilibrium__time_slice,
    cp1d::IMAS.core_profiles__profiles_1d,
    rho_fluxmatch::Real,
    iion::Integer)

    eq1d = eqt.profiles_1d

    rmin = 0.5 * (eq1d.r_outboard - eq1d.r_inboard)

    m_to_cm = 1e2
    Rmaj0 = eq1d.geometric_axis.r[1] * m_to_cm

    a = eqt.boundary.minor_radius * m_to_cm
    rho_cp = cp1d.grid.rho_tor_norm
    rho_eq = eq1d.rho_tor_norm
    gridpoint_eq = argmin(abs.(rho_eq .- rho_fluxmatch))
    gridpoint_cp = argmin(abs.(rho_cp .- rho_fluxmatch))
    eps = rmin[gridpoint_eq] * m_to_cm / Rmaj0

    q = eq1d.q[gridpoint_eq]

    ne = cp1d.electrons.density_thermal * 1e-6
    Te = cp1d.electrons.temperature[gridpoint_cp]
    Ti = cp1d.ion[iion].temperature

    rmin = IMAS.interp1d(rho_eq, rmin).(rho_cp)
    Rmaj = IMAS.interp1d(eq1d.rho_tor_norm, 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    drmaj = IMAS.gradient(rmin, Rmaj)
    shift = -drmaj[gridpoint_cp]

    dlntdr = -IMAS.gradient(rmin, Ti)[gridpoint_cp] / Ti[gridpoint_cp]
    dlntdr = dlntdr * a / m_to_cm
    dlnndr = -IMAS.gradient(rmin, ne)[gridpoint_cp] / ne[gridpoint_cp]
    dlnndr = dlnndr * a / m_to_cm
    Ti = Ti[gridpoint_cp]
    ne = ne[gridpoint_cp]
    k = 1.6e-12
    e = 4.8e-10

    mi = 1.6726e-24 * cp1d.ion[iion].element[1].a

    c_s = sqrt(k * Te / mi)
    loglam = 24.0 - log(sqrt(ne / Te))

    k0 = 0.66
    a0 = 1.03
    b0 = 0.31
    c0 = 0.74

    Zi = IMAS.avgZ(cp1d.ion[iion].element[1].z_n, Ti)
    alpha = Zi - 1.0
    nui = sqrt(2.0) * pi * ne * (Zi * e)^4 * loglam / sqrt(mi) / (k * Ti)^1.5
    nu = nui * a / c_s / sqrt(Ti / Ti) * (Ti / Te)^1.5
    nui_HH = nu * (4.0 / 3.0) / sqrt(2.0 * pi)
    nui_star_HH = nui_HH * (Rmaj0 / a) * abs(q / (sqrt(eps) * sqrt(eps) * sqrt(eps)))
    mu_star = (1.0 + 1.54 * alpha) * nui_star_HH

    CH_Bmag2inv_avg = ((1.0 + 1.5 * (eps * eps + eps * shift)
                        + 0.375 * eps * eps * eps * shift)
                       /
                       (1.0 + 0.5 * eps * shift))
    CH_Bmag2avg_inv = ((sqrt(1.0 - eps * eps) * (1.0 + 0.5 * eps * shift))
                       /
                       (1.0 + (shift / eps) * (sqrt(1.0 - eps * eps) - 1.0)))
    CH_I_div_psip = q / eps

    F2 = (0.5 / sqrt(eps)) * (CH_Bmag2inv_avg - CH_Bmag2avg_inv)

    K1 = (-alpha * (0.83 + 0.42 * alpha) / (0.58 + alpha)
          *
          (c0 * mu_star * sqrt(eps * eps * eps) * F2)
          /
          (1.0 + c0 * mu_star * sqrt(eps * eps * eps)))

    K2 = (
        (k0 * (1.0 + 1.54 * alpha)
         +
         (1.88 * sqrt(eps) - 1.54 * eps) * (1.0 + 3.75 * alpha)) * CH_Bmag2inv_avg
        /
        (1 + a0 * sqrt(mu_star) + b0 * mu_star)
        +
        k0 * sqrt(eps * eps * eps) * (c0 * c0 / b0) * mu_star * F2
        * (1.0 + 1.33 * alpha * (1.0 + 0.6 * alpha) / (1.0 + 1.79 * alpha))
        /
        (1.0 + c0 * sqrt(eps * eps * eps) * mu_star)
    )

    neo_rho_star_in = 0.001

    efluxi = (CH_I_div_psip^2
              *
              Ti / Te
              * (cp1d.ion[iion].element[1].a / Zi^2 * neo_rho_star_in^2 * sqrt(eps) * nui_HH)
              * ((K2 + K1) * dlntdr + K1 * dlnndr))

    qneo_gb = neo_rho_star_in^2

    sol = IMAS.flux_solution(0.0, 0.0, 0.0, efluxi / qneo_gb)
    return sol
end