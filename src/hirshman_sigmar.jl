using SpecialFunctions

function build_parameter_matrices(dd::IMAS.dd)
	eqt = dd.equilibrium.time_slice[]
	cp1d = dd.core_profiles.profiles_1d[1]

    m_to_cm = IMAS.gacode_units.m_to_cm

	rmin = IMAS.r_min_core_profiles(cp1d, eqt)
	a = rmin[end]

	num_ions = length(cp1d.ion)

    e = 4.8e-10
    k = 1.6e-12
   
	loglam = 24.0 .- log.(sqrt.((cp1d.electrons.density ./ 1e6) ./ (cp1d.electrons.temperature ./ 1e3)))

    n_norm = cp1d.electrons.density ./ 1e6 
    t_norm = cp1d.electrons.temperature ./ 1e3

    m_norm = 2.014102 # mass of deuterium in atomic mass units
    nu_norm = sqrt.(t_norm ./ m_norm) ./ a

	Z = Vector{Float64}(undef, num_ions+1)
	mass = Vector{Float64}(undef, num_ions+1)


    dens = zeros(Float64, length(cp1d.ion[1].density), num_ions+1)
	temp = zeros(Float64, length(cp1d.ion[1].temperature), num_ions+1)
    vth = zeros(Float64, length(cp1d.ion[1].temperature), num_ions+1)
    nu = zeros(Float64, length(cp1d.ion[1].temperature), num_ions+1)


    dlnndr = zeros(Float64, length(cp1d.ion[1].density), num_ions+1)
    dlntdr = zeros(Float64, length(cp1d.ion[1].temperature), num_ions+1)

	for i in 1:num_ions
		Z[i] = cp1d.ion[i].element[1].z_n
		mass[i] = cp1d.ion[i].element[1].a ./ m_norm

		dens[:,i] = cp1d.ion[i].density ./ 1e6 ./ n_norm
		temp[:,i] = cp1d.ion[i].temperature ./ 1e3 ./ t_norm

		nu[:,i] = (@. sqrt(2) * pi * dens[:,i] * Z[i]^4.0 * e^4.0 * loglam / sqrt(1.6726e-24 .* mass[i]) / (k * temp[:,i])^1.5) ./ nu_norm

		# gradients - these are finally actually right so don't change anything
		dlnndr[:,i] = -IMAS.calc_z(rmin / a , cp1d.ion[i].density)
		dlntdr[:,i] = -IMAS.calc_z(rmin / a, cp1d.ion[i].temperature)



        vth[:,i] = sqrt.(cp1d.ion[i].temperature ./ 1e3 ./ cp1d.ion[i].element[1].a)
	end

    # tacking on electron parameters at the end 
    Z[end] = -1.0
    mass[end] = 0.00054858 / m_norm # 0.00054858 is the mass of an electron in AMU
    dens[:,end] = cp1d.electrons.density ./ 1e6 ./ n_norm 
    temp[:,end] = cp1d.electrons.temperature ./ 1e3 ./ t_norm

    nu[:,end] = (@. sqrt(2) * pi * dens[:,end] * Z[end]^4.0 * e^4.0 * loglam / sqrt(1.6726e-24 .* mass[end]) / (k .* temp[:,end])^1.5) ./ nu_norm

    dlnndr[:,end] = -IMAS.calc_z(rmin / a, cp1d.electrons.density_thermal)
    dlntdr[:,end] = -IMAS.calc_z(rmin / a, cp1d.electrons.temperature) # this might still be wrong

    vth[:,end] = sqrt.(cp1d.electrons.temperature ./ 1e3 ./ 0.00054858)

    return Z, mass, dens, temp, nu, dlnndr, dlntdr, vth
end


function gauss_legendre(x1::Int, x2::Int, n::Int)
	eps = 1e-15
	pi = 3.141592653589793

	xm = 0.5 * (x2 + x1)
	xl = 0.5 * (x2 - x1)

	x = zeros(Float64, n)
	w = zeros(Float64, n)

	# Exception for n=1 is required:
	if n == 1
		x[1] = xm
		w[1] = 2.0 * xl
		return x, w
	end

	# Roots are symmetric. We only need to find half of them.
	m = (n + 1) / 2

	# Initialize to fail first do test
	z1 = -1.0
	pp = 0.0

	# Loop over roots.
	for i in 1:m
		i = Int(i)
		z = cos(pi * (i - 0.25) / (n + 0.5))

		while abs(z - z1) > eps
			p1 = 1.0
			p2 = 0.0

			for j in 1:n
				p3 = p2
				p2 = p1
				p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
                # println("The value of p1 is:", p1)
            end

			# p1 is the Legendre polynomial. Now compute its derivative, pp.
			pp = n * (z * p1 - p2) / (z * z - 1.0)
			z1 = z
			z = z1 - p1 / pp
		end

		x[i]     = xm - xl * z
		x[n+1-i] = xm + xl * z
		w[i]     = 2.0 * xl / ((1.0 - z * z) * pp * pp)
		w[n+1-i] = w[i]
	end

	return x, w
end

function gauss_integ(xmin, xmax, func, order, n_subdiv)

	x0, w0 = gauss_legendre(0, 1, order)

	dx = (xmax - xmin) / n_subdiv

	n_node = n_subdiv * order

	x = zeros(Float64, n_subdiv * order)
	w = zeros(Float64, n_subdiv * order)

	for i in 1:n_subdiv
		for j in 1:order
			p = (i - 1) * order + j
			x[p] = xmin + ((i - 1) + x0[j]) * dx
			w[p] = w0[j] * dx
		end
	end

	answer = 0.0
	for p in 1:n_node
		answer = answer + w[p] * func(x[p])
	end

	return answer
end

function get_coll_freqs(ir_loc, is_loc, js_loc, ene, dd::IMAS.dd)
    Z, mass, dens, temp, nu, dlnndr, dlntdr, vth = NEO.build_parameter_matrices(dd)
    
    fac = (1.0 * Z[js_loc]^2 / (1.0 * Z[is_loc])^2 * dens[:,js_loc][ir_loc]) / dens[:,is_loc][ir_loc]

	xa = sqrt(ene)
    xb = xa * (vth[:,is_loc][ir_loc]/vth[:,js_loc][ir_loc])

    if xb < 1e-4 
        nu_d = fac * (1.0/sqrt(pi)) * (4.0/3.0 * (vth[:,is_loc][ir_loc]/vth[:,js_loc][ir_loc]) - 4.0/15.0  * (vth[:,is_loc][ir_loc]/vth[:,js_loc][ir_loc])^3 * ene + 2.0/35.0  * (vth[:,is_loc][ir_loc]/vth[:,js_loc][ir_loc])^5 * ene^2 - 2.0/189.0 * (vth[:,is_loc][ir_loc]/vth[:,js_loc][ir_loc])^7 * ene^3)

    else
        Hd_coll = exp(-xb * xb) / (xb * sqrt(pi)) + (1 - (1 / (2 * xb * xb))) * erf(xb)
        Xd_coll = 1 / xa
        nu_d = fac * Hd_coll * Xd_coll
    end

	return nu_d
end

function myHSenefunc(x)
    cp1d = ddFunc.core_profiles.profiles_1d[1]
    eq = ddFunc.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    Z, mass, dens, temp, nu, dlnndr, dlntdr, vth = NEO.build_parameter_matrices(ddFunc)

    m_to_cm = IMAS.gacode_units.m_to_cm


	n_species = length(cp1d.ion)+1

	emin = 0.0
	emax = 16.0

	xa = 2.0 / (1.0 - sqrt(emin / emax))
	xb = -(1.0 + sqrt(emin / emax)) / (1.0 - sqrt(emin / emax))
	ene = emax * ((x - xb) / xa)^2
	de = 2.0 * sqrt(emax) / xa
	val = de * exp(-ene)

    rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    rmin = IMAS.r_min_core_profiles(cp1d, eqt) # the function r_min_core_profiles already does the conversion from m_to_cm
    q = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)

    eps = rmin[ir_global] / rmaj[ir_global]

	nu_d_tot = 0.0
	for js in 1:n_species
		nu_d = get_coll_freqs(ir_global,is_globalFunc,js,ene,ddFunc)
		nu_d_tot += nu_d * nu[:,is_globalFunc][ir_global]  
        # println("nu_d_tot = ", nu_d_tot)
	end

	ft_star = (3.0 * pi / 16.0) * eps^2 * vth[:,is_globalFunc][ir_global] * sqrt(2.0) * ene^1.5 / (rmaj[ir_global] * abs(q[ir_global]) * nu_d_tot)
    ftrap = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.trapped_fraction).(cp1d.grid.rho_tor_norm)
    ft_fac = 1.0 / (1.0 + ftrap[ir_global] / ft_star)

    if ietypeFunc == 1
	    myHSenefunc = val * nu_d_tot * ene * ft_fac
    elseif ietypeFunc == 2
        myHSenefunc = val * nu_d_tot * ene * ene * ft_fac
    elseif ietypeFunc == 3
        myHSenefunc = val * nu_d_tot * ene * ene * ene * ft_fac
    end

	return myHSenefunc

end

Bmag2_avg = 0.68

function compute_HS(ir, dd::IMAS.dd)
    global ddFunc = dd

    Z, mass, dens, temp, nu, dlnndr, dlntdr, vth = NEO.build_parameter_matrices(ddFunc)

	#uses mass, density, dlndndr, dlntdr, n_species, q, r 
	eqt = ddFunc.equilibrium.time_slice[]
    eq1d = eqt.profiles_1d
	cp1d = ddFunc.core_profiles.profiles_1d[1]

	n_species = length(cp1d.ion) + 1
    m_to_cm = IMAS.gacode_units.m_to_cm

    ftrap = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.trapped_fraction).(cp1d.grid.rho_tor_norm)

    rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
    rmin = IMAS.r_min_core_profiles(cp1d, eqt)


    rho = cp1d.grid.rho_tor_norm
    q = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)

    global ir_global = ir 

	Nx = 100
	integ_order = 8
	omega_fac = 1.0 / Bmag2_avg
	HS_I_div_psip = rmaj[ir] * q[ir] / rmin[ir] 

	nux0 = zeros(Float64, n_species)
	nux2 = zeros(Float64, n_species)
	nux4 = zeros(Float64, n_species)


	for is_global in 1:n_species
		for ietype in 1:3
            global ietypeFunc = ietype
            global is_globalFunc = is_global

			eii_val = gauss_integ(-1.0, 1.0, NEO.myHSenefunc, integ_order, Nx)
            # println("eii_val = ", eii_val)

			if ietype == 1
				nux0[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
			elseif ietype == 2
				nux2[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
			elseif ietype == 3
				nux4[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
			end
		end
	end

	sum_nm = 0.0
	for is_global in 1:n_species
		sum_nm += mass[is_global] * dens[:,is_global][ir] * nux0[is_global]
	end

	pflux_multi = zeros(Float64, n_species)
	eflux_multi = zeros(Float64, n_species)
	for is_global in 1:n_species
		A1 = -dlnndr[:,is_global][ir] + (1.5 * dlntdr[:,is_global][ir])
        # println("A1 = ", A1)
		A2 = -dlntdr[:,is_global][ir]
        # println("A2 = ", A2)

		pflux_multi[is_global] = 0.0
		eflux_multi[is_global] = 0.0

		L_a = nux0[is_global] * omega_fac * HS_I_div_psip^2 * rho[ir]^2 * ftrap[ir] * dens[:,is_global][ir] * temp[:,is_global][ir] * mass[is_global] / (Z[is_global] * Z[is_global] * 1.0)
        
		for js in 1:n_species

			L_b = nux0[js] * omega_fac * HS_I_div_psip^2 * rho[ir]^2 * ftrap[ir] * dens[:,js][ir] * temp[:,js][ir] * mass[js] / (Z[js] * Z[js] * 1.0)
            # println("L_b = ", L_b)

			if is_global == js
				L11 = -L_a * (sum_nm - mass[is_global] * dens[:,is_global][ir] * nux0[is_global]) / sum_nm
				L12 = L11 * (nux2[is_global] / nux0[is_global])
				L22 = (nux2[is_global] / nux0[is_global]) * (L12 + L_a * (nux2[is_global] / nux0[is_global] - nux4[is_global] / nux2[is_global]))
				L21 = L12
			else
				L11 = L_a * (Z[is_global] * temp[:,js][ir]) / (Z[js] * temp[:,is_global][ir]) * mass[js] * dens[:,js][ir] * nux0[js] / sum_nm
				L12 = (nux2[js] / nux0[js]) * L11
				L21 = (nux2[is_global] / nux0[is_global]) * Z[js] / (1.0 * Z[is_global]) * (Z[is_global] * dens[:,is_global][ir] * nux0[is_global] / sum_nm) * L_b
				L22 = (nux2[is_global] / nux0[is_global]) * L12

			end

			pflux_multi[is_global] += L11 * A1 + L12 * A2
			eflux_multi[is_global] += L21 * A1 + L22 * A2

		end


	end

	return  pflux_multi, eflux_multi

end