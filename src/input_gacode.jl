Base.@kwdef mutable struct InputGACODE
    NION::Union{Float64,Missing} = missing
    NEXP::Union{Int,Missing} = missing 
    RHO::Union{Vector{Float64},Missing} = missing 
    RMIN::Union{Vector{Float64},Missing} = missing 
    POLFLUX::Union{Vector{Float64},Missing} = missing 
    Q::Union{Vector{Float64},Missing} = missing 
    OMEGA0::Union{Vector{Float64},Missing} = missing 

    RMAJ::Union{Vector{Float64},Missing} = missing 
    ZMAG::Union{Vector{Float64},Missing} = missing
    KAPPA::Union{Vector{Float64},Missing} = missing
    DELTA::Union{Vector{Float64},Missing} = missing
    ZETA::Union{Vector{Float64},Missing} = missing
    SHAPE_COS0::Union{Vector{Float64},Missing} = missing
    SHAPE_COS1::Union{Vector{Float64},Missing} = missing
    SHAPE_COS2::Union{Vector{Float64},Missing} = missing
    SHAPE_COS3::Union{Vector{Float64},Missing} = missing
    SHAPE_SIN3::Union{Vector{Float64},Missing} = missing

    NE::Union{Vector{Float64},Missing} = missing
    TE::Union{Vector{Float64},Missing} = missing
    PTOT::Union{Vector{Float64},Missing} = missing
    ZEFF::Union{Vector{Float64},Missing} = missing
    NI::Union{Vector{Vector{Float64}},Missing} = missing
    TI::Union{Vector{Vector{Float64}},Missing} = missing

    JBS::Union{Vector{Float64},Missing} = missing
    JRF::Union{Vector{Float64},Missing} = missing
    JNB::Union{Vector{Float64},Missing} = missing
    JBSTOR::Union{Vector{Float64},Missing} = missing
    VTOR::Union{Vector{Vector{Float64}},Missing} = missing
    VPOL::Union{Vector{Vector{Float64}},Missing} = missing

    QOHM::Union{Vector{Float64},Missing} = missing
    QBEAME::Union{Vector{Float64},Missing} = missing
    QBEAMI::Union{Vector{Float64},Missing} = missing
    QRFE::Union{Vector{Float64},Missing} = missing
    QRFI::Union{Vector{Float64},Missing} = missing
    QFUSE::Union{Vector{Float64},Missing} = missing
    QFUSI::Union{Vector{Float64},Missing} = missing
    QSYNC::Union{Vector{Float64},Missing} = missing
    QBREM::Union{Vector{Float64},Missing} = missing
    QLINE::Union{Vector{Float64},Missing} = missing
    QEI::Union{Vector{Float64},Missing} = missing
    QIONE::Union{Vector{Float64},Missing} = missing
    QIONI::Union{Vector{Float64},Missing} = missing
    QCXI::Union{Vector{Float64},Missing} = missing

    QPAR_BEAM::Union{Vector{Float64},Missing} = missing
    QPAR_WALL::Union{Vector{Float64},Missing} = missing
    QMOM::Union{Vector{Float64},Missing} = missing

end

function InputGACODE(dd::IMAS.dd)
    input_gacode = InputGACODE()

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    eq1d = eqt.profiles_1d
    cp1d = dd.core_profiles.profiles_1d[]
    ions = cp1d.ion
    cs = dd.core_sources

    input_gacode.NION = length(ions)
    input_gacode.NEXP = length(cp1d.grid.rho_tor_norm)
    input_gacode.RHO = cp1d.grid.rho_tor_norm
    input_gacode.RMIN = IMAS.r_min_core_profiles(cp1d, eqt)
    input_gacode.POLFLUX = eq1d.psi .- eq1d.psi[1]
    input_gacode.Q = eq1d.q
    input_gacode.OMEGA0 = cp1d.rotation_frequency_tor_sonic

    input_gacode.RMAJ = 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)
    # input_gacode.ZMAG = 0.25 * (eq1d.squareness_lower_inner .+ eq1d.squareness_lower_outer .+ eq1d.squareness_upper_inner .+ eq1d.squareness_upper_outer) 
    input_gacode.KAPPA = eq1d.elongation
    input_gacode.DELTA = 0.5 * (eq1d.triangularity_lower .+ eq1d.triangularity_upper)
    input_gacode.ZETA = 0.25 * (eq1d.squareness_lower_inner .+ eq1d.squareness_lower_outer .+ eq1d.squareness_upper_inner .+ eq1d.squareness_upper_outer)


    input_gacode.NE = cp1d.electrons.density_thermal ./ 1e19 
    input_gacode.TE = cp1d.electrons.temperature ./ 1e3 
    input_gacode.PTOT = cp1d.electrons.density_thermal .* cp1d.electrons.temperature .* constants.e 


    input_gacode.JBS = cp1d.j_bootstrap ./ 1e6


    ohmic = findfirst(:ohmic, cs.source)
    input_gacode.QOHM = ohmic.profiles_1d[].electrons.energy ./ 1e6 # change this to qohmE
    
    nbi = findfirst(:nbi, cs.source)
    input_gacode.QBEAME = nbi.profiles_1d[].electrons.energy ./ 1e6
    input_gacode.QBEAMI = nbi.profiles_1d[].total_ion_energy ./ 1e6

    fusion = findfirst(:fusion, cs.source)
    input_gacode.QFUSE = fusion.profiles_1d[].electrons.energy ./ 1e6
    input_gacode.QFUSI = fusion.profiles_1d[].total_ion_energy ./ 1e6

    sync = findfirst(:synchrotron_radiation, cs.source)
    input_gacode.QSYNC = sync.profiles_1d[].electrons.energy ./ 1e6

    brem = findfirst(:bremsstrahlung, cs.source)
    input_gacode.QBREM = brem.profiles_1d[].electrons.energy ./ 1e6

    line = findfirst(:line_radiation, cs.source)
    input_gacode.QLINE = line.profiles_1d[].electrons.energy ./ 1e6 

    return input_gacode

end