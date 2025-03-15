using Comrade
using Pyehtim
using Krang
using Pigeons
import CairoMakie as CM

include(joinpath((@__DIR__),"models","JuKeBOX.jl"))
include(joinpath((@__DIR__),"models","modifiers.jl"))

seed = 1234
phasecal = true
ampcal = true
add_th_noise = true
scan_avg = true
fractional_noise = 0.01
#npix = 180
npix = 180
n_tempering_levels = 20
#fovx = μas2rad(120.0)
#fovy = μas2rad(120.0)
fovx = μas2rad(120.0)
fovy = μas2rad(120.0)

function ModifiedJuKeBOX(θ, meta) 
    RenormalizedFlux(modify(JuKeBOX(θ), Stretch((θ.m_d), (θ.m_d)), Rotate(θ.pa)), 0.6)
end

inobs = ehtim.obsdata.load_uvfits(joinpath(dirname(@__DIR__), "data", "withbhex.uvfits"))
obs = inobs
obs = scan_average(obs).flag_uvdist(uv_min=0.2e9)
obs = obs.add_fractional_noise(fractional_noise)

dvis = extract_table(obs, Visibilities())
dvisamp = extract_table(obs, VisibilityAmplitudes())
dcphase = extract_table(obs, ClosurePhases(;snrcut=3.0))
dlcamp = extract_table(obs, LogClosureAmplitudes(;snrcut=3.0))
#fig = CM.Figure(;size=(800, 600));
#axisfields(fig[1,1], dvis, :uvdist, :measurement)
#axisfields(fig[1,2], dvisamp, :uvdist, :measurement)
#axisfields(fig[2,1], dcphase, :uvdist, :measurement)
#axisfields(fig[2,2], dlcamp, :uvdist, :measurement)
#fig 
using VLBIImagePriors
using Distributions

grid = imagepixels(fovx, fovy, npix, npix; executor=ThreadsEx())
priorJB = (
    m_d = Uniform(μas2rad(1.5), μas2rad(8.0)), 
    spin = Uniform(-1.0, -0.01),
    θo = Uniform(1.0, 40.0),
    θs = Uniform(40.0, 90.0),
    pa = Uniform(-π, 0),
    rpeak = Uniform(1.0, 10.0),
    p1 = Uniform(0.1, 10.0),
    p2 = Uniform(1.0, 10.0),
    χ = Uniform(-π, π),
    ι = Uniform(-π/2, π/2),
    βv = Uniform(0.0, 0.99),
    σ = Uniform(-1.0, 5.0),
    η = Uniform(-π, π),
)

skymJB = SkyModel(ModifiedJuKeBOX, priorJB, grid)
intmodel = Comrade.IdealInstrumentModel()
postJB = VLBIPosterior(skymJB, intmodel, dlcamp, dcphase)#; admode=set_runtime_activity(Enzyme.Reverse))

