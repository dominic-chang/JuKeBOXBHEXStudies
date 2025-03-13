#This script fits images of GRMHD simulations with a blured JuKeBOX model along with the blurring Kernel.
#It uses a NxCorr as a discriminant between the image and the model, which is optimized over globally using a differental evolution algorithm provided by `OptimizationBBO`. 
#The script is parallelized over the blur parameter which controls the amount of blurring applied to the GRMHD image.

using Pkg; Pkg.activate(dirname(@__DIR__));
Pkg.instantiate();
using Revise 
using VIDA
using Plots
using StatsBase
using Comrade
using OptimizationBBO
using ImageFiltering
using Krang
using ComradeBase
using CairoMakie
include(joinpath((@__DIR__), "models", "JuKeBOX.jl"))
include(joinpath((@__DIR__), "models", "defaults.jl"))
include(joinpath((@__DIR__), "vidawrappers.jl"))

inbase = abspath((@__DIR__), "..", "data", "GRMHD")
#List of GRMHD images to fit
files = joinpath.(Ref(inbase),filter(x-> match(r".*fits",x) != nothing, readdir(inbase)))

maxevals = 40_000 # Number of evaluations for the optimizer

function dual_cone_blur(θ) 
    return smoothed(modify(JuKeBOX(θ), Stretch(μas2rad(θ.m_d), μas2rad(θ.m_d)), Rotate(θ.pa*π/180), Shift(μas2rad(θ.x0),μas2rad(θ.y0))), μas2rad(θ.fwhm)/(2√(2*log(2))))
end

lower_dc = (
    m_d=1.5e0,
    spin=-1.00e0,
    pa=0.0e0,
    θo=0.0e0,
    θs=20e0,
    rpeak=1e0,
    p1=0.01e0,
    p2=0.1e0,
    χ=-1.0e0π,
    ι=0e0,
    βv=0.01e0,
    σ=-1.0e0,
    η=-1e0π,
    fwhm=0.0e0,
    x0=-10.0e0,
    y0=-10.0e0,
)

upper_dc = (
    m_d=(8.0e0),
    spin=-.01e0,
    pa=180.0e0,
    θo=40.0e0,
    θs=90e0,
    rpeak=8e0,
    p1=20.0e0,
    p2=10.0e0,
    χ=1.0e0π,
    ι=π / 2.0e0,
    βv=0.9e0,
    σ=3.0e0,
    η=1e0π,
    fwhm=10.0e0,
    x0=10.0e0,
    y0=10.0e0,
)

dists = VIDA._distize(lower_dc, upper_dc)
p_sample = begin 
    pkeys = Tuple(keys(dists))
    vals = [rand(dists[k]) for k in pkeys]
    NamedTuple{pkeys}(vals)
end
p_sample
npix = 240
fov = 120.0
grid = imagepixels(μas2rad(fov), μas2rad(fov), npix, npix)
intmap = intensitymap(dual_cone_blur(p_sample), grid) 
intmap |> imageviz

# Fit GRMHD images blurred with a FWHM of `blur`
Threads.@threads for file in files
    model_info_list = ModelInfo(dual_cone_blur, "JBOX", lower_dc, upper_dc, file, 0e0)
    imgfit64(model_info_list, maxevals;will_overwrite=true)
end
