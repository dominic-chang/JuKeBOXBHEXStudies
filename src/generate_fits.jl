# This script generates a fits file for a JuKeBOX model
using Pkg; Pkg.activate(dirname(@__DIR__));
using Revise
using Krang
using VIDA
using CairoMakie
include(joinpath((@__DIR__), "models", "JuKeBOX.jl"))

source = "JuKeBOX"
ra = 12.5137287172 # Right ascension 
dec = 12.3911232392 # Declination
mjd = 58228 # Modified Julian Date
rf = 227071000000.0 # frequency

fov = μas2rad(120) # In μas
npix = 400
θ = (
    m_d=3.5e0, # Mass to distance ratio in μas
    spin=-0.95e0, # Black hole spin
    pa=0.0e0, # spin position angle
    θo=17.0, # Observer inclination
    θs=75, # Cone opening angle
    rpeak=3.0, # Characteristic radius of emissivity profile
    p1=4.0, # Outer power-law index of emissivity profile
    p2=4.0, # Inner power-law index of emissivity profile
    χ=-π/2, # Azimuthal angle of fluid flow in ZAMO
    ι=0.6, # Inclination of magnetic field in ZAMO
    βv=0.9, # Speed of fluid flow in ZAMO
    σ=0.7, # spectral index of the electron distribution
    η=-1e0π, # Azimuthal angle of the magnetic field in ZAMO
)

mdl = modify(JuKeBOX(θ), Stretch(μas2rad(θ.m_d), μas2rad(θ.m_d)), Rotate(θ.pa*π/180))

grid = VIDA.imagepixels(
    fov, 
    fov, 
    npix, 
    npix, 
    executor=ThreadsEx();
    header = ComradeBase.MinimalHeader(
        source,
        (ra, dec, mjd, rf)...,
    ),
    )
intmap = VIDA.intensitymap(mdl, grid)

imageviz(intmap)

outpath = joinpath(dirname(@__DIR__), "data", "fits")
mkpath(outpath)
save_fits(joinpath(outpath, "JuKeBOX.fits"), intmap)