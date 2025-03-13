struct JuKeBOX{T,F} <: ComradeBase.AbstractModel
    spin::T
    θo::T
    scene::F
end

function JuKeBOX(θ::NamedTuple)
    (; spin, θo, θs, rpeak, p1, p2, χ, ι, βv, σ, η) = θ
    T = typeof(spin)


    magfield1 = Krang.SVector(sin(ι) * cos(η), sin(ι) * sin(η), cos(ι))
    magfield2 = Krang.SVector(-sin(ι) * cos(η), -sin(ι) * sin(η), cos(ι))
    vel = Krang.SVector(βv, T(π / 2), χ)
    nlist = (0, 1)
    material1 = Krang.ElectronSynchrotronPowerLawIntensity(
        magfield1...,
        vel...,
        σ,
        rpeak,
        p1,
        p2,
        nlist,
    )
    material2 = Krang.ElectronSynchrotronPowerLawIntensity(
        magfield2...,
        vel...,
        σ,
        rpeak,
        p1,
        p2,
        nlist,
    )

    geometry1 = Krang.ConeGeometry(θs * π / 180)
    geometry2 = Krang.ConeGeometry((1.0f0π - θs * π / 180))

    scene =
        Krang.Scene((Krang.Mesh(geometry1, material1), Krang.Mesh(geometry2, material2)))

    return JuKeBOX(spin, θo, scene)
end

@inline function ComradeBase.intensity_point(m::JuKeBOX{T}, p) where {T}
    (; X, Y) = p
    (; scene, θo) = m

    pix = Krang.IntensityPixel(Krang.Kerr(m.spin), -X, Y, θo * T(π / 180))
    ans = Krang.render(pix, scene)
    return ans + eps(T)
end
