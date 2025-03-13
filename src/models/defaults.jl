abstract type AbstractModelInfo end
struct ModelInfo{M,T,F} <: AbstractModelInfo
    model::M
    name::String
    lower::NamedTuple{T}
    upper::NamedTuple{T}
    file::Union{String, Vector{String}}
    blur::F
    function ModelInfo(model::M, name::String, lower::NamedTuple{T}, upper::NamedTuple{T}, file::Union{String,Vector{String}}, blur::F=zero(F)) where {M,T,F}
        new{M,T,F}(model, name, lower, upper, file, blur)
    end
end

