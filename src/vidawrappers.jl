using OptimizationOptimisers

function blur_fwhm(img::ComradeBase.IntensityMap, fwhm)
    dx, dy = Float64.(rad2μas.([pixelsizes(inimg)...])) 
    σ_px = fwhm./(2*sqrt(2*log(2e0)))./abs(dx)
    σ_py = fwhm./(2*sqrt(2*log(2e0)))./abs(dy)

    # Now I need to pick my kernel size. I am going out to 5σ for the
    # gaussian kernel. I have to add one for the convolution to play nice
    nkern = Int(floor(σ_px)*10 + 1)
    # Note I tried to use IIRGaussian but it wasn't accurate enough for us.
    return imfilter(img,
                    Kernel.gaussian((σ_py, σ_px),(nkern,nkern)),
                    Fill(0.0, img),
                    Algorithm.FFT()Algorithm.FFT()
                )
end

function imgfit64(info::ModelInfo, maxevals::Int; will_overwrite::Bool=false)
    (; model, name, lower, upper, file, blur) = info
    iobase = joinpath(pwd(), "runs", "image_domain", split(file, "/")[end], name *"_$blur")
    if (!will_overwrite) && isdir(iobase)
        return
    end

    mkpath(iobase)
    println("Saving output at :" * iobase)
    println("lower :$lower")
    println("upper :$upper")

    #T = Float64
    #inimg = VIDA.load_image(file)
    inimg = regrid(VIDA.load_image(file), imagepixels((μas2rad(80.0)), (μas2rad(80.0)), 240, 240))
    dx, dy = ([pixelsizes(inimg)...])
    #g = axisdims(inimg)
    vals = zeros(Float64, size(inimg))

    if iszero(blur)
        vals = Float64.(parent(inimg))
    else
        σ_px = μas2rad(blur) / (2 * sqrt(2 * log(2))) / abs(dx)
        σ_py = μas2rad(blur) / (2 * sqrt(2 * log(2))) / abs(dy)

        nkern = Int(floor(σ_px) * 10 + 1)
        vals = imfilter(parent(inimg),
            Kernel.gaussian((σ_py, σ_px), (nkern, nkern)),
            Fill(0e0, inimg),
            Algorithm.FFT()
        )
    end
    #g64 = RectiGrid((X=map(Float64, (g.X)), Y=map(Float64, (g.Y))), executor=ThreadsEx())
    g64 = rebuild(typeof(axisdims(inimg)), dims(inimg), ThreadsEx())
    img = max.(IntensityMap(vals, g64), 0.0) #|>MtlArray

    nx = VIDA.NxCorr(img)
    lower = map(Float64, lower)
    upper = map(Float64, upper)

    println(typeof(model))  
    println(typeof(lower))
    println(typeof(upper))
    prob = VIDAProblem(nx, model, lower, upper)
    xopt, opt_temp, divmin = vida(prob, BBO_adaptive_de_rand_1_bin(); maxiters=maxevals, verbose=true)


    fig = triptic(img, opt_temp)

    newkeys = fieldnames(typeof(xopt))
    opt = getfield.(Ref(xopt), newkeys)
    intmap = intensitymap(model(xopt), g64)
    intmap |> Plots.plot
    VIDA.save_fits(iobase * "/best.fits", intmap)

    fileout = open(iobase * "/best_nxcorr.txt", "w")

    divmin = divergence(nx, opt_temp)
    write(fileout, "nxcorr = " * string(exp(-divmin)) * "\n")
    write(fileout, "upper = " * string(upper) * "\n")
    write(fileout, "lower = " * string(lower) * "\n")
    write(fileout, "best_fit = " * string(NamedTuple{newkeys}(opt)))
    close(fileout)

    save(joinpath(iobase, "tripic.png"), fig)

end