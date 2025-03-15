preamble_path = joinpath((@__DIR__),"self_fit_preamble.jl")

include(preamble_path)

using Pigeons 

params= (
    exec = "srun",
    submit = `sbatch`,
    del = `scancel`,
    directive = "#SBATCH",
    job_name = "--job-name=",
    output_file = "-o ",
    error_file = "-e ",
    submit_dir = "\$SLURM_SUBMIT_DIR",
    job_status = `squeue --job`,
    job_status_all = `squeue -u`,
    ncpu_info = `sinfo`
)

Pigeons.add_custom_submission_system(params)

function Pigeons.resource_string(m::Pigeons.MPIProcesses, ::Val{:custom})
    return """
    #SBATCH -t $(m.walltime)
    #SBATCH --ntasks=$(m.n_mpi_processes)
    #SBATCH --cpus-per-task=$(m.n_threads)
    #SBATCH --mem-per-cpu=2gb
    """
end

settings = Pigeons.MPISettings(;
submission_system=:custom, 
add_to_submission = [
    "#SBATCH -p blackhole",
    ], 
    environment_modules=["intel","intelmpi"]
)
Pigeons.setup_mpi(settings)

#pt = PT(abspath(dirname(@__DIR__), ".." ,"results", "latest"))
#pt = increment_n_rounds!(pt, 4)
#
#result = pigeons(pt.exec_folder, 
#    on = Pigeons.MPIProcesses(
#        n_mpi_processes = n_tempering_levels,
#    	walltime="10-00:00:00",
#        n_threads = 48,
#        dependencies = [
#            Pigeons, # <- Pigeons itself can be skipped, added automatically
#	    preamble_path
#        ],
#        mpiexec_args=`--mpi=pmi2`
#    )
#)


pt = Pigeons.pigeons(
    target=ascube(postJB), 
    #reference=log_prior; 
    record = [traces, round_trip, Pigeons.timing_extrema], 
    checkpoint=true, 
    n_chains=n_tempering_levels, 
    on = Pigeons.MPIProcesses(
        n_mpi_processes = n_tempering_levels,
    	walltime="14-00:00:00",
        n_threads = 24,
        dependencies = [
            Pigeons, # <- Pigeons itself can be skipped, added automatically
	    preamble_path
        ],
        mpiexec_args=`--mpi=pmi2`
    ),
    n_rounds=20
)
#pt = Pigeons.PT("/n/home06/dochang/2018JuKeBOXFit/results/all/2025-01-23-19-43-36-EBUhaFJN")
#result = pigeons(pt.exec_folder,
#    on = Pigeons.MPIProcesses(
#        n_mpi_processes = n_tempering_levels,
#        walltime="10-00:00:00",
#        n_threads = 24,
#        dependencies = [
#            Pigeons, # <- Pigeons itself can be skipped, added automatically
#           preamble_path
#        ],
#        mpiexec_args=`--mpi=pmi2`
#    )
#)


using MCMCChains, PairPlots, CairoMakie
samples = Pigeons.PT("/n/home06/dochang/2018JuKeBOXFit/results/all/2025-03-07-13-57-26-fU4QsXkC")#results/all/2025-02-11-16-26-25-7dGWgmJ0")#results/all/2025-02-10-11-14-56-mSKJXywd")
#pigeons(samples)
#chain = Chains(samples)
chain = sample_array(ascube(postJB), samples)
chain.sky.p2
#chain = filter(x->x.sky.pa <0.0, chain)
scatter(chain.sky.m_d |> rad2μas)
scatter(chain.sky.spin)
scatter(chain.sky.θo)
scatter(chain.sky.pa * 180/π)

pairplot(chain.sky)
#ess_rhat(chain)
using Accessors

curr = rand(chain)
curr.sky.p1
curr.sky.p2
curr.sky.rpeak
curr.sky.m_d |> rad2μas
curr.sky.pa * 180/π
curr.sky.θs
mdl = skymodel(postJB, curr)    
intensitymap(mdl, grid)|> imageviz
log.(max.(intensitymap(mdl, grid), 1e-5))|> imageviz


#@reset curr.sky.pa = 3.0/180*π
@reset curr.sky.θo = 80.0
@reset curr.sky.θs = 90.0
mdl = skymodel(postJB, curr)    
intensitymap(mdl, grid)|> imageviz
log.(max.(intensitymap(mdl, grid), 1e-5))|> imageviz

begin
    #curr = prior_sample(postJB)
    curr = rand(chain)
    mdl = skymodel(postJB, curr)    
    fig = intensitymap(mdl, grid)|> imageviz
    println(curr.sky.pa * 180/π)
    println(curr.sky.p2)
    println(curr.sky.m_d |> rad2μas)
display(fig)
end
#@reset curr.sky.rpeak = 89.0

mdl = skymodel(postJB, curr)
log.(max.(intensitymap(mdl, grid), 1e-5))|> imageviz

mdl.model.transform[1].α |> rad2μas
mdl.model.model.θo .- 163

#mdl.model.model.scene[1].geometry.opening_angle  * 180/π
#mdl.model.model.scene[2].geometry.opening_angle  * 180/π
##Accessors.@reset mdl.model.model.θo = 27.0
#fnames = fieldnames(typeof(mdl.model.model.scene[1].material))
#args = (getfield.(Ref(mdl.model.model.scene[1].material), fnames))
#args = (args[1]..., args[2]..., args[3:end]...)
##@set mdl.model.model.scene[1].material = Krang.ElectronSynchrotronPowerLawIntensity(args...)
##@reset mdl.model.model.scene[1].material.p1 = 0.5
##@set mdl.model.model.scene[2].material.p1 = 0.5

##@reset mdl.model.model.scene[1].material.R = 1.0
#mdl.model.model.scene[1].material.R
#mdl.model.model.scene[1].material.p2
intensitymap(mdl, grid)|> imageviz
smooth(intensitymap(mdl, grid), μas2rad(15.0)/(2.0*√2.0*log(2.0))) |> imageviz
map(x->x[1], chi2.(Ref(postJB), chain) ) ./ length(dlcamp) |> CairoMakie.scatter
map(x->x[2], chi2.(Ref(postJB), chain) ) ./ length(dcphase) |> CairoMakie.scatter


