# Found Graph Data and Planted Vertex Covers

This Julia software accompanies the following paper:

- [Found Graph Data and Planted Vertex Covers](http://www.cs.cornell.edu/~arb/papers/fgd-pvc-NIPS-2018.pdf).
  Austin R. Benson and Jon Kleinberg. In *Proceedings of Neural Information Processing (NeurIPS)*, 2018.

This code is designed to reproduce the results in the paper as well as provide some basic code which could be used by others to test out their algorithms on planted vertex covers.

The data used in the paper is also available [here](http://www.cs.cornell.edu/~arb/data/index.html#pvc).

### Setup

Download the software and the email-Enron, email-W3C, call-Reality, and text-Reality datasets.

```bash
git clone https://github.com/arbenson/FGDnPVC
```

The code is written in Julia 1.0. To run everything, you need the following packages:

```julia
using Pkg
Pkg.add("FileIO")
Pkg.add("JLD2")
Pkg.add("LightGraphs")
Pkg.add("PyPlot")
Pkg.add("MathProgBase")
Pkg.add("Gurobi")
Pkg.add("ScikitLearn")
```

Gurobi requires a license, which is free for academics and students. We only use the Gurobi library to find the minimum vertex cover size, so you can still run several parts of the code without it.

### Figure 2: 1-hop neighborhoods covering 2-hop neighborhoods

This code reproduces Figure 2, which shows improved bounds from the maximal matching 2-approximation algorithm for minimum vertex cover.

Note: this function is multithreaded, so you can run with, e.g., JULIA_NUM_THREADS=4 julia.

```julia
include("neighborhoods_minimum_vc_bounds.jl")
# Read graph from data/static/CollegeMsg.txt and 
# creates output file output/CollegeMsg-neighborhood-stats.mat
neighborhoods_minimum_VC_bounds("CollegeMsg")
neighborhoods_minimum_VC_bounds("CollegeMsg", niter=30)  # same but with 30 approximations

# Same thing for Caida and Gnutella data
neighborhoods_minimum_VC_bounds("as-caida20071105")
neighborhoods_minimum_VC_bounds("p2p-Gnutella31")

# Generate plots in Figure 1.
include("paper_plots.jl")
neighborhoods_plot("CollegeMsg")        # --> CollegeMsg-neighborhoods.png
neighborhoods_plot("as-caida20071105")  # --> as-caida20071105-neighborhoods.png
neighborhoods_plot("p2p-Gnutella31")    # --> p2p-Gnutella31-neighborhoods.png
```

### Table 1: Summary statistics of the datasets

This reproduces results in Table 1. Note that you need a Gurobi license in order to use this script, as it uses the Gurobi linear integer programming solver. They provide [free academic licenses](https://user.gurobi.com/download/licenses/free-academic).

```julia
include("summary_statistics.jl")
summary_statistics("email-Enron")
summary_statistics("email-W3C")
summary_statistics("call-Reality")
summary_statistics("text-Reality")
```

### Figure 4: recovery performance experiments

This reproduces results in Figure 4. Again, it is multi-threaded.

```julia
include("temporal_analysis.jl")
# collect data
recovery_over_time("text-Reality")  # --> output/text-Reality-temporal-perf-stats.jld2

# make plot
include("paper_plots.jl")
recovery_plots("text-Reality")  # --> text-Reality-temporal.eps
```

We used separate software for belief propagation and the Path-Core scores. The pre-computed values are stored in `output/dataset-temporal-perf-stats-FULL.jld2`, where dataset is, e.g., email-Enron.

```julia
include("paper_plots.jl")
recovery_plots("email-Enron", full=true)  # --> email-Enron-temporal-FULL.eps
```

### Table 2: timing experiments

This code runs the timing experiments, using Julia's `@time` macro.

```julia
include("timing_analysis.jl")

# (may want to run each command twice for warm-up comparison)
timing_experiment("text-Reality", "degree");
timing_experiment("text-Reality", "UMVC");
timing_experiment("text-Reality", "betweenness");
timing_experiment("text-Reality", "BorgattiEverett");
```