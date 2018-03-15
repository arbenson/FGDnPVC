using MAT
using PyPlot

""" 
recovery_plots
--------------

Makes a plot similar to those appearing in Figure 2 of the paper. This function
relies on the file output/$dataset-temporal-perf-stats.mat, which is
produced with the recovery_over_time() function.

recovery_plots(dataset::String, full::Bool=false)

Input parameters:
- dataset::String: dataset name
- full::Bool: if set to true, then uses pre-computed data which includes
              belief propagation and Path-Core scores.

Writes file $dataset-temporal.eps if full is false or $dataset-temporal-FULL.eps
if full is set to true.
"""
function recovery_plots(dataset::String, full::Bool=false)
    data = matread("output/$dataset-temporal-perf-stats.mat")
    if full
        data = matread("output/$dataset-temporal-perf-stats-FULL.mat")
    end
    ps = data["ps"]
    days = cumsum(data["interval"] * ones(length(ps)))
    fsz = 14
    
    close()
    figure(figsize=(16, 3.75))

    # precision @ core size
    subplot(131)
    plot(days, data["upb_pacs"], color="k",       lw=3.5,  linestyle="-",  label="upper bound")
    plot(days, data["mvc_pacs"], color="#1b9e77", lw=2,    linestyle="-",  label="UMVC")
    plot(days, data["deg_pacs"], color="#d95f02", lw=2,    linestyle=":",  label="Degree")
    plot(days, data["btw_pacs"], color="#7570b3", lw=1.25, linestyle="-",  label="Betweenness")
    if full
        plot(days, data["bpr_pacs"], color="#e6ab02", lw=0.75, linestyle=":",  label="BP")
        plot(days, data["pcs_pacs"], color="#e7298a", lw=1.25, linestyle="--", label="PC")
    end
    plot(days, data["bev_pacs"], color="#66a61e", lw=1.5,  linestyle="-.", label="BE")
    ax = gca()    
    ax[:tick_params]("both", labelsize=fsz-3, length=4, width=1.25)
    if dataset == "email-Enron"
        legend(fontsize=fsz-2, frameon=false, bbox_to_anchor=(0.25, 0.3),
               ncol=2, columnspacing=0.5, labelspacing=0.0)
    end
    xlabel("Days", fontsize=fsz)
    ylabel("Precision at core size", fontsize=fsz)
    title(dataset, fontsize=fsz)

    # area under the precision-recall curve
    subplot(132)
    plot(days, data["upb_auprc"], color="k",       lw=3.5,  linestyle="-",  label="upper bound")
    plot(days, data["mvc_auprc"], color="#1b9e77", lw=2,    linestyle="-",  label="UMVC")
    plot(days, data["deg_auprc"], color="#d95f02", lw=2,    linestyle=":",  label="Degree")
    plot(days, data["btw_auprc"], color="#7570b3", lw=1.25, linestyle="-",  label="Betweenness")
    if full
        plot(days, data["bpr_auprc"], color="#e6ab02", lw=0.75, linestyle=":",  label="BP")
        plot(days, data["pcs_auprc"], color="#e7298a", lw=1.25, linestyle="--", label="PC")
    end
    plot(days, data["bev_auprc"], color="#66a61e", lw=1.5,  linestyle="-.", label="BE")
    ax = gca()    
    ax[:tick_params]("both", labelsize=fsz-3, length=4, width=1.25)
    if dataset == "email-Enron"
        legend(fontsize=fsz-2, frameon=false, bbox_to_anchor=(0.25, 0.3),
               ncol=2, columnspacing=0.5, labelspacing=0.0)
    end
    xlabel("Days", fontsize=fsz)
    ylabel("Area under P-R curve", fontsize=fsz)
    title(dataset, fontsize=fsz)

    # Network statistics
    subplot(133)
    core_present = data["frac_core"] * data["ncore"]
    nodes_present = data["frac_nodes"] * data["nnodes"]
    frac_in_core = core_present ./ nodes_present
    plot(days, data["frac_core"],  color="#e41a1c", lw=2.0, linestyle="-",  label="core nodes")
    plot(days, data["frac_nodes"], color="#377eb8", lw=1.0, linestyle="-",  label="total nodes")
    plot(days, data["frac_edges"], color="#4daf4a", lw=1.5, linestyle=":",  label="total edges")
    plot(days, frac_in_core,       color="#984ea3", lw=1.5, linestyle="--", label="nodes in core")
    if dataset == "email-Enron"
        legend(fontsize=fsz-2, frameon=false, bbox_to_anchor=(0.48, 0.09), labelspacing=0.25)
    end
    ax = gca()    
    ax[:tick_params]("both", labelsize=fsz-3, length=4, width=1.25)
    xlabel("Days", fontsize=fsz)
    ylabel("Fraction", fontsize=fsz)
    title("$dataset", fontsize=fsz)

    tight_layout()
    if   full; savefig("$dataset-temporal-FULL.eps")
    else       savefig("$dataset-temporal.eps")
    end
    show()
end

""" 
neighborhoods_plot
------------------

Makes a plot similar to those appearing in Figure 1 of the paper. This function
relies on the file output/$dataset-neighborhood-stats.mat, which is produced
with the neighborhoods_minimum_VC_bounds() function.

neighborhoods_plot(dataset::String)

Input parameters:
- dataset::String: dataset name

Produces file $dataset-neighborhoods.png.
"""
function neighborhoods_plot(dataset::String)
    close()
    data = matread("output/$dataset-neighborhood-stats.mat")
    counted = data["counted"]
    keep = counted .== 1
    N1 = data["N1"][keep]
    N2 = data["N2"][keep]

    u = data["u"][keep]
    l = data["l"][keep]
    improved_bounds = max.(min.((N1 - u + 2) .* l, (N1 + 1).^2 / 4 + N1), N1)
    x = collect(minimum(N1):maximum(N1))
    x_ub = (x + 1).^2 / 4 + x

    fsz=18
    plot(x, x_ub, color="k", lw=1.5, label="(k+1)^2 / 4 + k")
    scatter(N1, N2,              marker=".", alpha=0.2, s=7, label="node")
    scatter(N1, improved_bounds, marker="s", alpha=0.2, s=5, label="improved bound")
    plot(x, x, color="k", lw=0.75, label="k")
    n = length(counted)
    plot(x, n * ones(length(x)), color="k", lw=1.5, linestyle="--", label="# of nodes")
    ax = gca()
    ax[:set_xscale]("log")
    ax[:set_yscale]("log")
    xlabel("1-hop neighborhood size (k)", fontsize=fsz)
    ylabel("2-hop neighborhood size", fontsize=fsz)
    title(dataset, fontsize=fsz)
    ax = gca()    
    ax[:tick_params]("both", labelsize=fsz, length=6, width=2)
    ax[:tick_params]("both", which="minor", length=3, width=1.25)    
    if dataset == "CollegeMsg"
        leg = legend(fontsize=fsz-3, frameon=false, bbox_to_anchor=(0.5, 0.37),
                     labelspacing=0.0, markerscale=4)
    end
    tight_layout()
    savefig("$dataset-neighborhoods.png", dpi=500)
    show()
end
