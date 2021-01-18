# Lists of available plots divided in categories.

acf_plot_list <- list("acf", "acf_bar")

areas_plot_list <- list("areas", "areas_ridges")

dens_plot_list <- list("dens", "dens_chains", "dens_overlay")

hex_plot_list <- list("hex")

hist_plot_list <- list("hist", "hist_by_chain")

intervals_plot_list <- list("intervals")

neff_plot_list <- list("neff", "neff_hist")

pairs_plot_list <- list("pairs")

parcoord_plot_list <- list("parcoord")

recover_plot_list <- list("recover_hist", "recover_intervals", "recover_scatter")

rhat_plot_list <- list("rhat", "rhat_hist")

scatter_plot_list <- list("scatter")

trace_plot_list <- list("trace", "trace_highlight")

violin_plot_list <- list("violin")

all_plots_list <- list(
  acf = acf_plot_list,
  areas = areas_plot_list,
  dens = dens_plot_list,
  hex = hex_plot_list,
  hist = hist_plot_list,
  intervals = intervals_plot_list,
  neff = neff_plot_list,
  pairs = pairs_plot_list,
  parcoord = parcoord_plot_list,
  recover = recover_plot_list,
  rhat = rhat_plot_list,
  scatter = scatter_plot_list,
  trace = trace_plot_list,
  violin = violin_plot_list,
  combo = "combo"
)
