mkdir output/figs/manuscript/$(date +%F)  # create a folder with today's date
cp -a output/figs/adjusted_diversity_indices.png output/figs/manuscript/$(date +%F)/Figure_1.png
cp -a output/figs/adjusted_SD_C34_abund.png output/figs/manuscript/$(date +%F)/Figure_2.png
cp -a output/figs/C43ratio_NP_partial_regression_plot.png output/figs/manuscript/$(date +%F)/Figure_3.png
cp -a output/figs/daily_env_var.png output/figs/manuscript/$(date +%F)/Figure_S1.png
cp -a output/figs/fig_relativeabund_dom_sub.png output/figs/manuscript/$(date +%F)/Figure_S2.png
cp -a output/figs/FACE_vegetation_CO2_Scatter.png output/figs/manuscript/$(date +%F)/Figure_S3.png
cp -a output/figs/adjusted_c43_ratio.png output/figs/manuscript/$(date +%F)/Figure_S4a.png
cp -a output/figs/adjusted_sd_ratio.png output/figs/manuscript/$(date +%F)/Figure_S4b.png
cp -a output/figs/adjusted_diversity_indices.pdf output/figs/manuscript/$(date +%F)/Figure_1.pdf
cp -a output/figs/adjusted_SD_C34_abund.pdf output/figs/manuscript/$(date +%F)/Figure_2.pdf
cp -a output/figs/C43ratio_NP_partial_regression_plot.pdf output/figs/manuscript/$(date +%F)/Figure_3.pdf
