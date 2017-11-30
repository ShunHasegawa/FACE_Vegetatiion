mkdir output/figs/manuscript/$(date +%F)  # create a folder with today's date
cp -a output/figs/adjusted_diversity_indices.png output/figs/manuscript/$(date +%F)/fig1.png
cp -a output/figs/adjusted_SD_C34_abund.png output/figs/manuscript/$(date +%F)/fig2.png
cp -a output/figs/C43ratio_NP_partial_regression_plot.png output/figs/manuscript/$(date +%F)/fig3.png
cp -a output/figs/daily_env_var.png output/figs/manuscript/$(date +%F)/figS1.png
cp -a output/figs/fig_relativeabund_dom_sub.png output/figs/manuscript/$(date +%F)/figS2.png
cp -a output/figs/FACE_vegetation_CO2_Scatter.png output/figs/manuscript/$(date +%F)/figS3.png
cp -a output/figs/adjusted_c43_ratio.png output/figs/manuscript/$(date +%F)/figS4a.png
cp -a output/figs/adjusted_sd_ratio.png output/figs/manuscript/$(date +%F)/figS4b.png
