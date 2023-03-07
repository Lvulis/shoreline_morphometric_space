# shoreline_morphometric_space
This repository contains the code to reproduce key elements of the shoreline analysis done in Vulis et al., (2023). A single example delta is given to show necessary file structure.

The following workflow is used to reproduce the results:
- Water occurrence maps for each delta are first obtained from the Global Surface Water dataset via Google Earth Engine using `GSW_monthly_history.R` or `GSW_occurrence.R`. A shapefile containing delta outlines is in the "./analysis/globaldata/delta_outlines" folder.
- Water masks and shorelines are created using `maskmaking.R`. This is meant to be run in semi-interactive mode as it's necessary to choose a cardinal direction when walking on the delta. This uses the ROAM package available on [https://github.com/lvulis/ROAM](GitHub) which will have to be installed from source.
- Run `main.R` to process the shorelines into spatial-series, compute relevant metrics, and plot figures. 