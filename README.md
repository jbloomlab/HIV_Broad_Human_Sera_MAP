# HIV Broad Human Sera MAP README

This repository converts an analysis of MAP data using broad human anti-HIV sera and HIV Envelope strain BG505, performed and written by Adam Dingens, into interactive plots. This repository was made by Caelan Radford. 

## Organization

The `adingens_analysis` folder contains the analysis written by Adam Dingens. Caelan did not modify its contents. This folder contains its own documentation, which should be viewed for details on that analysis. 

The `generate_interactive_plots.ipynb` converts the results from Adam Dingens' analysis into interactive plots from the [*polyclonal*](https://jbloomlab.github.io/polyclonal/) software. Because these plots are produced outside of a pipeline like the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline), they will not contain as much information as plots produced by that pipeline, such as [this](https://dms-vep.github.io/HIV_Envelope_BF520_DMS_CD4bs_sera/IDC508_escape_plot.html) one. 

The `interactive_plots` folder contains the interactive plots produced in html files. 

The `results_from_lentivirus_dms` folder contains interactive plots produced from lentivirus deep mutational scanning using the BF520 strain for some of these same sera. You can view the pre-print describing this work [here](https://www.biorxiv.org/content/10.1101/2023.03.23.533993v1.abstract).