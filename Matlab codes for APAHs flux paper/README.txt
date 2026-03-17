This repository contains the model code and input parameter table used in the manuscript: “A Major Missing Source of Alkylated Polycyclic Aromatic Hydrocarbon Emissions from Oil Spill Volatilization”

The file APAHs15 & oil calculated 3.0 20251023.xlsx contains a summary of all parameters used in the model. Corresponding tables and references can be found in the Supporting Information (SI) of the article.

The model was developed in MATLAB (version 2023a). The calculation workflow is as follows:

1. Oil_APAHs_data_import.m: Used to import all parameters to be calculated.

2. Monte_carlo_main_20260120.m: Main Monte Carlo loop code, used to perform calculations after importing parameters. The number of Monte Carlo runs can be set by the user. Please save the resulting data matrix file from the run for subsequent use by the plotting codes.

3. With the resulting data matrix file loaded, run the subsequent plotting codes to generate figures. Note that the figure numbers do not strictly correspond to those in the article.

Notes: 
- All model assumptions and parameter definitions are documented in the manuscript and its Supporting Information. 
- For questions regarding the code or input files, please contact the corresponding author.