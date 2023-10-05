## Cosmic Ray Impacts on Turbulence

This repo contains scripts supporting the research published in Bustard and Oh 2022 and Bustard and Oh 2023 probing the interplay between cosmic rays and magnetized turbulence.

#### Instructions
To get all the necessary dependencies, start by running `pip install -r requirements.txt`

Next, check out the Jupyter notebook `CR_Turb_Notebook.ipynb`, which is a mixture of code snippets, figures, and explanations -- essentially an interactive summary of Bustard and Oh 2022, 2023. From there, feel free to explore any other Python scripts! 

#### Other Scripts and Details

`HodgeHelmholtz.py` is useful if you have any simulation volume and want to decompose it into compressive and solenoidal motions. The script will print out the ratio of solenoidal to compressive power, and it will generate spectra for each component. 

To run, there are three required command line options you'll need to set: file_path = the path to your files in HDF5 format or any other format that yt can read; files = the file names; num_cells = number of cells along an axis of the simulation volume (we assume a 3D cubic box with equal number of cells in each direction). Everything else should be taken care of by yt. 

To run: `python3 HodgeHelmholtz.py --file_path "/../files/" --files "cr.out1.004*" --num_cells 512 `

The `utils` directory contains a few functions that come up often, namely for taking Fourier transforms of input quantities to create spectra as a function of scale, for instance to get kinetic energy as a function of scale. Within this directory, `HH_utils` contains some functions that help decompose the simulation volume into solenoidal and compressive components.

`PDF_density_multiplot.py` and `PDF_ecr_multiplot.py` show some examples of using yt to generate volume-weighted probability density functions of gas density and cosmic ray energy, respectively. 

`deviations_lognormal.py` further quantifies the standard deviation of quantities such as density, cosmic ray energy, etc. relative to their mean values in log space -- we work in log space because PDFs in compressive turbulence are lognormal, and quantifications like this of $\delta \rho / \rho$, $\delta v/v$, etc. are common in the literature. Namely, in compressive, hydrodynamic turbulence, there's an equivalence between these quantities, but in the presence of cosmic rays, this equivalence is broken. 

`plotFcr_plotly.py` plots fractions of the turbulent energy input that go into cosmic ray energization, cosmic ray energy loss, and dissipation at the grid scale. This script (which is quite gross) uses Plotly to generate Figure 2 in Bustard and Oh 2023 

`multiplot.py` creates Figure 5 in Bustard and Oh 2022, showing density fluctuations, the misalignment between magnetic field and cosmic ray pressure gradient, and the cosmic ray scale length. 