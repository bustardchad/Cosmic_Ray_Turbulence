## Cosmic Ray Impacts on Turbulence

This repo contains scripts supporting the research published in Bustard and Oh 2022 and Bustard and Oh 2023 probing the interplay between cosmic rays and magnetized turbulence.

#### Instructions
To get all the necessary dependencies, start by running `pip install -r requirements.txt`

Next, check out the Jupyter notebook `CR_Turb_Notebook.ipynb`, which is a mixture of code snippets, figures, and explanations -- essentially an interactive summary of Bustard and Oh 2022, 2023

From there, feel free to explore any other Python scripts. In particular, HodgeHelmholtz.py is useful if you have any simulation volume and want to decompose it into compressive and solenoidal motions. The script will print out the ratio of solenoidal to compressive power, and it will generate spectra for each component. To run, there are three required command line options you'll need to set: file_path = the path to your files in HDF5 format or any other format that yt can read; files = the file names; num_cells = number of cells along an axis of the simulation volume (we assume a 3D cubic box with equal number of cells in each direction). Everything else should be taken care of by yt. E.g. to run, type 

`python3 HodgeHelmholtz.py --file_path "/../files/" --files "cr.out1.004*" --num_cells 512 `