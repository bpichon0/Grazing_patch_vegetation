# Code for running analyses of grazing effect on global drylands

> [!NOTE]
> Contact: Benoît Pichon, **benoit.pichon0@gmail.com**


This repository contains the code used to perform the analyses for both main text and supplementary information. All the code was made on R (*v4.1.0*).

<p align="center">
    <img src="https://github.com/bpichon0/Grazing_patch_vegetation/blob/master/Fig_overview.jpg" width="800">
</p>


## `Installing R dependencies`

> [!IMPORTANT]  
> To install the packages needed for the analyses and create the folder architecture used to save the data sets, please load the file `Structure_grazing_function.R` using: 

```R
source(Structure_grazing_function.R)
```
> If not all the packages are already installed, install them using the commented code line in L15 of `Structure_grazing_function.R`.

## `Replicating the analyses`

People interested in directly replicating the figures can go to "**Replicating the figures**" since we provide all the generated data.
The main script (`Structure_grazing_main.R`) is organized in independant chunk of codes (Steps 1 to 6) that can be seen easily with Rstudio by pressing *Alt+O*. It allows to replicate all statistical analyses of the paper.
It includes the scripts for the simulations.
Note that all images are available in the Data/Images folder. The image corresponding to each line of the dataset can be found using the ID and Sub_id columns of the dataset.

## `Replicating the figures`

The file `Make_figs.R` plot the figures. As we provide all the data to replicate each figure, everything can be runed without replicating the different analyses. The file is organized in different chunks of code, that can be seen by pressing *Alt+O*.

