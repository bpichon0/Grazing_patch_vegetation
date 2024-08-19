# Code for running analyses of grazing effect on global drylands

> [!NOTE]
> Contact: Benoît Pichon, **benoit.pichon0@gmail.com**


This repository contains the code used to perform the analyses for both main text and supplementary information. All the code was made on R (*v4.1.0*).

<p align="center">
    <img src="https://github.com/bpichon0/Grazing_patch_vegetation/blob/master/Fig.overview.jpg" width="800">
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
The main script (`Structure_grazing_function.R`) is organized in independant chunk of codes that can be seen easily with Rstudio by pressing *Alt+O*. It allows to replicate all statistical analyses of the paper.

### Step 4: Comparison to other models 

Last, to replicate the comparison of our approach with other models, we provide a complete framework through a *bash (.sh)* file in `./Data/Model_confirmation_Guichard` and `./Data/Model_confirmation_Kefi`.
The framework draws the parameters, performs the simulations with the two models, infers the parameters, and estimates the distance to the tipping point by our approach and in the two models (mussel-bed and dryland vegetation). 
Then, to postprocess these simulations, you can run Step 7 of the `ABC_drylands_main.R` file.

## `Replicating the figures`

The file `Make_figs.R` plot the figures. As we provide all the data to replicate each figure, everything can be runed without replicating the different analyses. The file is organized in different chunks of code, that can be seen by pressing *Alt+O*.

