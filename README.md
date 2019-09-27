# `OncoCast` R package
Ensemble learner for survival prediction and stratification with improved interface. Allows for multiple machine learning to be applied to a high-dimensional data set with survival outcome (including left-truncated data).

## Installing `OncoCast`

`OncoCast` has some dependencies that will be installed if they are not found in your library. When installing from github a prompt in the console may ask if you want to install the binaries for curl v4.0 instead of 3.3. There is no need to update it for the package to work properly.

```{r}
install.packages("devtools")
devtools::install_github("AxelitoMartin/OncoCast",build_vignettes = T, build_manual=T)
```
The vignette will need a minute to build when installing the package.

## Using `OncoCast` 

We recommend first time users to go through the tutorial which can be found using the command:

```{r}
browseVignettes("OncoCast")
```
Once the package has been installed. This will open your default web browser, select HTML and the tutorial will appear.
