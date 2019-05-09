# ppfoodwebs: Empirical pitcher plant inquiline food web models.

## Summary

* Models of are naturally occurring ecosystems with co-evoled
  organisms and well defined boundaries are difficult to quantify.
* The inquiline communities of the purple pitcher plant (*Sarracenia
  purpurea*) provide a food web whose size permits replicated
  sampling of this "microecosystem". 
* This project synthesizes work on the structure and dynamics of the
  pitcher plant microecosystem [@Mouquet, @Ellison, @Gotelli, @Baker,
  @Lau] to generate over 3000 empirical models based on sampling of
  food web organism abundances in the field over the course of several
  field seasons.

## Microecosystem Models

Details about the pitcher plant microecosystem see [@Baker, @Baiser,
@Lau].
* For background on the system see [@Ellison].
* Sampling methods of the inquline abundances are detailed in
  [@Gotelli]. 

## Loading Models in R

To load the pitcher plant inquiline food web models, just run the
following from the `src` directory:

```{r eval = FALSE}
source("load_ppnets.R")
```

This will run the model calculations and produce two objects: 

* *pp.nets*: a list of flow matrices (carbon).
* *pp.info*: sampling details for each model.


## Software Dependencies

* R Packages: *Rgraphviz*, *Hmisc*, *enaR*, *ggplot2*, *igraph*,
  *intergraph*, *reshape*, *statnet*, *sna*.
* Package statnet.common requires update R version

# References

