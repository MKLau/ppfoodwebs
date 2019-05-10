ppfoodwebs: Empirical pitcher plant inquiline food web models
=============================================================

-   Models of are naturally occurring ecosystems with co-evoled
    organisms and well defined boundaries are difficult to quantify.
-   The inquiline communities of the purple pitcher plant (*Sarracenia
    purpurea*) provide a food web whose size permits replicated sampling
    of this "microecosystem".
-   This project synthesizes work on the structure and dynamics of the
    pitcher plant microecosystem to generate over 3000 empirical models
    based on sampling of food web organism abundances in the field over
    the course of several field seasons.

Microecosystem Models
=====================

Details about the pitcher plant microecosystem see (Butler, Gotelli, and
Ellison 2008, Ellison and Gotelli (2009), Baiser, Ardeshiri, and Ellison
(2011)).

Loading Models in R
===================

To load the pitcher plant inquiline food web models, just run the
following from the `src` directory:

    source("load_ppnets.R")

This will run the model calculations and produce two objects:

-   *pp.nets*: a list of flow matrices (currently carbon only).
-   *pp.info*: sampling details for each model.

Software Dependencies
=====================

-   R Packages: *Rgraphviz*, *Hmisc*, *enaR*, *ggplot2*, *igraph*,
    *intergraph*, *reshape*, *statnet*, *sna*.
-   Package statnet.common requires update R version

References
==========

Baiser, Benjamin, Roxanne S. Ardeshiri, and Aaron M. Ellison. 2011.
“Species richness and trophic diversity increase decomposition in a
co-evolved food web.” Edited by Simon Thrush. *PLoS One* 6 (5). Public
Library of Science: e20672.
doi:[10.1371/journal.pone.0020672](https://doi.org/10.1371/journal.pone.0020672).

Butler, J L, N J Gotelli, and A M Ellison. 2008. “Linking the brown and
green: nutrient transformation and fate in the ∖emph{Sarracenia}
microecosystem.” *Ecology* 89: 898–904.

Ellison, A M, and N J Gotelli. 2009. “Energetics and the evolution of
carnivorous plants - Darwin’s ‘most wonderful plants in the world’.”
*Journal of Experimental Botany* 60: 19–42.
