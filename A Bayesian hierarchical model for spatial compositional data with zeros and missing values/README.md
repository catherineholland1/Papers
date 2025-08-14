# A Bayesian hierarchical model for spatial compositional data with zeros and missing values #

Code and supplementary material for paper *A Bayesian hierarchical model for spatial compositional data with zeros and missing values* submitted to Spatial Statistics.

As we do not have permissions to share the data used for the work, simulated data has been created to mimic the structure of the data used for this work. This is located within the Data folder.

### Abstract

Compositional data –- non-negative parts of some whole –- arise in many spatial contexts across the physical and social sciences. However, modelling such data is challenging due to the need to account for both spatial and compositional structures simultaneously. Reviewing existing approaches for spatial compositional data, usually based on log-ratio transformations, we identify a lack of a general framework that accommodates both zeros and missing components in the compositions, allows for flexible variance structures, and pools spatial information across zero and non-zero data. 

To address this, we present a Bayesian hierarchical approach based on the Generalised-Dirichlet-Multinomial family of distributions. Here, spatial dependence is modelled via flexible latent effects, such as two-dimensional penalised regression splines, supporting a variety of spatial structures and spatially heterogeneous variance patterns. We apply our method to tree species canopy cover data from a UK woodland, where most observations are zeros, such that traditional log-ratio approaches may be unsuitable. Through an out-of-sample prediction experiment, we demonstrate improved performance over simpler alternatives, both in reconstructing missing components within partially observed compositions and in predicting at unseen locations.
