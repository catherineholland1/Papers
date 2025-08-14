# A Bayesian hierarchical mixture model for classifying compositional data with structural zeros #

Code and supplementary material for paper *A Bayesian hierarchical mixture model for classifying compositional data with structural zeros* submitted to Computational Statistics \& Data Analysis.

As we do not have permissions to share the data used for the work, simulated data has been created to mimic the structure of the data used for this work. This is located within the Data folder.

### Abstract

Compositional data takes the form of parts of some whole, consisting of sets of non-negative components, often expressed as proportions, counts or other non-negative values, which are difficult to analyse using traditional statistical techniques. Compositional data analysis has emerged as a powerful tool in forensic glass analysis as it can help account for the dependencies between the chemical elements in glass. However, compositional data can often include many zero values and exhibit a multilevel hierarchical structure. This presence of zeros renders the main technique of using a log-ratio transformation unsuitable. 

To address this, we propose a flexible integrated clustering approach within a Bayesian hierarchical model for compositional data with a large proportion of structural zeros in the compositions and a multilevel hierarchical structure. We apply our methodology to a forensic elemental glass database where the interest lies in classifying the type of glass each fragment belongs to - a common task in forensic data science. Through a out-of-sample classification performance via five-fold cross-validation, we demonstrate improved performance over less flexible alternatives which require manual or expert knowledge. Additionally, we implement our methodology using NIMBLE - a package for which allows for flexible implementation of Bayesian models - which results in computationally efficient frameworks, making it practical for many real-world applications.
