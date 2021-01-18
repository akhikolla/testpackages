# BeviMed 5.1

* New function, `subset_variants`, which retains only variants with data bearing upon pathogenicity.
* Return posterior mean of `omega` even when not explicitly sampled in `summary.BeviMed_m`.

# BeviMed 5.0

* `bevimed_polytomous` function added which enables application of BeviMed across multiple association models. 
* `BeviMed` objects now more general, representing results of inference with respect to the baseline model `gamma = 0` and an arbitrary number of alternative association models - typically, one for each mode of inheritance. The `$moi` slot has been replaced with `$models`.
* `prob_pathogenic` now returns a list when broken down by mode of inheritance/model.

# BeviMed 4.3

* Make BeviMed work smoothly when number of individuals or number of variants is 0.
* Retain names of variants from columns of original allele count matrix.
* Improvements to guide, with more detail on model selection.

# BeviMed 4.2

* Fixed bug in calculation of expected number of explaining variants by only including those with pathogenic configurations.

# BeviMed 4.0

* Previous `bevimed` function now replaced by `bevimed_m`, with the `_m` indicating that it conditions on mode of inheritance. 
* `bevimed` now integrates over indicator of association (gamma) and mode of inheritance (m), allowing user to specify priors on probability of association and probability of dominance.
* The `BeviMed` class object has been replaced by `BeviMed_m`, and a new `BeviMed` class has been introduced for inference with respect to all models: gamma 0 and gamma 1 under each mode of inheritance.
* A new vignette with more detail called `BeviMed Guide` which relates the package to the paper.
* Names used for summary statistics in summary objects have changed, see function help pages for details on current names. 
* `print`ing a `BeviMed` object now shows conditional probabilities of pathogenicity for each mode of inheritance, and expected explained cases and expected explaining variants shown too.
* Bug fixed in adaptive tuning for omega and phi proposals.

# BeviMed 3.0

* Re-naming of parameters in `bevimed` function to match the names of variables in the paper (under submission). 
* The allele count matrix `G` should now be supplied as a matrix with rows corresponding to individuals, not variants.
* `expected_explained` and `explaining_variants` functions have been added, respectively computing the expected number of cases with their disease explained by the given variants, and expected number of pathogenic variants present amongst cases.
