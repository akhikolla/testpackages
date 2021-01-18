# sboost 0.1.1

## Major Changes

* Classifier output from *sboost()* now includes *right_categories* column which is similar to *left_categories* but is associated with the outcomes in *right* column. When assessing this new classifier, if a categorical input cannot be found in either *right_categories* or *left_categories* (i.e. was not found in training data) the vote for this feature will now be 0. (Before this, if an input was not found in *left_categories*, it was assumed that the input would be associated with the *right* outcome.)

* There is a new optional parameter in the *sboost()* and *validate()* functions called *verbose*. The default value for *verbose* is FALSE, and there is no change from previous versions when *verbose* = FALSE. If *verbose* is set to TRUE, a progress bar will appear in the console for each classifier that is created.

## Minor Changes
* The Description of the package in the DESCRIPTION file now contains a reference to Freund and Schapire's paper on AdaBoost.