#
#  ccdrAlgorithm-messages.R
#  ccdrAlgorithm
#
#  Created by Bryon Aragam (local) on 11/20/16.
#  Copyright (c) 2014-2017 Bryon Aragam. All rights reserved.
#

#
# PACKAGE CCDRALGORITHM: Messages
#
#   CONTENTS:
#       max_nodes_warning
#

### These warnings are all internal to this package and hence
###  do not need to be exported

### User inputs invalid data object
max_nodes_warning <- function(numnode){
    msg <- "This dataset contains more than %d variables -- in order to
            run CCDr on this dataset, please download the source, increase
            _MAX_CCS_ARRAY_SIZE_ to at least %d, and re-build the package
            from source. If you have any trouble, please contact the
            maintainer."
    stop(sprintf(msg, MAX_CCS_ARRAY_SIZE(), numnode))
}
