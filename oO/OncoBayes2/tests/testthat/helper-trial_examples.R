## Examples used for blrm_trial tests

library(tibble)
library(tidyr)
library(dplyr)

examples <- list(

  # Single-agent example -------------------------------------------------------
  single_agent = list(
    histdata = dplyr::bind_rows(list(
      tibble::tibble(
        group_id = "hist_A",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble::tibble(
        group_id = "hist_B",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      )
    )),
    dose_info = dplyr::bind_rows(list(
      tibble::tibble(
        group_id = "cur_A",
        drug1 = 1:2
      ),
      tibble::tibble(
        group_id = "cur_B",
        drug1 = 1:2
      )
    )),
    drug_info = tibble::tibble(
      drug_name = "drug1",
      dose_ref  = 1,
      dose_unit = "ngogn",
      reference_p_dlt = 0.1
    )
  ),

  # Combo2 example -------------------------------------------------------
  combo2 = list(
    histdata = dplyr::bind_rows(list(
      tibble::tibble(
        group_id = "hist_A",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble::tibble(
        group_id = "hist_B",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        num_patients = 3,
        num_toxicities = 0
      )
    )),
    dose_info = dplyr::bind_rows(list(
      tibble::tibble(
        group_id = "cur_A",
        drug1 = 1:2,
        drug2 = 100 * (1:2)
      ),
      tibble::tibble(
        group_id = "cur_B",
        drug1 = 1:2,
        drug2 = 100 * (1:2)
      )
    )),
    drug_info = tibble::tibble(
      drug_name = c("drug1", "drug2"),
      dose_ref  = c(1, 100),
      dose_unit = c("ngogn", "potrzebie"),
      reference_p_dlt=0.3
    )
  ),

  # Combo3 example -------------------------------------------------------
  combo3 = list(
    histdata = dplyr::bind_rows(list(
      tibble::tibble(
        group_id = "hist_A",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        drug3 = 1000 * (1:2),
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble::tibble(
        group_id = "hist_B",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        drug3 = 1000 * (1:2),
        num_patients = 3,
        num_toxicities = 0
      )
    )),
    dose_info = dplyr::bind_rows(list(
      tibble::tibble(
        group_id = "cur_A",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        drug3 = 1000 * (1:2)
      ),
      tibble::tibble(
        group_id = "cur_B",
        drug1 = 1:2,
        drug2 = 100 * (1:2),
        drug3 = 1000 * (1:2)
      )
    )),
    drug_info = tibble::tibble(
      drug_name = paste0("drug", 1:3),
      dose_ref  = 10 ^ c(0, 2, 3),
      dose_unit = c("ngogn", "potrzebie", "blintz"),
      reference_p_dlt = 0.1
    )
  ),

  single_drug_with_strata = list(
    histdata = dplyr::bind_rows(list(
      tibble::tibble(
        group_id = "hist_A",
        stratum_id = "strat_A",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble::tibble(
        group_id = "hist_B",
        stratum_id = "strat_A",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      ),
      tibble::tibble(
        group_id = "hist_C",
        stratum_id = "strat_B",
        drug1 = 1:2,
        num_patients = 3,
        num_toxicities = 0
      )
    )),
    dose_info = dplyr::bind_rows(list(
      tibble::tibble(
        group_id = "cur_A",
        stratum_id = "strat_A",
        drug1 = 1:2
      ),
      tibble::tibble(
        group_id = "cur_B",
        stratum_id = "strat_B",
        drug1 = 1:2
      )
    )),
    drug_info = tibble::tibble(
      drug_name = "drug1",
      dose_ref  = 1,
      dose_unit = "ngogn"
    )
  ),
  
  multi_drug_single_group = list(
    histdata = tibble::tibble(
      group_id = "single_group",
      A = 1,
      B = 1,
      C = 1,
      num_patients = 3,
      num_toxicities = 0
    ),
    dose_info = tibble::tibble(
      group_id = "single_group",
      A = 1,
      B = 1,
      C = 1
    ),
    drug_info = tibble::tibble(
      drug_name = c("A", "B", "C"),
      dose_ref  = 1,
      dose_unit = "ngogn"
    )    
  )
)
