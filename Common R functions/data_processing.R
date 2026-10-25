# Data-processing functions from the original graph_functions.R.
# Retained here for this standalone project; the source file is unchanged.
library(dplyr)
library(tidyr)

#' Normalize column names for robust matching
norm_col <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

#' Tolerant numeric coercion
to_numeric_tolerant <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.logical(x)) return(as.numeric(x))
  if (is.factor(x)) x <- as.character(x)

  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x[tolower(x) %in% c("na", "n/a", "null", ".", "-", "—", "–", "NA")] <- NA_character_
  x <- gsub("\u00A0", "", x)
  x <- gsub(",", "", x)
  suppressWarnings(as.numeric(x))
}
#' Normalize numeric or textual treatment identifiers
#'
#' Numeric identifiers are stored in a canonical character form (e.g. 2.0
#' becomes "2"), textual missing-value markers become NA, and genuine text
#' labels are retained. Retaining labels is important for mixed Excel columns:
#' coercing the whole column to numeric would silently discard text-labelled
#' treatment arms.
normalize_treatment_ids <- function(x) {
  x_chr <- trimws(as.character(x))
  missing <- is.na(x) | x_chr == "" |
    tolower(x_chr) %in% c("na", "n/a", "null", ".", "-", "—", "–")
  x_chr[missing] <- NA_character_

  x_num <- suppressWarnings(as.numeric(x_chr))
  numeric_id <- !is.na(x_num)
  if (any(numeric_id & abs(x_num - round(x_num)) >= 1e-8)) {
    bad <- unique(x_chr[numeric_id & abs(x_num - round(x_num)) >= 1e-8])
    stop("Numeric treatment identifiers must be integers: ", paste(bad, collapse = ", "))
  }
  x_chr[numeric_id] <- as.character(as.integer(round(x_num[numeric_id])))
  x_chr
}
#' @title Create a Trial-by-Treatment Edge-Indicator Matrix
#' @description Detects "long" or "wide" format and returns the bipartite
#'   incidence matrix B, where B[i, j] is 1 when trial i includes treatment j.
#' @return A binary matrix whose rows are trials and columns are treatments.

create_edge_list <- function(nma_i, study_col = "StudyID", long_format_trt_col = "t", long_format_n_col = "n") {
  
  # 💡【全局清理空列】
  is_completely_empty <- function(v) {
    v_chr <- trimws(as.character(v))
    all(is.na(v_chr) | v_chr == "" | tolower(v_chr) %in% c("na", "n/a", "null", ".", "-", "—", "–"))
  }
  
  nma_i <- nma_i[, !vapply(nma_i, is_completely_empty, logical(1)), drop = FALSE]
  if (ncol(nma_i) == 0) stop("The entire dataset is empty or only contains NA/missing values.")
  
  cols_raw  <- names(nma_i)
  cols_norm <- norm_col(cols_raw)
  
  # 1. 识别 StudyID
  pref_norm <- norm_col(study_col)
  cand_idx <- which(cols_norm == pref_norm | grepl("^studyid\\d*$", cols_norm) | grepl("^id\\d*$", cols_norm))
  if (length(cand_idx) == 0) stop("StudyID/ID column not found.")
  study_col_resolved <- cols_raw[cand_idx[1]]
  
  # 2. 识别 Format (Long vs Wide)
  long_idx <- which(cols_norm == norm_col(long_format_trt_col))
  wide_t_idx <- which(grepl("^t\\d+$", cols_norm))
  
  trial_arms <- NULL

  # Long format (one row per trial arm)
  if (length(long_idx) > 0) {
    t_col <- cols_raw[long_idx[1]]

    trial_arms <- data.frame(
      trial_index = as.integer(factor(nma_i[[study_col_resolved]], levels = unique(nma_i[[study_col_resolved]]))),
      StudyID = nma_i[[study_col_resolved]],
      treatment = normalize_treatment_ids(nma_i[[t_col]])
    ) %>%
      drop_na(trial_index, treatment)

  # Wide format (one row per trial, with t1, t2, ... columns)
  } else if (length(wide_t_idx) >= 2) {
    wide_cols_t <- cols_raw[wide_t_idx]
    wide_subset <- data.frame(
      trial_index = seq_len(nrow(nma_i)),
      StudyID = nma_i[[study_col_resolved]],
      nma_i[wide_cols_t],
      check.names = FALSE
    )

    # pivot_longer() cannot combine numeric and character treatment columns.
    # Normalize every treatment column first, while preserving genuine labels.
    wide_subset[wide_cols_t] <- lapply(
      wide_subset[wide_cols_t],
      normalize_treatment_ids
    )

    trial_arms <- wide_subset %>%
      pivot_longer(cols = all_of(wide_cols_t), names_to = "arm", values_to = "treatment") %>%
      drop_na(treatment)
  } else {
    stop("Could not determine data format (missing long 't' or wide 't1,t2,...').")
  }

  if (nrow(trial_arms) == 0) return(matrix(integer(0), nrow = 0, ncol = 0))
  trial_arms <- trial_arms %>%
    distinct(trial_index, treatment, .keep_all = TRUE) %>%
    group_by(trial_index) %>%
    filter(n() >= 2) %>%
    ungroup()

  if (nrow(trial_arms) == 0) return(matrix(integer(0), nrow = 0, ncol = 0))

  trial_info <- trial_arms %>% distinct(trial_index, StudyID) %>% arrange(trial_index)
  treatments <- unique(trial_arms$treatment)
  treatment_numbers <- suppressWarnings(as.numeric(treatments))
  treatments <- treatments[order(
    is.na(treatment_numbers),
    ifelse(is.na(treatment_numbers), Inf, treatment_numbers),
    treatments
  )]
  B <- matrix(
    0L,
    nrow = nrow(trial_info),
    ncol = length(treatments),
    dimnames = list(
      make.unique(as.character(trial_info$StudyID)),
      as.character(treatments)
    )
  )
  B[cbind(match(trial_arms$trial_index, trial_info$trial_index),
          match(trial_arms$treatment, treatments))] <- 1L
  attr(B, "StudyID") <- trial_info$StudyID

  return(B)
}
