# Joint (M, N, S) fixed-edge rewiring simulation functions.


sample_fixed_edges <- function(M, N, total_edges) {
  if (total_edges < 0L || total_edges > M * N) {
    stop("total_edges must be between 0 and M * N.")
  }
  B <- matrix(0L, nrow = M, ncol = N)
  if (total_edges > 0L) {
    B[sample.int(M * N, total_edges, replace = FALSE)] <- 1L
  }
  B
}

is_valid_rewired_network <- function(B) {
  if (any(rowSums(B) < 2L) || any(colSums(B) < 1L)) return(FALSE)
  igraph::is_connected(
    igraph::graph_from_biadjacency_matrix(B, directed = FALSE)
  )
}

fix_trial_degrees <- function(B) {
  for (target_row in sample(seq_len(nrow(B)))) {
    while (sum(B[target_row, ]) < 2L) {
      donor_rows <- which(rowSums(B) > 2L)
      movable_edges <- do.call(rbind, lapply(donor_rows, function(donor_row) {
        columns <- which(B[donor_row, ] == 1L & B[target_row, ] == 0L)
        if (length(columns) == 0L) return(NULL)
        cbind(row = donor_row, column = columns)
      }))
      if (is.null(movable_edges) || nrow(movable_edges) == 0L) return(NULL)
      edge <- movable_edges[sample.int(nrow(movable_edges), 1L), ]
      B[edge["row"], edge["column"]] <- 0L
      B[target_row, edge["column"]] <- 1L
    }
  }
  B
}

fix_treatment_degrees <- function(B) {
  for (target_column in sample(seq_len(ncol(B)))) {
    if (sum(B[, target_column]) >= 1L) next
    donor_columns <- which(colSums(B) > 1L)
    movable_edges <- do.call(rbind, lapply(donor_columns, function(donor_column) {
      rows <- which(B[, donor_column] == 1L)
      if (length(rows) == 0L) return(NULL)
      cbind(row = rows, column = donor_column)
    }))
    if (is.null(movable_edges) || nrow(movable_edges) == 0L) return(NULL)
    edge <- movable_edges[sample.int(nrow(movable_edges), 1L), ]
    B[edge["row"], edge["column"]] <- 0L
    B[edge["row"], target_column] <- 1L
  }
  B
}

connect_rewired_components <- function(B) {
  M <- nrow(B)
  for (step in seq_len(M + ncol(B) - 1L)) {
    graph <- igraph::graph_from_biadjacency_matrix(B, directed = FALSE)
    components <- igraph::components(graph)
    if (components$no == 1L) return(B)
    edge_ends <- igraph::ends(graph, igraph::E(graph), names = FALSE)
    edge_components <- components$membership[edge_ends[, 1L]]
    nonbridge_edges <- setdiff(
      seq_len(igraph::ecount(graph)), as.integer(igraph::bridges(graph))
    )
    if (length(nonbridge_edges) == 0L) return(NULL)
    first_edge <- sample(nonbridge_edges, 1L)
    other_edges <- which(edge_components != edge_components[first_edge])
    if (length(other_edges) == 0L) return(NULL)
    second_edge <- sample(other_edges, 1L)

    orient_edge <- function(vertices) {
      if (vertices[1L] <= M) {
        c(row = vertices[1L], column = vertices[2L] - M)
      } else {
        c(row = vertices[2L], column = vertices[1L] - M)
      }
    }
    edge_one <- orient_edge(edge_ends[first_edge, ])
    edge_two <- orient_edge(edge_ends[second_edge, ])
    B[edge_one["row"], edge_one["column"]] <- 0L
    B[edge_two["row"], edge_two["column"]] <- 0L
    B[edge_one["row"], edge_two["column"]] <- 1L
    B[edge_two["row"], edge_one["column"]] <- 1L
  }
  NULL
}

rewire_until_valid <- function(B, max_attempts = 50L) {
  original_edges <- sum(B)
  if (is_valid_rewired_network(B)) {
    return(list(B = B, attempts = 0L, initially_valid = TRUE))
  }
  for (attempt in seq_len(max_attempts)) {
    candidate <- fix_trial_degrees(B)
    if (is.null(candidate)) next
    candidate <- fix_treatment_degrees(candidate)
    if (is.null(candidate)) next
    candidate <- connect_rewired_components(candidate)
    if (is.null(candidate)) next
    if (sum(candidate) == original_edges && is_valid_rewired_network(candidate)) {
      return(list(B = candidate, attempts = attempt, initially_valid = FALSE))
    }
  }
  NULL
}

fit_MNS_copula <- function(empirical_MNS) {
  required <- c("M", "N", "S")
  if (!all(required %in% names(empirical_MNS))) {
    stop("empirical_MNS must contain M, N and S.")
  }
  empirical_MNS <- empirical_MNS[, required, drop = FALSE]
  if (any(!is.finite(as.matrix(empirical_MNS)))) {
    stop("empirical_MNS must contain only finite M, N and S values.")
  }

  pseudo_observations <- copula::pobs(as.matrix(empirical_MNS))
  fitted <- copula::fitCopula(
    copula::normalCopula(dim = 3L, dispstr = "un"),
    data = pseudo_observations,
    method = "mpl"
  )

  list(
    fit = fitted,
    margins = lapply(empirical_MNS, as.integer),
    coefficients = stats::coef(fitted),
    family = "Three-dimensional Gaussian copula with empirical margins"
  )
}

sample_MNS_copula <- function(copula_model, n = 1L) {
  if (is.null(copula_model$fit) || is.null(copula_model$margins)) {
    stop("copula_model must be returned by fit_MNS_copula().")
  }
  simulated_uniforms <- copula::rCopula(
    n,
    methods::slot(copula_model$fit, "copula")
  )
  empirical_quantile <- function(values, probabilities) {
    as.integer(stats::quantile(
      values, probs = probabilities, type = 1L, names = FALSE
    ))
  }
  data.frame(
    M = empirical_quantile(copula_model$margins$M, simulated_uniforms[, 1L]),
    N = empirical_quantile(copula_model$margins$N, simulated_uniforms[, 2L]),
    S = empirical_quantile(copula_model$margins$S, simulated_uniforms[, 3L])
  )
}

check_failure_reasons <- function(B) {
  M <- nrow(B)
  S <- sum(B)
  disconnected <- S == 0L || !igraph::is_connected(
    igraph::graph_from_biadjacency_matrix(B, directed = FALSE)
  )
  failures <- c(
    disconnected = disconnected,
    S_lt_2M = S < 2L * M,
    trial_degree_lt_2 = any(rowSums(B) < 2L),
    treatment_degree_lt_1 = any(colSums(B) < 1L)
  )
  list(valid = !any(failures), failures = failures)
}

try_one_MNS <- function(M, N, S, allocation_attempts = 10L,
                        rewiring_attempts = 50L) {
  failure_counts <- c(
    disconnected = 0L, S_lt_2M = 0L,
    trial_degree_lt_2 = 0L, treatment_degree_lt_1 = 0L
  )
  for (allocation_attempt in seq_len(allocation_attempts)) {
    B <- sample_fixed_edges(M, N, S)
    check <- check_failure_reasons(B)
    if (check$valid) {
      return(list(accepted = TRUE, B = B,
                  allocation_attempt = allocation_attempt,
                  rewiring_attempts = 0L, failure_counts = failure_counts))
    }
    failure_counts <- failure_counts + as.integer(check$failures)
    rewired <- rewire_until_valid(B, max_attempts = rewiring_attempts)
    if (!is.null(rewired)) {
      return(list(accepted = TRUE, B = rewired$B,
                  allocation_attempt = allocation_attempt,
                  rewiring_attempts = rewired$attempts,
                  failure_counts = failure_counts))
    }
  }
  list(accepted = FALSE, B = NULL, allocation_attempt = allocation_attempts,
       rewiring_attempts = rewiring_attempts, failure_counts = failure_counts)
}

simulate_copula_MNS_rewiring <- function(
    copula_model, n_qualified = 300L,
    allocation_attempts_per_MNS = 10L, max_rewiring_attempts = 50L,
    max_MNS_attempts = 100000L, seed = 123L) {
  set.seed(seed)
  networks <- vector("list", n_qualified)
  summary_rows <- vector("list", n_qualified)
  failure_totals <- c(
    disconnected = 0L, S_lt_2M = 0L,
    trial_degree_lt_2 = 0L, treatment_degree_lt_1 = 0L
  )
  tuple_condition_totals <- c(
    M_lt_2 = 0L, N_lt_4 = 0L, S_lt_2M = 0L, S_lt_N = 0L,
    S_lt_M_plus_N_minus_1 = 0L, S_gt_MN = 0L
  )
  rejected_tuple_condition_totals <- tuple_condition_totals
  allocation_attempts_total <- 0L
  accepted <- 0L
  MNS_attempt <- 0L
  infeasible_MNS <- 0L
  generation_attempts_exhausted <- 0L
  while (accepted < n_qualified) {
    MNS_attempt <- MNS_attempt + 1L
    if (MNS_attempt > max_MNS_attempts) {
      stop("Reached max_MNS_attempts after accepting ", accepted, " networks.")
    }
    parameters <- sample_MNS_copula(copula_model, n = 1L)
    M <- as.integer(parameters$M)
    N <- as.integer(parameters$N)
    S <- as.integer(parameters$S)
    tuple_conditions <- c(
      M_lt_2 = M < 2L,
      N_lt_4 = N < 4L,
      S_lt_2M = S < 2L * M,
      S_lt_N = S < N,
      S_lt_M_plus_N_minus_1 = S < M + N - 1L,
      S_gt_MN = S > M * N
    )
    tuple_condition_totals <-
      tuple_condition_totals + as.integer(tuple_conditions)
    minimum_edges <- max(2L * M, N, M + N - 1L)
    if (M < 1L || N < 1L || S < minimum_edges || S > M * N) {
      infeasible_MNS <- infeasible_MNS + 1L
      rejected_tuple_condition_totals <-
        rejected_tuple_condition_totals + as.integer(tuple_conditions)
      next
    }
    result <- try_one_MNS(
      M, N, S, allocation_attempts_per_MNS,
      max_rewiring_attempts
    )
    allocation_attempts_total <-
      allocation_attempts_total + result$allocation_attempt
    failure_totals <- failure_totals + result$failure_counts
    if (!result$accepted) {
      generation_attempts_exhausted <- generation_attempts_exhausted + 1L
      rejected_tuple_condition_totals <-
        rejected_tuple_condition_totals + as.integer(tuple_conditions)
      next
    }
    accepted <- accepted + 1L
    networks[[accepted]] <- result$B
    summary_rows[[accepted]] <- data.frame(
      simulation_id = accepted, MNS_attempt = MNS_attempt,
      M = nrow(result$B), N = ncol(result$B), S = sum(result$B),
      density = sum(result$B) / length(result$B),
      allocation_attempt = result$allocation_attempt,
      rewiring_attempts = result$rewiring_attempts
    )
  }
  names(networks) <- paste0("copula_MNS_", seq_len(n_qualified))
  final_failure_totals <- Reduce(
    `+`, lapply(networks, function(B) {
      as.integer(check_failure_reasons(B)$failures)
    })
  )
  names(final_failure_totals) <- names(failure_totals)
  list(
    networks = networks,
    summary = do.call(rbind, summary_rows),
    failure_summary = data.frame(
      reason = c(
        names(failure_totals), "infeasible_copula_MNS",
        "generation_attempts_exhausted"
      ),
      count = c(
        as.integer(failure_totals), infeasible_MNS,
        generation_attempts_exhausted
      )
    ),
    rejection_summary = data.frame(
      reason = c(
        "structurally_infeasible_MNS", "generation_attempts_exhausted"
      ),
      count = c(infeasible_MNS, generation_attempts_exhausted)
    ),
    tuple_condition_summary = data.frame(
      condition = names(tuple_condition_totals),
      all_sampled_MNS = as.integer(tuple_condition_totals),
      rejected_MNS = as.integer(rejected_tuple_condition_totals)
    ),
    allocation_failure_summary = data.frame(
      reason = names(failure_totals),
      count = as.integer(failure_totals)
    ),
    final_network_failure_summary = data.frame(
      reason = names(final_failure_totals),
      count = as.integer(final_failure_totals)
    ),
    allocation_attempts_total = allocation_attempts_total,
    MNS_attempts = MNS_attempt,
    rejected_MNS_attempts = MNS_attempt - n_qualified,
    n_qualified = n_qualified,
    simulation_schema_version = 3L,
    copula_family = copula_model$family,
    copula_coefficients = copula_model$coefficients
  )
}

# Tracked variants used by the separate rewiring-effort simulation. One edge
# is counted whenever an existing edge is moved to a new matrix position. A
# component-joining double-edge swap is one rewiring step and two rewired
# edges. Counts describe the successful repair path only; discarded repair
# attempts are not included.

fix_trial_degrees_tracked <- function(B) {
  edges_rewired <- 0L
  for (target_row in sample(seq_len(nrow(B)))) {
    while (sum(B[target_row, ]) < 2L) {
      donor_rows <- which(rowSums(B) > 2L)
      movable_edges <- do.call(rbind, lapply(donor_rows, function(donor_row) {
        columns <- which(B[donor_row, ] == 1L & B[target_row, ] == 0L)
        if (length(columns) == 0L) return(NULL)
        cbind(row = donor_row, column = columns)
      }))
      if (is.null(movable_edges) || nrow(movable_edges) == 0L) return(NULL)
      edge <- movable_edges[sample.int(nrow(movable_edges), 1L), ]
      B[edge["row"], edge["column"]] <- 0L
      B[target_row, edge["column"]] <- 1L
      edges_rewired <- edges_rewired + 1L
    }
  }
  list(B = B, edges_rewired = edges_rewired)
}

fix_treatment_degrees_tracked <- function(B) {
  edges_rewired <- 0L
  for (target_column in sample(seq_len(ncol(B)))) {
    if (sum(B[, target_column]) >= 1L) next
    donor_columns <- which(colSums(B) > 1L)
    movable_edges <- do.call(rbind, lapply(donor_columns, function(donor_column) {
      rows <- which(B[, donor_column] == 1L)
      if (length(rows) == 0L) return(NULL)
      cbind(row = rows, column = donor_column)
    }))
    if (is.null(movable_edges) || nrow(movable_edges) == 0L) return(NULL)
    edge <- movable_edges[sample.int(nrow(movable_edges), 1L), ]
    B[edge["row"], edge["column"]] <- 0L
    B[edge["row"], target_column] <- 1L
    edges_rewired <- edges_rewired + 1L
  }
  list(B = B, edges_rewired = edges_rewired)
}

connect_rewired_components_tracked <- function(B) {
  M <- nrow(B)
  edges_rewired <- 0L
  rewiring_steps <- 0L
  for (step in seq_len(M + ncol(B) - 1L)) {
    graph <- igraph::graph_from_biadjacency_matrix(B, directed = FALSE)
    components <- igraph::components(graph)
    if (components$no == 1L) {
      return(list(
        B = B,
        edges_rewired = edges_rewired,
        rewiring_steps = rewiring_steps
      ))
    }
    edge_ends <- igraph::ends(graph, igraph::E(graph), names = FALSE)
    edge_components <- components$membership[edge_ends[, 1L]]
    nonbridge_edges <- setdiff(
      seq_len(igraph::ecount(graph)), as.integer(igraph::bridges(graph))
    )
    if (length(nonbridge_edges) == 0L) return(NULL)
    first_edge <- sample(nonbridge_edges, 1L)
    other_edges <- which(edge_components != edge_components[first_edge])
    if (length(other_edges) == 0L) return(NULL)
    second_edge <- sample(other_edges, 1L)

    orient_edge <- function(vertices) {
      if (vertices[1L] <= M) {
        c(row = vertices[1L], column = vertices[2L] - M)
      } else {
        c(row = vertices[2L], column = vertices[1L] - M)
      }
    }
    edge_one <- orient_edge(edge_ends[first_edge, ])
    edge_two <- orient_edge(edge_ends[second_edge, ])
    B[edge_one["row"], edge_one["column"]] <- 0L
    B[edge_two["row"], edge_two["column"]] <- 0L
    B[edge_one["row"], edge_two["column"]] <- 1L
    B[edge_two["row"], edge_one["column"]] <- 1L
    edges_rewired <- edges_rewired + 2L
    rewiring_steps <- rewiring_steps + 1L
  }
  NULL
}

rewire_until_valid_tracked <- function(B, max_attempts = 50L) {
  original_edges <- sum(B)
  if (is_valid_rewired_network(B)) {
    return(list(
      B = B, attempts = 0L, initially_valid = TRUE,
      rewiring_steps = 0L, edges_rewired_total = 0L,
      edges_rewired_fix_trial = 0L,
      edges_rewired_fix_treatment = 0L,
      edges_rewired_fix_connected = 0L
    ))
  }
  for (attempt in seq_len(max_attempts)) {
    trial_repair <- fix_trial_degrees_tracked(B)
    if (is.null(trial_repair)) next
    treatment_repair <- fix_treatment_degrees_tracked(trial_repair$B)
    if (is.null(treatment_repair)) next
    connection_repair <- connect_rewired_components_tracked(
      treatment_repair$B
    )
    if (is.null(connection_repair)) next
    candidate <- connection_repair$B
    if (sum(candidate) == original_edges && is_valid_rewired_network(candidate)) {
      trial_edges <- trial_repair$edges_rewired
      treatment_edges <- treatment_repair$edges_rewired
      connected_edges <- connection_repair$edges_rewired
      return(list(
        B = candidate, attempts = attempt, initially_valid = FALSE,
        rewiring_steps = trial_edges + treatment_edges +
          connection_repair$rewiring_steps,
        edges_rewired_total = trial_edges + treatment_edges + connected_edges,
        edges_rewired_fix_trial = trial_edges,
        edges_rewired_fix_treatment = treatment_edges,
        edges_rewired_fix_connected = connected_edges
      ))
    }
  }
  NULL
}

try_one_MNS_tracked <- function(M, N, S, allocation_attempts = 10L,
                                rewiring_attempts = 50L) {
  for (allocation_attempt in seq_len(allocation_attempts)) {
    B <- sample_fixed_edges(M, N, S)
    if (check_failure_reasons(B)$valid) {
      return(list(
        accepted = TRUE, B = B,
        allocation_attempt = allocation_attempt,
        rewiring_attempts = 0L, rewiring_steps = 0L,
        edges_rewired_total = 0L,
        edges_rewired_fix_trial = 0L,
        edges_rewired_fix_treatment = 0L,
        edges_rewired_fix_connected = 0L
      ))
    }
    rewired <- rewire_until_valid_tracked(
      B, max_attempts = rewiring_attempts
    )
    if (!is.null(rewired)) {
      return(c(
        list(
          accepted = TRUE, B = rewired$B,
          allocation_attempt = allocation_attempt
        ),
        rewired[c(
          "attempts", "rewiring_steps", "edges_rewired_total",
          "edges_rewired_fix_trial", "edges_rewired_fix_treatment",
          "edges_rewired_fix_connected"
        )] |>
          setNames(c(
            "rewiring_attempts", "rewiring_steps", "edges_rewired_total",
            "edges_rewired_fix_trial", "edges_rewired_fix_treatment",
            "edges_rewired_fix_connected"
          ))
      ))
    }
  }
  list(accepted = FALSE, B = NULL)
}

simulate_copula_MNS_rewiring_tracked <- function(
    copula_model, n_qualified = 10000L,
    allocation_attempts_per_MNS = 10L, max_rewiring_attempts = 50L,
    max_MNS_attempts = 200000L, seed = 124L) {
  set.seed(seed)
  networks <- vector("list", n_qualified)
  summary_rows <- vector("list", n_qualified)
  accepted <- 0L
  MNS_attempt <- 0L
  structurally_infeasible <- 0L
  generation_attempts_exhausted <- 0L
  while (accepted < n_qualified) {
    MNS_attempt <- MNS_attempt + 1L
    if (MNS_attempt > max_MNS_attempts) {
      stop("Reached max_MNS_attempts after accepting ", accepted, " networks.")
    }
    parameters <- sample_MNS_copula(copula_model, n = 1L)
    M <- as.integer(parameters$M)
    N <- as.integer(parameters$N)
    S <- as.integer(parameters$S)
    minimum_edges <- max(2L * M, N, M + N - 1L)
    if (M < 1L || N < 1L || S < minimum_edges || S > M * N) {
      structurally_infeasible <- structurally_infeasible + 1L
      next
    }
    result <- try_one_MNS_tracked(
      M, N, S,
      allocation_attempts = allocation_attempts_per_MNS,
      rewiring_attempts = max_rewiring_attempts
    )
    if (!result$accepted) {
      generation_attempts_exhausted <- generation_attempts_exhausted + 1L
      next
    }
    accepted <- accepted + 1L
    networks[[accepted]] <- result$B
    summary_rows[[accepted]] <- data.frame(
      simulation_id = accepted,
      MNS_attempt = MNS_attempt,
      M = nrow(result$B), N = ncol(result$B), S = sum(result$B),
      density = sum(result$B) / length(result$B),
      allocation_attempt = result$allocation_attempt,
      rewiring_attempts = result$rewiring_attempts,
      rewiring_steps = result$rewiring_steps,
      edges_rewired_total = result$edges_rewired_total,
      edges_rewired_fix_trial = result$edges_rewired_fix_trial,
      edges_rewired_fix_treatment = result$edges_rewired_fix_treatment,
      edges_rewired_fix_connected = result$edges_rewired_fix_connected
    )
  }
  names(networks) <- paste0("copula_MNS_tracked_", seq_len(n_qualified))
  list(
    networks = networks,
    summary = do.call(rbind, summary_rows),
    n_qualified = n_qualified,
    MNS_attempts = MNS_attempt,
    rejected_MNS_attempts = MNS_attempt - n_qualified,
    structurally_infeasible = structurally_infeasible,
    generation_attempts_exhausted = generation_attempts_exhausted,
    simulation_schema_version = 1L,
    copula_family = copula_model$family,
    copula_coefficients = copula_model$coefficients
  )
}
