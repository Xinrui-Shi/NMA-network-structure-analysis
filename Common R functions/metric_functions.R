# Metric calculations used by the empirical, null-model, and comparison reports.
# Deliberately omitted: neighbourhoods, centrality, degree distributions, and
# all other metrics not plotted in the report.

compute_bipartite_degree_assortativity <- function(B) {
  B <- 1L * (as.matrix(B) != 0)
  incidence_index <- which(B != 0L, arr.ind = TRUE)
  if (nrow(incidence_index) == 0L) return(NA_real_)

  trial_endpoint_degrees <- rowSums(B)[incidence_index[, "row"]]
  treatment_endpoint_degrees <- colSums(B)[incidence_index[, "col"]]

  # Each trial-treatment incidence contributes exactly one ordered pair:
  # (degree of the trial, degree of the treatment). Pearson correlation is
  # undefined when either endpoint-degree vector is constant.
  if (
    length(unique(trial_endpoint_degrees)) <= 1L ||
      length(unique(treatment_endpoint_degrees)) <= 1L
  ) {
    return(NA_real_)
  }

  value <- stats::cor(trial_endpoint_degrees, treatment_endpoint_degrees)
  if (is.finite(value)) value else NA_real_
}

compute_network_metrics <- function(B, network_id = NA_character_) {
  if (is.null(B) || length(B) == 0L || nrow(B) == 0L || ncol(B) == 0L) {
    stop("B must be a non-empty trial-by-treatment matrix.")
  }
  B <- 1L * (as.matrix(B) != 0)
  B <- B[rowSums(B) >= 2L, , drop = FALSE]
  B <- B[, colSums(B) > 0L, drop = FALSE]
  if (nrow(B) == 0L || ncol(B) == 0L) stop("No valid trial-treatment data remain.")

  g_bipartite <- igraph::graph_from_biadjacency_matrix(B, directed = FALSE)
  bi_components <- igraph::components(g_bipartite)
  g_bipartite_giant <- igraph::induced_subgraph(
    g_bipartite,
    which(bi_components$membership == which.max(bi_components$csize))
  )
  trial_degrees <- rowSums(B)
  treatment_degrees <- colSums(B)
  mean_trial_degree <- mean(trial_degrees)
  mean_treatment_degree <- mean(treatment_degrees)

  weighted_projection <- crossprod(B)
  diag(weighted_projection) <- 0
  unweighted_projection <- 1L * (weighted_projection > 0L)
  g_unipartite <- igraph::graph_from_adjacency_matrix(
    unweighted_projection, mode = "undirected", diag = FALSE
  )
  g_weighted <- igraph::graph_from_adjacency_matrix(
    weighted_projection, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  uni_components <- igraph::components(g_unipartite)
  g_unipartite_giant <- igraph::induced_subgraph(
    g_unipartite,
    which(uni_components$membership == which.max(uni_components$csize))
  )
  weighted_degrees <- igraph::strength(
    g_weighted, mode = "all", weights = igraph::E(g_weighted)$weight
  )
  unweighted_degrees <- igraph::degree(g_unipartite, mode = "all")
  mean_weighted_degree <- mean(weighted_degrees)
  mean_unweighted_degree <- mean(unweighted_degrees)

  incidence_index <- which(B != 0L, arr.ind = TRUE)
  tnet_edges <- cbind(
    top_node = incidence_index[, "row"],
    bottom_node = incidence_index[, "col"]
  )

  safe_assortativity <- function(graph) {
    degrees <- igraph::degree(graph)
    if (igraph::ecount(graph) == 0L || length(unique(degrees)) <= 1L) {
      return(NA_real_)
    }
    value <- igraph::assortativity_degree(graph, directed = FALSE)
    if (is.finite(value)) value else NA_real_
  }

  network_type <- if (!igraph::is_connected(g_unipartite)) {
    "disconnected"
  } else if (igraph::is_tree(g_unipartite)) {
    if (max(igraph::degree(g_unipartite)) == igraph::vcount(g_unipartite) - 1L) {
      "star"
    } else {
      "tree"
    }
  } else {
    "looped"
  }

  data.frame(
    network_id = as.character(network_id),
    N_treatments = ncol(B),
    N_total_trials = nrow(B),
    number_of_arms = igraph::ecount(g_bipartite),
    number_of_distinct_treatment_comparisons = igraph::ecount(g_unipartite),
    mean_bipartite_trial_degree = mean_trial_degree,
    mean_bipartite_treatment_degree = mean_treatment_degree,
    mean_unipartite_weighted_degree = mean_weighted_degree,
    mean_unipartite_unweighted_degree = mean_unweighted_degree,
    density = igraph::edge_density(g_unipartite),
    bi_density = sum(B) / (nrow(B) * ncol(B)),
    avg_path_length = if (igraph::vcount(g_unipartite_giant) > 1L) {
      igraph::mean_distance(g_unipartite_giant, directed = FALSE)
    } else NA_real_,
    diameter = igraph::diameter(
      g_unipartite, directed = FALSE, weights = NA
    ),
    bi_apl = if (igraph::vcount(g_bipartite_giant) > 1L) {
      igraph::mean_distance(g_bipartite_giant, directed = FALSE)
    } else NA_real_,
    unipartite_c = igraph::transitivity(g_unipartite, type = "global"),
    bipartite_c = tnet::reinforcement_tm(tnet_edges),
    assortativity = safe_assortativity(g_unipartite),
    bi_assortativity = compute_bipartite_degree_assortativity(B),
    deg_irreg = if (is.finite(mean_weighted_degree) && mean_weighted_degree > 0) {
      mean((weighted_degrees - mean_weighted_degree)^2) / mean_weighted_degree^2
    } else NA_real_,
    bi_trial_deg_irreg = if (
      is.finite(mean_trial_degree) && mean_trial_degree > 0
    ) {
      mean((trial_degrees - mean_trial_degree)^2) / mean_trial_degree^2
    } else NA_real_,
    bi_treatment_deg_irreg = if (
      is.finite(mean_treatment_degree) && mean_treatment_degree > 0
    ) {
      mean((treatment_degrees - mean_treatment_degree)^2) /
      mean_treatment_degree^2
    } else NA_real_,
    type = network_type,
    check.names = FALSE
  )
}

compute_metrics_for_networks <- function(networks) {
  if (is.null(names(networks))) names(networks) <- as.character(seq_along(networks))
  dplyr::bind_rows(lapply(seq_along(networks), function(i) {
    compute_network_metrics(networks[[i]], names(networks)[i])
  }))
}

compute_edge_multiplicity <- function(B, network_id = NA_character_) {
  if (is.null(B) || length(B) == 0L || nrow(B) == 0L || ncol(B) == 0L) {
    stop("B must be a non-empty trial-by-treatment matrix.")
  }
  B <- 1L * (as.matrix(B) != 0)
  weighted_projection <- crossprod(B)
  multiplicities <- weighted_projection[
    upper.tri(weighted_projection) & weighted_projection > 0L
  ]
  if (length(multiplicities) == 0L) {
    return(data.frame(
      network_id = character(), multi = integer(), count = integer()
    ))
  }
  multiplicity_counts <- as.data.frame(
    table(multiplicities), stringsAsFactors = FALSE
  )
  data.frame(
    network_id = as.character(network_id),
    multi = as.integer(as.character(multiplicity_counts$multiplicities)),
    count = as.integer(multiplicity_counts$Freq),
    stringsAsFactors = FALSE
  )
}

compute_edge_multiplicity_for_networks <- function(networks) {
  if (is.null(names(networks))) names(networks) <- as.character(seq_along(networks))
  dplyr::bind_rows(lapply(seq_along(networks), function(i) {
    compute_edge_multiplicity(networks[[i]], names(networks)[i])
  }))
}
