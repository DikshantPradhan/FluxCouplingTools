% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gurobi_tools.R
\name{GRB_generate_set_lists_cluster}
\alias{GRB_generate_set_lists_cluster}
\title{Function for PACT method}
\usage{
GRB_generate_set_lists_cluster(model_og, suppression_idxs = -1,
  reaction_indexes = c(), compare_known_init_sets = FALSE,
  optimize_suppr = FALSE, optimize_rxns = FALSE, cores = 1,
  avoid_idxs = c(), file_output = NULL)
}
\arguments{
\item{model_og}{model upon which to run PACT}

\item{suppression_idxs}{list of integers indicating indexes of reactions to iteratively block and calculate coupling for. DEFAULT is -1, which indicates that all reactions should be suppressed}

\item{reaction_indexes}{list of integers indicating indexes of reaction to check coupling for. DEFAULT is empty list, which indicates that all reactions in the model should be considered}

\item{compare_known_init_sets}{boolean indicating whether or not to calculate R0 sets in order to run further optimizations during flux coupling and PACT. DEFAULT: FALSE}

\item{optimize_suppr}{boolean indicating whether or not to optimize the reactions to suppress based on R0 sets. DEFAULT: FALSE}

\item{optimize_rxns}{boolean indicating whether or not to optimize the reactions couple. DEFAULT: FALSE}

\item{cores}{integer indicating the number of cores to run on, if intending to parallelize}

\item{avoid_idxs}{list of integers indicating indexes to specifically avoid suppressing during PACT method}

\item{file_output}{filename to output coupling vector to}
}
\value{
a coupling list: list of list of integers which each indicate the couplings induced by each suppression. The index position of each list indicates the index of the reaction which was suppressed and the integers in eacch list indicate the positions of TRUE values in the coupling array
}
\description{
Function for PACT method
}
\examples{

}
