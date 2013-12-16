\name{gaga}
\alias{gaga}
\title{Genetic Algorithm for Generating Ancestry (GAGA)}
\usage{
gaga(observations, annotations, min_number_of_clones, max_number_of_clones,
  pop_size = 100, mutation_rate = 0.8, iterations = 1000,
  parthenogenesis = 2, nroot = 0, contamination = 0)
}
\arguments{
  \item{observations}{Observation data frame.  Rows
  represent the proportion of cells containing a mutation,
  columns represent discrete samples separated by time or
  space}

  \item{annotations}{Annotation data frame.  Each row
  corresponds to a mutation ??????}

  \item{min_number_of_clones}{The minimum number of clones
  to be considered}

  \item{max_number_of_clones}{The maximum number of clones
  to be considered}

  \item{pop_size}{The number of individuals in each
  generation}

  \item{mutation_rate}{The likelihood of each individual
  undergoing mutation per generation}

  \item{iterations}{The number of generations that will
  occur}

  \item{parthenogenesis}{The number of best-fitness
  individuals allowed to survive each generation}

  \item{nroot}{Number of roots the phylogeny is expected to
  have}

  \item{contamination}{Is the input contaminated?  If set
  to 1, an extra clone is created in which to place
  contaminants}
}
\value{
Description of the returned object
}
\description{
Use a genetic algorithm to find the relationships between
the values in an input file.  Generates number of clones,
the phylogenetic relationship between them and the
proportion of each clone that each sample is composed of.
}
\examples{
gaga("BYB1-G07_pruned.txt", "BYB-G07_output_pruned_contaminated", "BYB1-G07_anno_pruned.txt", 6, 6, iterations=10, contamination=1)
}
\author{
Alex Murison and Christopher Wardell
\email{Alexander.Murison@icr.ac.uk}
}

