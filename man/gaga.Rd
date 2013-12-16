\name{gaga}
\alias{gaga}
\title{Genetic Algorithm for Generating Ancestry (GAGA)}
\usage{
gaga(observations, annotations, number_of_clones, pop_size = 100,
  mutation_rate = 0.8, iterations = 1000,
  stoppingCriteria = round(iterations/10), parthenogenesis = 2, nroot = 0,
  contamination = 0, seed)
}
\arguments{
  \item{observations}{Observation data frame.  Rows
  represent the proportion of cells containing a mutation,
  columns represent discrete samples separated by time or
  space}

  \item{annotations}{Annotation data frame.  Each row
  corresponds to a mutation ??????}

  \item{number_of_clones}{The integer number of clones to
  be considered}

  \item{pop_size}{The number of individuals in each
  generation}

  \item{mutation_rate}{The likelihood of each individual
  undergoing mutation per generation}

  \item{iterations}{The number of generations that will
  occur}

  \item{stoppingCriteria}{The number of consecutive
  generations without improvement that will stop the
  algorithm. Default value is 10\% of iterations.}

  \item{parthenogenesis}{The number of best-fitness
  individuals allowed to survive each generation}

  \item{nroot}{Number of roots the phylogeny is expected to
  have}

  \item{contamination}{Is the input contaminated?  If set
  to 1, an extra clone is created in which to place
  inferred contaminants}

  \item{seed}{An integer vector to set the random number
  generator.  Allows runs of the algorithm to be repeated}
}
\value{
Returns an object of class ga \code{\link{ga-class}}
}
\description{
Use a genetic algorithm to find the relationships between
the values in an input file.  Generates number of clones,
the phylogenetic relationship between them and the
proportion of each clone that each sample is composed of.
}
\examples{
# Load the included synthetic example data
data("gaga_synthetic_data","gaga_synthetic_data_annotation")
x=gaga(gaga_synthetic_data, gaga_synthetic_data_annotation, number_of_clones=6, iterations=100)
}
\author{
Alex Murison \email{Alexander.Murison@icr.ac.uk} and
Christopher Wardell \email{Christopher.Wardell@icr.ac.uk}
}
\seealso{
\code{\link{ga-class}}, \code{\link{ga}}
}

