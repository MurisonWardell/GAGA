\name{gaga}
\alias{gaga}
\title{Genetic Algorithm for Generating Ancestry (GAGA)}
\usage{
gaga(observations, number_of_clones, pop_size = 100, mutation_rate = 0.8,
  iterations = 3000, stoppingCriteria = round(iterations/5),
  parthenogenesis = 2, nroot = 1, contamination = 0)
}
\arguments{
  \item{observations}{Observation data frame where each row
  represents an SNV and each column (except the first)
  represents a discrete sample separated by time or space.
  Note that the first row is a header containing the name
  of each sample and the first column must contain a unique
  name for each SNV.  Every value must be a proportion
  between 0 and 1}

  \item{number_of_clones}{An integer number of clones to be
  considered}

  \item{pop_size}{The number of individuals in each
  generation}

  \item{mutation_rate}{The likelihood of each individual
  undergoing mutation per generation}

  \item{iterations}{The maximum number of generations to
  run}

  \item{stoppingCriteria}{The number of consecutive
  generations without improvement that will stop the
  algorithm. Default value is 20\% of total number of
  iterations.}

  \item{parthenogenesis}{The number of best-fitness
  individuals allowed to survive each generation}

  \item{nroot}{Number of roots in phylogenies.  Either
  specify an integer equal to or less than the number of
  clones or 0.  When nroot=0, a random integer between 1
  and the number of clones is generated for each phylogeny}

  \item{contamination}{Is the input contaminated?  If set
  to 1, an extra (empty) clone is created which allows the
  sum of the clones to be less than 1.  Analagous to normal
  cell contamination in a cancer sample}
}
\value{
Returns an object of class ga \code{\link{ga-class}}.  Note
that the number of clones and number of cases are stored in
the unused min and max slots of the output object.
}
\description{
Use a genetic algorithm to find the relationships between
the values in an input file - the package was written to
deal with single nucleotide variants (SNVs) in mixtures of
cancer cells, but it will work with any mixture.  It will
calculate appropriate phylogenetic relationships between
clones them and the proportion of each clone that each
sample is composed of.  For detailed usage, please read the
accompanying vignette.
}
\details{
Do we need a details section?
}
\examples{
## Load the included synthetic example data
data("gaga_synthetic_data")

## Run gaga using 6 clones and default parameters
solution=gaga(gaga_synthetic_data, number_of_clones=6)

## Output the solution in the current working directory
gagaReport(gaga_synthetic_data,solution)
}
\author{
Alex Murison \email{Alexander.Murison@icr.ac.uk} and
Christopher Wardell \email{Christopher.Wardell@icr.ac.uk}
}
\seealso{
\code{\link{ga-class}}, \code{\link{ga}},
\code{\link{gagaReport}}
}

