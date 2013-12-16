\name{gagaReport}
\alias{gagaReport}
\title{View and create output of GAGA obects}
\usage{
gagaReport()
}
\arguments{
  \item{o_file_prefix}{Output file prefix - SHOULD BE MOVED
  TO REPORT OUTPUT FUNCTION}

  \item{contamination}{Is the input contaminated?  If set
  to 1, an extra clone is created in which to place
  contaminants}
}
\value{
Description of the returned object
}
\description{
Description here
}
\examples{
gaga("BYB1-G07_pruned.txt", "BYB-G07_output_pruned_contaminated", "BYB1-G07_anno_pruned.txt", 6, 6, iterations=10, contamination=1)
}
\author{
Alex Murison and Christopher Wardell
\email{Alexander.Murison@icr.ac.uk}
}

