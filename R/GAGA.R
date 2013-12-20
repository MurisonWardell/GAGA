#' Genetic Algorithm for Generating Ancestry (GAGA)
#' 
#' Use a genetic algorithm to find the relationships between the values in an input file - the package was written to deal with single
#' nucleotide variants (SNVs) in mixtures of cancer cells, but it will work with any mixture.  It will calculate appropriate phylogenetic
#' relationships between clones them and the proportion of each clone that each sample is composed of.  For detailed usage, please
#' read the accompanying vignette.
#' 
#' @param observations              Observation data frame where each row represents an SNV and each column (except the first) represents
#' a discrete sample separated by time or space.  Note that the first row is a header containing the name of each sample and the first
#' column must contain a unique name for each SNV.  Every value must be a proportion between 0 and 1
#' @param number_of_clones          An integer number of clones to be considered
#' @param pop_size                  The number of individuals in each generation
#' @param mutation_rate             The likelihood of each individual undergoing mutation per generation
#' @param iterations                The maximum number of generations to run
#' @param stoppingCriteria          The number of consecutive generations without improvement that will stop the algorithm.
#' Default value is 20\% of total number of iterations.
#' @param parthenogenesis           The number of best-fitness individuals allowed to survive each generation
#' @param nroot                     Number of roots in phylogenies.  Either specify an integer equal to or less than the number of clones
#' or 0.  When nroot=0, a random integer between 1 and the number of clones is generated for each phylogeny
#' @param contamination             Is the input contaminated?  If set to 1, an extra (empty) clone is created which allows the sum of
#' the clones to be less than 1.  Analagous to normal cell contamination in a cancer sample
#' @details Do we need a details section?
#' @return Returns an object of class ga \code{\link{ga-class}}.  Note that the number of clones and number of cases are 
#' stored in the unused min and max slots of the output object.
#' @export
# @references GAGA paper here!
#' @author Alex Murison \email{Alexander.Murison@@icr.ac.uk} and Christopher Wardell \email{Christopher.Wardell@@icr.ac.uk}
#' @seealso \code{\link{ga-class}}, \code{\link{ga}}, \code{\link{gagaReport}}
#' @examples
#' ## Load the included synthetic example data
#' data("gaga_synthetic_data")
#' 
#' ## Run gaga using 6 clones and default parameters
#' solution=gaga(gaga_synthetic_data, number_of_clones=6)
#' 
#' ## Output the solution in the current working directory
#' gagaReport(gaga_synthetic_data,solution)

gaga<-function(observations, number_of_clones, pop_size=100, mutation_rate=0.8, iterations=3000,
               stoppingCriteria=round(iterations/5), parthenogenesis=2,nroot=1, contamination=0) {

  
  ##############################
  ## Start internal functions ##
  ##############################
  
  ## Crossover as gabin_SpCrossover but do not select a crossover point within a phylogeny
  phylo_cross<-function(object, parents) {
    fitness <- object@fitness[parents]
    parents <- object@population[parents, , drop = FALSE]
    n <- ncol(parents)
    children <- matrix(NA, nrow = 2, ncol = n)
    fitnessChildren <- rep(NA, 2)
    crossOverPoint <- sample((number_of_clones+1):n, size = 1)
    if (crossOverPoint == 0) {
      children[1:2, ] <- parents[2:1, ]
      fitnessChildren[1:2] <- fitness[2:1]
    }
    else if (crossOverPoint == n) {
      children <- parents
      fitnessChildren <- fitness
    }
    else {
      children[1, ] <- c(parents[1, 1:crossOverPoint], parents[2, 
                                                               (crossOverPoint + 1):n])
      children[2, ] <- c(parents[2, 1:crossOverPoint], parents[1, 
                                                               (crossOverPoint + 1):n])
      fitnessChildren <- NA
    }
    out <- list(children = children, fitness = fitnessChildren)
    return(out)
  }
  ## /end of phylo_cross

  
  ## Calculate fitness of phylogeny
  fit_phylogeny<-function(input) {
    
    ### Work out the phylogeny matrix
    dep_mat<-generate_phylo_matrix(input[1:number_of_clones],number_of_clones)
    
    ### Calculate the proportions of each clone in each case
    proportion_matrix<-matrix(input[(number_of_clones+1):((number_of_cases*pseudo_number_of_clones)+number_of_clones)], 
                              ncol=pseudo_number_of_clones, nrow=number_of_cases, byrow=TRUE)
    for (rows in 1:nrow(proportion_matrix)) {
      scale<-sum(proportion_matrix[rows,])
      proportion_matrix[rows,]<-proportion_matrix[rows,]/scale
    }
    
    tot_score<-0
    
    # For each mutation, obtain the clones it appears in from the phylogeny matrix
    # Calculate the expected frequency in each case by summing clone frequencies for each mutation in each case
    prediction_matrix<-matrix(ncol=number_of_cases, nrow=number_of_mutations)
    rownum<-1
    for (mutation in input[((number_of_clones)+(pseudo_number_of_clones*number_of_cases)+1):length(input)]) {
      affected_clones<-dep_mat[mutation,]
      for (case in 1:number_of_cases) {
        score<-0
        for (clone in 1:number_of_clones) {         
          score<-score+(affected_clones[clone]*proportion_matrix[case,clone]) ## 16% of total load
        }
        prediction_matrix[rownum,case]<-score
        
      }
      rownum<-rownum+1
    }
    
    # calculate the deviance for each case/mutation from the observations
    # return -1 times the total deviance for all mutations
    vals<-prediction_matrix-observation_matrix ## 11% of total load
    for (plink in 1:nrow(vals)) {
      tot_score<-tot_score+sum(abs(vals[plink,])) ## 62% of total load
    }
    return(-1*tot_score)
    
  }
  # /end of fit_phylogeny
  
  ## Generates a population of individuals
  generate_phylogeny_aware_population<-function(object) {
    #print("popgen")
    population <- matrix(NA, nrow = object@popSize, ncol = object@nBits)
    for (j in 1:object@popSize) {
      population[j, ] <- generate_phylogeny_aware_individual()
    }
    return(population)
  }
  # /end of generate_phylogeny_aware_population
  
  ## Generates an individual
  generate_phylogeny_aware_individual<-function() {
    #  print("indivgen")
    
    # work out how many proportions
    clo_cas<-(pseudo_number_of_clones*number_of_cases)
    #work out how many mutations
    how_big<-number_of_clones+clo_cas+number_of_mutations
    # create our new individual
    new_individual<-vector(length=how_big)
    
    #Assign parents from 0 to number of clones
    new_individual[1:number_of_clones]<-new_phylo(nroot,number_of_clones)
    # assign proportions between 0 and 5
    new_individual[(number_of_clones+1):(number_of_clones+clo_cas)]<-round(runif(clo_cas,0,5))
    # assign mutation first appears as 1:number of clones
    new_individual[(number_of_clones+clo_cas+1):how_big]<-round(runif(number_of_mutations,0.5,(number_of_clones+0.49)))
    # return new individual
    return(new_individual)
  }
  # /end of generate_phylogeny_aware_individual
  
  ## Mutation function
  phylo_mutate<-function(object, parent) {
    
    #  print("mutator")
    # Turn the string into a vector
    victim_ready <- parent <- as.vector(object@population[parent, ])
    #<-as.numeric(unlist(strsplit(unsuspecting_victim, "")))
    
    # work out where we are going to mutate - 
    # 0 a proportion of clonal population
    # 1 mutation presence/absence in a clone
    # Randomly select with equal probability of both
    what<-round(runif(1,0,2))
    if (what == 0) {
      # Mutate Proportions
      # Randomly select a proportion
      where<-round(runif(1,(number_of_clones+0.5),((number_of_cases*pseudo_number_of_clones)+number_of_clones+0.49)))
      new<-victim_ready[where]
      # if it's already 0 increase it by 1
      if (new==0) {
        new<-new+1
        victim_ready[where]<-new
      }
      else { 
        # Otherwise randomly increment or decrement by 1
        do_what<-round(runif(1,0,1))
        if (do_what==1) {new<-new+1}
        if (do_what==0) {new<-new-1}
        victim_ready[where]<-new
      }
      
    }
    else if (what == 1 ) {
      # Mutate parents
      
      victim_ready[1:number_of_clones]<-new_phylo(nroot,number_of_clones)
    }
    else if (what == 2 ) {
      # Mutate a first_appears call
      where<-round(runif(1,((number_of_cases*pseudo_number_of_clones)+number_of_clones+0.5),length(victim_ready)))
      victim_ready[where]<-round(runif(1,0.5,(number_of_clones+0.49)))  
    }
    
    # Return it
    return(victim_ready)
  }
  # /end of phylo_mutate
  
  
  ############################
  ## End internal functions ##
  ############################
  
  
  #########################
  ## Start of main logic ##
  #########################
  
  library(GA)
  library(compiler)
  
  ## Get Observation File and associated values
  ## We ignore the first column, as it is the object names
  observation_matrix<-observations[,2:ncol(observations)] 
  number_of_mutations<-as.numeric(nrow(observation_matrix))
  number_of_cases<-as.numeric(ncol(observation_matrix))

  if (contamination == 1) {
    pseudo_number_of_clones<-number_of_clones+1
    numBits<-(pseudo_number_of_clones*number_of_cases)+(number_of_mutations+number_of_clones)    
  } else {
    numBits<-(number_of_clones*number_of_cases)+(number_of_mutations+number_of_clones)
    pseudo_number_of_clones<-number_of_clones
  }
  
  ### Compile the genetic algorithm - this increases speed significantly
  setCompilerOptions(suppressAll=TRUE)
  ga=cmpfun(ga)
    
  ### RUN THE GENETIC ALGORITHM
  goo<-ga(type="binary", 
          fitness = fit_phylogeny, 
          population = generate_phylogeny_aware_population,
          selection = gabin_tourSelection,
          crossover = phylo_cross,
          mutation = phylo_mutate,
          popSize = pop_size,
          nBits = numBits,
          pmutation=mutation_rate,
          maxiter=iterations,
          elitism = parthenogenesis,
          run=stoppingCriteria
  )
  ## Add annotation to the object
  goo@names=as.character(observations[,1])
    
  ## Add number of clones to min slot.  Note that this is an incorrect use of the slot
  goo@min = number_of_clones
  ## Add number of cases (e.g. timepoints) to max slot.  Note that this is an incorrect use of the slot
  goo@max= number_of_cases

  ## Define a new class and create the object
  #setClass("gaga", contains="ga",slots = c(id = "character"))
  
  
  #object <- new("ga", call = call, type = type, min = min, 
  #              max = max, nBits = nBits, names = if (is.null(names)) 
  #                character()
  #              else names, popSize = popSize, iter = 0, run = 1, maxiter = maxiter, 
  #              suggestions = suggestions, population = matrix(), elitism = elitism, 
  #              pcrossover = pcrossover, pmutation = pmutation, fitness = Fitness, 
  #              best = bestEval, mean = meanEval, bestSol = bestSol)
  
  
  return(goo)

  #######################
  ## End of main logic ##
  #######################
  
}
