###############################################################################
## Genetic Algorithm for Generating Ancestry  (GAGA)                         ##
## Alex Murison & Christopher Wardell                                        ##
## Institute of Cancer Research, Sutton, London, UK                          ##
## testing Version 0.2 - 13/12/2013                                          ##
## e-mail: Alexander.Murison@icr.ac.uk                                       ##
###############################################################################


#' Genetic Algorithm for Generating Ancestry (GAGA)
#' 
#' Use a genetic algorithm to find the relationships between the values in an input file.  Generates number of clones, the
#' phylogenetic relationship between them and the proportion of each clone that each sample is composed of.
#' 
#' @param observations              Observation data frame.  Rows represent the proportion of cells containing
#' a mutation, columns represent discrete samples separated by time or space
#' @param annotations               Annotation data frame.  Each row corresponds to a mutation ??????
#' @param min_number_of_clones      The minimum number of clones to be considered
#' @param max_number_of_clones      The maximum number of clones to be considered
#' @param pop_size                  The number of individuals in each generation
#' @param mutation_rate             The likelihood of each individual undergoing mutation per generation
#' @param iterations                The number of generations that will occur
#' @param parthenogenesis           The number of best-fitness individuals allowed to survive each generation
#' @param nroot                     Number of roots the phylogeny is expected to have
#' @param contamination             Is the input contaminated?  If set to 1, an extra clone is created in which to place contaminants
#' @return Description of the returned object
#' @export
# @seealso \code{\link{fermat.test}}  #### other functions; e.g. gaga report output
# @references GAGA paper here!
#' @author Alex Murison and Christopher Wardell \email{Alexander.Murison@@icr.ac.uk}
#' @examples
#' gaga("BYB1-G07_pruned.txt", "BYB-G07_output_pruned_contaminated", "BYB1-G07_anno_pruned.txt", 6, 6, iterations=10, contamination=1)

gaga<-function(observations, annotations, min_number_of_clones, max_number_of_clones, 
               pop_size=100, mutation_rate=0.8, iterations=1000, parthenogenesis=2, nroot=0, contamination=0) {

  ##############################
  ## Start internal functions ##
  ##############################
  
  ## Generates a phylogeny with one and only one root node
  new_phylo<-function(nroot) {
    inp<-vector(length=number_of_clones)
    
    if (nroot==1) { 
      i<-1
    } else { 
      i<-round(runif(1,0.5, number_of_clones+0.49)) 
    }
    
    inp[1:i]<-0
    
    last_gen_start<-1
    last_gen_end<-i
    spaces_left<-number_of_clones-i
    
    while (spaces_left > 0) {
      
      length_next_gen<-round(runif(1,0.5, spaces_left+0.49))
      current_sample<-sample(rep(last_gen_start:last_gen_end, number_of_clones))
      #<-sample(last_gen_start:last_gen_end, size=number_of_clones, replace=TRUE)
      inp[(last_gen_end+1):(last_gen_end+length_next_gen)]<-current_sample[1:length_next_gen]
      last_gen_start<-last_gen_end+1
      last_gen_end<-last_gen_end+length_next_gen
      spaces_left<-spaces_left-length_next_gen
      
    }
    ## /end of new_phylo
    
    return(inp)
  }
  
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
  
  
  ## Creates a matrix where rows indicate the clones in which a mutation must occur.
  ## Identify the root nodes, set these to have 1 in the diagonal corresponding to their column
  ## For each root node, find their children, copy the parents column and add 1 to their diagonal
  ## find the children of those children and repeat until matrix constructed
  generate_phylo_matrix<-function(input) {
    pmat<-matrix(rep(0, (number_of_clones*number_of_clones)), ncol=number_of_clones, byrow=TRUE)
    next_generation<-which(input==0)
    new_next_generation<-c()  
    for (node in next_generation) {
      pmat[node,node]<-1
      new_next_generation<-c(node, new_next_generation) 
    }
    next_generation<-new_next_generation
    while (length(next_generation>0)) {
      new_next_generation<-c() 
      for (i in next_generation) {
        getme<-which(input==i)
        for (j in getme) {
          pmat[,j]<-pmat[,i]
          pmat[j,j]<-1
          new_next_generation<-c(j, new_next_generation) 
        }
      }
      next_generation<-new_next_generation
    }
    return(pmat)
  }
  # /end of generate_phylo_matrix
  
  ## Calculate fitness of phylogeny
  fit_phylogeny<-function(input) {
    
    ### Work out the phylogeny matrix
    dep_mat<-generate_phylo_matrix(input[1:number_of_clones])
    
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
    new_individual[1:number_of_clones]<-new_phylo(nroot)
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
      
      victim_ready[1:number_of_clones]<-new_phylo(nroot)
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
  library(graph)
  library(heatmap.plus)
  library(png)
  
  # Get Observation File and associated values
  #observation_matrix<-read.table(observations, sep="\t", stringsAsFactors=FALSE, header=T)
  observation_matrix<-observations
  number_of_mutations<-as.numeric(nrow(observation_matrix))
  number_of_cases<-as.numeric(ncol(observation_matrix))
  
    
  for (number_of_clones in min_number_of_clones:max_number_of_clones) {
    
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
            elitism = parthenogenesis
    )
    
    return(goo)
    
  }

  #######################
  ## End of main logic ##
  #######################
  
}
