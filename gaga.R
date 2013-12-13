###############################################################################
## Genetic Algorithm for Generating Ancestry  (Gaga)                         ##
## Alex Murison & Christopher Wardell                                        ##
## Institute of Cancer Research, Sutton, London, UK                          ##
## testing Version 0.2 - 13/12/2013                                          ##
## e-mail: Alexander.Murison@icr.ac.uk                                       ##
###############################################################################


gaga<-function(ob_file, o_file_prefix, annotations_file, min_number_of_clones, max_number_of_clones, 
               pop_size=100, mutation_rate=0.8, iterations=1000, parthenogenesis=2, loopy=1, nroot=0, contamination=0) {
  
  new_phylo<-function(nroot) {
    # Generates a phylogeny with one and only one root node
    
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
    
    return(inp)
  }
  
  phylo_cross<-function(object, parents) {
    # Crossover as gabin_SpCrossover but do not select a crossover point within a phylogeny
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
  
  generate_phylo_matrix<-function(input) {
    
    
    ## Creates a matrix where rows indicate the clones in which a mutation must occur.
    # Identify the root nodes, set these to have 1 in the diagonal corresponding to their column
    # For each root node, find their children, copy the parents column and add 1 to their diagonal
    # find the children of those children and repeat until matrix constructed,
    
    
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
  
  generate_phylogeny_aware_population<-function(object) {
    #print("popgen")
    population <- matrix(NA, nrow = object@popSize, ncol = object@nBits)
    for (j in 1:object@popSize) {
      population[j, ] <- generate_phylogeny_aware_individual()
    }
    return(population)
    
    
  }
  
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
  
  write_the_output<-function() {
    
    ## Load the required packages
    library(heatmap.plus)
    library(png)
    
    #print("nyah")
    # Set ICR Colour Scheme
    buylrd=c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026") # blue yellow red 
    primary=c("#C9DD03","#FFD602","#F9A100","#EE7EA6","#A71930","#616365") # primary ICR colours
    secondary=c("#726E20","#6E273D","#F00034","#ADAFAF","#003D4C") # secondary ICR colours
    colours=c(primary,secondary[c(1,2,3,5)]) # exclude ICR light grey
    
    #Read in the ICR Branding
    ima<-readPNG("ICR.png")  
    
    ## Load annotation data
    snvs=read.table(file=annotations_file,header=TRUE,stringsAsFactors=FALSE)
    #print("nyah")
    print("Outputting the highest scoring solutions")
    
    top_number<-min(5,nrow(goo@solution))
    for (faint_o_ateb in 1:top_number) {
      #best_result<-matrix(goo@solution[faint_o_ateb,], ncol=number_of_clones, byrow=TRUE)
      ateb<-goo@solution[faint_o_ateb,]
      phylogeny_matrix<-generate_phylo_matrix(ateb[1:number_of_clones])
      proportion_matrix<-matrix(ateb[(number_of_clones+1):((number_of_cases*pseudo_number_of_clones)+number_of_clones)], 
                                ncol=pseudo_number_of_clones, nrow=number_of_cases, byrow=TRUE)
      phylogeny_matrix<-rbind(phylogeny_matrix, ateb[1:number_of_clones])
      #if(contamination==1) {cbind(phylogeny_matrix, c(-1)}
      
      for (rows in 1:nrow(proportion_matrix)) {
        scale<-sum(proportion_matrix[rows,])
        proportion_matrix[rows,]<-proportion_matrix[rows,]/scale
      }
      presence_matrix<-matrix(nrow=number_of_mutations, ncol=number_of_clones)
      for (rows in 1:number_of_mutations) {
        presence_matrix[rows,]<-phylogeny_matrix[ateb[((number_of_cases*pseudo_number_of_clones)+number_of_clones+rows)],]
      }
      write.table(ateb,paste(o_file_prefix, ".",number_of_clones,"clones.complete.solution.",faint_o_ateb,"rep", loop_value, ".txt", sep=""),sep="\t", row.names=FALSE)
      write.table(phylogeny_matrix,paste(o_file_prefix, ".",number_of_clones,"clones.phylogeny.solution.",faint_o_ateb,"rep", loop_value, ".txt", sep=""),sep="\t", row.names=FALSE)
      write.table(proportion_matrix,paste(o_file_prefix, ".",number_of_clones,"clones.proportions.solution.",faint_o_ateb,"rep", loop_value, ".txt", sep=""),sep="\t", row.names=FALSE)
      write.table(presence_matrix,paste(o_file_prefix, ".",number_of_clones,"clones.mutation.matrix.",faint_o_ateb,"rep", loop_value, ".txt", sep=""),sep="\t", row.names=FALSE)
      
      clones=presence_matrix
      colnames(clones)=1:ncol(clones) 
      
      
      ## Iterate through presence_matrix data and replace values with hex colours
      for(i in 1:ncol(clones)){
        ## Use the modulo remainder NOT i, so we can recycle the same colour object
        ## regardless of the number of presence_matrix
        recycler=ifelse(i%%length(colours)==0,length(colours),i%%length(colours))
        clones[clones[,i]=="1",i]=colours[recycler]
        clones[clones[,i]=="0",i]=secondary[4] # background is ICR light grey
      }
      
      
    }
    
    
    pdf(paste(o_file_prefix,".solution.",faint_o_ateb, ".clone.",number_of_clones,"rep", loop_value, ".heatmap_and_proportions.pdf", sep=""), paper="a4")
    ## Produce heatmap
    #x=as.matrix(snvs[,colnames(observation_matrix)])
    heatmap.plus(as.matrix(observation_matrix),Colv=NA,col=buylrd,scale="none",Rowv=NULL,labRow=snvs$names,RowSideColors=clones,cexRow=0.9)
    
    props = proportion_matrix
    rownames(props)=colnames(observation_matrix)
    barplot(t(props),col=colours,border=NA,xlab="Sample name",ylab="Proportion of cells")
    frame()
    lim<-par()
    rasterImage(ima, lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
    
    dev.off()
    
    
    
    #  print("Outputting Convergence plot")
    png(paste(o_file_prefix, ".",number_of_clones,"rep", loop_value, ".convergence.png", sep=""))
    plot(goo)
    dev.off()
    
    write.table(goo@best,paste(o_file_prefix, ".",number_of_clones,"clones.", loop_value, ".reps.highest.score.solution.",faint_o_ateb,".txt", sep=""),sep="\t", row.names=FALSE)
    
    
    print(summary(goo))
    
  }
  
  
  library(GA)
  
  # Get Observation File and associated values
  observation_matrix<-read.table(ob_file, sep="\t", stringsAsFactors=FALSE, header=T)
  number_of_mutations<-as.numeric(nrow(observation_matrix))
  number_of_cases<-as.numeric(ncol(observation_matrix))
  
  
  
  
  for (loop_value in 1:loopy) {
    
    for (number_of_clones in min_number_of_clones:max_number_of_clones) {
      
      if (contamination == 1) {
        pseudo_number_of_clones<-number_of_clones+1
        numBits<-(pseudo_number_of_clones*number_of_cases)+(number_of_mutations+number_of_clones)    
      } else {
        numBits<-(number_of_clones*number_of_cases)+(number_of_mutations+number_of_clones)
        pseudo_number_of_clones<-number_of_clones
      }
      
      ### Compile the genetic algorithm
      library(compiler)
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
      
      write_the_output()
      
    }
    
  }
  
}

## Example command:
gaga("BYB1-G07_pruned.txt", "BYB-G07_output_pruned_contaminated", "BYB1-G07_anno_pruned.txt", 6, 6, iterations=10, contamination=1, loopy=1)

################################
## Phylogeny plots begin here ##
################################
library(graph)
number_of_clones=16

## Create new phylogeny
p=new_phylo(0)

## Get phylogeny matrix from AM's function
pm=generate_phylo_matrix(p)

## These matrices are not compliant with Rgraphviz, so we edit them
## Matrices must have 0 diag, sum(column)<=1, keep leftmost value
diag(pm)=0
for(i in 1:nrow(pm)){
  for(j in ncol(pm):1){
    if(sum(pm[,j])>1){
      pm[i,j]=0
    }
  }
}
## Set lineage names and cast to graphNEL object
rownames(pm)=LETTERS[1:nrow(pm)]
colnames(pm)=rownames(pm)
pg=as(pm,"graphNEL")

## Set up plotting parameters and set colours
primary=c("#C9DD03","#FFD602","#F9A100","#EE7EA6","#A71930","#616365") # primary ICR colours
secondary=c("#726E20","#6E273D","#F00034","#ADAFAF","#003D4C") # secondary ICR colours
colours=c(primary,secondary[c(1,2,3,5)]) # exclude ICR light grey
nAttrs = list() ## per node attributes
nAttrs$fillcolor = rep(NA,nrow(pm))
## Iterate through nAttrs$fillcolor and replace values with hex colours
for(i in 1:nrow(pm)){
  ## Use the modulo remainder NOT i, so we can recycle the same colour object
  recycler=ifelse(i%%length(colours)==0,length(colours),i%%length(colours))
  nAttrs$fillcolor[i]=colours[recycler]
}
names(nAttrs$fillcolor) = LETTERS[1:nrow(pm)]
nAttrs$color = rep(NA,nrow(pm))
names(nAttrs$color) = LETTERS[1:nrow(pm)] ## Assign NA color to shape to ensure clean fill
plot(pg, nodeAttrs = nAttrs, attrs=list(edge = list(color=secondary[4],lwd=4))) # "attrs" sets global attributes

##############################
## Phylogeny plots end here ##
##############################
