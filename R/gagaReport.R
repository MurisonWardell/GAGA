#' View and create output of GAGA objects
#' 
#' Description here
#' 
#' @param gagaInput             Raw data put into gaga() 
#' @param gagaOutput            Object of class ga produced by gaga()
#' @param outType               Type of output desired
#' @param yRange                Y-axis range when plotting fitness of individuals
#' @param output_file_prefix    Optional prefix for output files
#' @return Description of the returned object
#' @export
# @seealso \code{\link{fermat.test}}  #### other functions; e.g. gaga report output
# @references GAGA paper here!
#' @author Alex Murison and Christopher Wardell \email{Alexander.Murison@@icr.ac.uk}
#' @examples
#' gagaReport(blah,UNKNOWN)


gagaReport<-function(gagaInput,gagaOutput,outType="complete",yRange=c(-250,0),output_file_prefix="") {
  
  #################
  ## Input check ##
  #################
  
  ## Check that the outType argument is correctly specified
  if(!(outType %in% c("complete","fitness","heatmap","phylogeny","proportion"))){
    stop("outType incorrect: must be one of the following: \"complete\",\"fitness\",\"heatmap\",\"phylogeny\",\"proportion\"")
  }
  
  #####################
  ## End input check ##
  #####################
  
  ## Load libraries
  library(graph)
  library(heatmap.plus)
  library(png)
  
  ##################################################
  ## Fetch, calculate and otherwise set variables ##
  ##################################################
  
  ## Set observation matrix for heatmap
  observation_matrix = gagaInput
  
  ## Get number of mutations from gagaOutput
  number_of_mutations=length(gagaOutput@names)
  
  ## Get number of clones from gagaOutput
  number_of_clones=gagaOutput@min
  
  ## Get number of cases (e.g. timepoints) from gagaOutput
  number_of_cases=gagaOutput@max
  
  ## We can now work out contamination status.  This is the reverse of calculating numBits in gaga()
  if(((gagaOutput@nBits-number_of_clones-number_of_mutations)/number_of_cases)==number_of_clones){
    contamination=0
  }else{
    contamination=1
  }
    
  ## Now we now the contamination status, we can calculate the pseudo_number_of_clones
  if (contamination == 1) {
    pseudo_number_of_clones<-number_of_clones+1
  } else {
    pseudo_number_of_clones<-number_of_clones
  }
      
  ## Set colour schemes
  heatmapColours=c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026") # blue yellow red 
  primary=c("#C9DD03","#FFD602","#F9A100","#EE7EA6","#A71930","#616365") # primary ICR colours
  secondary=c("#726E20","#6E273D","#F00034","#ADAFAF","#003D4C") # secondary ICR colours
  colours=c(primary,secondary[c(1,2,3,5)]) # exclude ICR light grey
  
  ######################
  ## End of variables ##
  ######################

  
  ########################################
  ## Begin looping over TOP 5 solutions ##
  ########################################
   
  ## Look at top solution(s) and parse into usable objects
  ## NOTE: if there are several equally good solutions, it will print (up to) the top five solutions
  top_number<-min(5,nrow(gagaOutput@solution))
  for (number_of_answers in 1:top_number) {
    
    #######################################
    ## Begin plotting / producing output ##
    #######################################
    
    answers<-gagaOutput@solution[number_of_answers,]
    phylogeny_matrix<-generate_phylo_matrix(answers[1:number_of_clones],number_of_clones)
    phylogeny_matrix<-rbind(phylogeny_matrix, answers[1:number_of_clones])
    #if(contamination==1) {cbind(phylogeny_matrix, c(-1)}
    colnames(phylogeny_matrix)=LETTERS[1:pseudo_number_of_clones]
    
    proportion_matrix<-matrix(answers[(number_of_clones+1):((number_of_cases*pseudo_number_of_clones)+number_of_clones)], 
                              ncol=pseudo_number_of_clones, nrow=number_of_cases, byrow=TRUE)
    colnames(proportion_matrix)=LETTERS[1:pseudo_number_of_clones]    
    for (rows in 1:nrow(proportion_matrix)) {
      scale<-sum(proportion_matrix[rows,])
      proportion_matrix[rows,]<-proportion_matrix[rows,]/scale
    }
    presence_matrix<-matrix(nrow=number_of_mutations, ncol=number_of_clones)
    for (rows in 1:number_of_mutations) {
      presence_matrix[rows,]<-phylogeny_matrix[answers[((number_of_cases*pseudo_number_of_clones)+number_of_clones+rows)],]
    }
    colnames(presence_matrix)=LETTERS[1:pseudo_number_of_clones]    
    
    ## Produce output text files only if outType is complete
    if(outType %in% c("complete")){
      message("Outputting the highest scoring solutions as text files")
      write.table(answers,paste0(output_file_prefix,number_of_clones,"clones.complete.solution",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE,col.names=FALSE)
      write.table(phylogeny_matrix,paste0(output_file_prefix,number_of_clones,"clones.phylogeny.solution.",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE)
      write.table(proportion_matrix,paste0(output_file_prefix,number_of_clones,"clones.proportions.solution.",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE)
      write.table(presence_matrix,paste0(output_file_prefix,number_of_clones,"clones.mutation.matrix.",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE)
      write.table(gagaOutput@best,paste0(output_file_prefix,number_of_clones,"clones.highest.score.solution.",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE,col.names=FALSE)
    }
  
    ## Output fitness convergence plot
    if(outType %in% c("complete","fitness")){
      message("Outputting fitness convergence plot")
      if(outType=="complete"){ png(paste0(output_file_prefix,number_of_clones,"clones.fitnessconvergence.png")) }
      plot(gagaOutput,ylim=yRange)
      if(outType=="complete"){ dev.off() }
    }
      
    ## Produce heatmap
    if(outType %in% c("complete","heatmap")){
      message("Outputting heatmap")
      ## Iterate through presence_matrix data and replace values with hex colours
      clones=presence_matrix
      colnames(clones)=1:ncol(clones) 
      for(i in 1:ncol(clones)){
        ## Use the modulo remainder NOT i, so we can recycle the same colour object
        ## regardless of the dimensions of presence_matrix
        recycler=ifelse(i%%length(colours)==0,length(colours),i%%length(colours))
        clones[clones[,i]=="1",i]=colours[recycler]
        clones[clones[,i]=="0",i]=secondary[4] # background is ICR light grey
      }
      if(outType=="complete"){ png(paste0(output_file_prefix,number_of_clones,"clones.solution",number_of_answers,"of",top_number,".heatmap.png")) }
      heatmap.plus(as.matrix(observation_matrix),Colv=NA,col=heatmapColours,scale="none",Rowv=NULL,labRow=gagaOutput@names,RowSideColors=clones,cexRow=0.9)
      if(outType=="complete"){ dev.off() }
    }  
    
    ## Produce proportion stacked barplot
    if(outType %in% c("complete","proportion")){
      message("Outputting proportion barplot")
      rownames(proportion_matrix)=colnames(observation_matrix)
      if(outType=="complete"){ png(paste0(output_file_prefix,number_of_clones,"clones.solution",number_of_answers,"of",top_number,".proportions.png")) }
      barplot(t(proportion_matrix),col=colours,border=NA,xlab="Sample name",ylab="Proportion of cells")
      if(outType=="complete"){ dev.off() }
    }
    
    ## Produce phylogeny plot
    if(outType %in% c("complete","phylogeny")){
      message("Outputting phylogeny plot")
          
      ## The phylogeny matrices are not compliant with Rgraphviz, so we edit them:
      ## 1.) We remove the bottom row (which is the original encoded string from the solution)
      pm=phylogeny_matrix[1:ncol(phylogeny_matrix),]
      ## 2.) Matrices must have 0 diag, sum(column)<=1, keeping the leftmost value
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
      
      ## Set up plotting parameters
      nAttrs = list() ## per node attributes
      nAttrs$fillcolor = rep(NA,nrow(pm))
      for(i in 1:nrow(pm)){
        ## Use the modulo remainder NOT i, so we can recycle the same colour object
        recycler=ifelse(i%%length(colours)==0,length(colours),i%%length(colours))
        nAttrs$fillcolor[i]=colours[recycler]
      }
      names(nAttrs$fillcolor) = LETTERS[1:nrow(pm)]
      nAttrs$color = rep(NA,nrow(pm))
      names(nAttrs$color) = LETTERS[1:nrow(pm)] ## Assign NA color to shape to ensure clean fill
      
      ## We store the original margins before plotting as stupid RGraphViz screws them up
      old.par <- par(mar = c(0, 0, 0, 0))
      if(outType=="complete"){ png(paste0(output_file_prefix,number_of_clones,"clones.solution",number_of_answers,"of",top_number,".phylogeny.png")) }
      plot(pg, nodeAttrs = nAttrs, attrs=list(edge = list(color=secondary[4],lwd=4))) # "attrs" sets global attributes
      if(outType=="complete"){ dev.off() }
      par(old.par)
    }
    
    #####################################
    ## End plotting / producing output ##
    #####################################
    
  }  
  
  ######################################
  ## End looping over TOP 5 solutions ##
  ######################################
  
}
# /end of gagaReport



