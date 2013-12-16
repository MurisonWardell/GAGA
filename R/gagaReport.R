#' View and create output of GAGA objects
#' 
#' Description here
#' 
#' @param gagaInput             Raw data put into gaga() 
#' @param gagaOutput            Object of class ga produced by gaga()
#' @param yRange                Y-axis range when plotting fitness of individuals
#' @return Description of the returned object
#' @export
# @seealso \code{\link{fermat.test}}  #### other functions; e.g. gaga report output
# @references GAGA paper here!
#' @author Alex Murison and Christopher Wardell \email{Alexander.Murison@@icr.ac.uk}
#' @examples
#' gagaReport(blah,UNKNOWN)


## Writes output 
gagaReport<-function(gagaInput,gagaOutput,yRange=c(-250,0)) {
  
  ## Set observation matrix for heatmap
  observation_matrix = gagaInput
  
  ## Get number of mutations from gagaOutput
  number_of_mutations=length(gagaOutput@names)
  
  ## Get number of clones from gagaOutput
  number_of_clones=gagaOutput@min
  
  ## Get number of cases (e.g. timepoints) from gagaOutput
  number_of_cases=gagaOutput@max
  
  ## We can now work out if contamination was present using the two above values and the length of solution
  if(((length(gagaOutput@solution)-number_of_clones-number_of_mutations)/number_of_cases)==number_of_clones){
    contamination=0
  }else{
    contamination=1
  }
    
  if (contamination == 1) {
    pseudo_number_of_clones<-number_of_clones+1
    #numBits<-(pseudo_number_of_clones*number_of_cases)+(number_of_mutations+number_of_clones)    
  } else {
  #  numBits<-(number_of_clones*number_of_cases)+(number_of_mutations+number_of_clones)
    pseudo_number_of_clones<-number_of_clones
  }
  
  
  
  
  # Set ICR Colour Scheme - user configurable?
  buylrd=c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026") # blue yellow red 
  primary=c("#C9DD03","#FFD602","#F9A100","#EE7EA6","#A71930","#616365") # primary ICR colours
  secondary=c("#726E20","#6E273D","#F00034","#ADAFAF","#003D4C") # secondary ICR colours
  colours=c(primary,secondary[c(1,2,3,5)]) # exclude ICR light grey
  
  ## Load annotation data
  #snvs=read.table(file=annotations_file,header=TRUE,stringsAsFactors=FALSE)
  snvs=gagaOutput@names
  print("Outputting the highest scoring solutions")
  
  ## Prints results if there is more than one top solution - should we just take the first one instead?
  top_number<-min(5,nrow(gagaOutput@solution))
  for (number_of_answers in 1:top_number) {
    #best_result<-matrix(gagaOutput@solution[number_of_answers,], ncol=number_of_clones, byrow=TRUE)
    answers<-gagaOutput@solution[number_of_answers,]
    phylogeny_matrix<-generate_phylo_matrix(answers[1:number_of_clones],number_of_clones)
    proportion_matrix<-matrix(answers[(number_of_clones+1):((number_of_cases*pseudo_number_of_clones)+number_of_clones)], 
                              ncol=pseudo_number_of_clones, nrow=number_of_cases, byrow=TRUE)
    phylogeny_matrix<-rbind(phylogeny_matrix, answers[1:number_of_clones])
    #if(contamination==1) {cbind(phylogeny_matrix, c(-1)}
    
    for (rows in 1:nrow(proportion_matrix)) {
      scale<-sum(proportion_matrix[rows,])
      proportion_matrix[rows,]<-proportion_matrix[rows,]/scale
    }
    presence_matrix<-matrix(nrow=number_of_mutations, ncol=number_of_clones)
    for (rows in 1:number_of_mutations) {
      presence_matrix[rows,]<-phylogeny_matrix[answers[((number_of_cases*pseudo_number_of_clones)+number_of_clones+rows)],]
    }
    #write.table(answers,paste(output_file_prefix, ".",number_of_clones,"clones.complete.solution.",number_of_answers,"rep", loop_value, ".txt", sep=""),sep="\t", row.names=FALSE)
    #write.table(phylogeny_matrix,paste(output_file_prefix, ".",number_of_clones,"clones.phylogeny.solution.",number_of_answers,"rep", loop_value, ".txt", sep=""),sep="\t", row.names=FALSE)
    #write.table(proportion_matrix,paste(output_file_prefix, ".",number_of_clones,"clones.proportions.solution.",number_of_answers,"rep", loop_value, ".txt", sep=""),sep="\t", row.names=FALSE)
    #write.table(presence_matrix,paste(output_file_prefix, ".",number_of_clones,"clones.mutation.matrix.",number_of_answers,"rep", loop_value, ".txt", sep=""),sep="\t", row.names=FALSE)
    
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
  
  
  ## Produce heatmap
  heatmap.plus(as.matrix(observation_matrix),Colv=NA,col=buylrd,scale="none",Rowv=NULL,labRow=gagaOutput@names,RowSideColors=clones,cexRow=0.9)
  props = proportion_matrix
  rownames(props)=colnames(observation_matrix)
  barplot(t(props),col=colours,border=NA,xlab="Sample name",ylab="Proportion of cells")
  #frame()
  #lim<-par()
  #rasterImage(lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
    
  
  #  print("Outputting Convergence plot")
  #png(paste(output_file_prefix, ".",number_of_clones,"rep", loop_value, ".convergence.png", sep=""))
  plot(gagaOutput,ylim=yRange)
  #dev.off()
  
  #write.table(gagaOutput@best,paste(output_file_prefix, ".",number_of_clones,"clones.", loop_value, ".reps.highest.score.solution.",number_of_answers,".txt", sep=""),sep="\t", row.names=FALSE)
  
  
  print(summary(gagaOutput))
  
}
# /end of gagaReport


# ################################
# ## Phylogeny plots begin here ##
# ################################
# 
# number_of_clones=16
# 
# ## Create new phylogeny
# p=new_phylo(0)
# 
# ## Get phylogeny matrix from AM's function
# pm=generate_phylo_matrix(p)
# 
# ## These matrices are not compliant with Rgraphviz, so we edit them:
# ## Matrices must have 0 diag, sum(column)<=1, keep leftmost value
# diag(pm)=0
# for(i in 1:nrow(pm)){
#   for(j in ncol(pm):1){
#     if(sum(pm[,j])>1){
#       pm[i,j]=0
#     }
#   }
# }
# ## Set lineage names and cast to graphNEL object
# rownames(pm)=LETTERS[1:nrow(pm)]
# colnames(pm)=rownames(pm)
# pg=as(pm,"graphNEL")
# 
# ## Set up plotting parameters and set colours
# primary=c("#C9DD03","#FFD602","#F9A100","#EE7EA6","#A71930","#616365") # primary ICR colours
# secondary=c("#726E20","#6E273D","#F00034","#ADAFAF","#003D4C") # secondary ICR colours
# colours=c(primary,secondary[c(1,2,3,5)]) # exclude ICR light grey
# nAttrs = list() ## per node attributes
# nAttrs$fillcolor = rep(NA,nrow(pm))
# ## Iterate through nAttrs$fillcolor and replace values with hex colours
# for(i in 1:nrow(pm)){
#   ## Use the modulo remainder NOT i, so we can recycle the same colour object
#   recycler=ifelse(i%%length(colours)==0,length(colours),i%%length(colours))
#   nAttrs$fillcolor[i]=colours[recycler]
# }
# names(nAttrs$fillcolor) = LETTERS[1:nrow(pm)]
# nAttrs$color = rep(NA,nrow(pm))
# names(nAttrs$color) = LETTERS[1:nrow(pm)] ## Assign NA color to shape to ensure clean fill
# plot(pg, nodeAttrs = nAttrs, attrs=list(edge = list(color=secondary[4],lwd=4))) # "attrs" sets global attributes
# 
# ##############################
# ## Phylogeny plots end here ##
# ##############################


