plotPCAmeta=function(GE,metadata,variable,filepath,discretizemethod)
{
  
  
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    
    
    temp=metadata[[i]]
    temp=temp[,variable]
    
    if(variable=="age:"|variable=="Age.at.Diagnosis.in.Days")
    {
      temp=discretizeAge(temp,discretizemethod)
      
      
    }
    
   
    
    
    #print(GE[i])
    #print(temp[,variable])
    x=plotPCA(filepath,GE[[i]],temp,variable)
    my.plots[[i]]=x[1]
    
  }
  
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}

discretizeAge=function(values,method)
{
  
  #if using SLE
  # age=as.numeric(gsub("(age: (([0-9]*[.])?[0-9]+))","\\2",values))
  #if using AML data
  age=values
  print(age)
  #For AML
  discretizedAgeFrequency=discretize(age,method = method, breaks = 3,  labels = c("infants", "child","AYA"))
  #For SLE
  #discretizedAgeFrequency=discretize(age,method = method, breaks = 4,  labels = c("young", "youth","adult","senior"))
  
  return(discretizedAgeFrequency)
}


plotPCA= function(filepath,GeneExpressionMatrixlocal,Groups,variable){
  
  myplots=NULL
  #pdf(filepath)
  # log transform 
  # apply PCA - scale. = TRUE is highly 
  # advisable, but default is FALSE. 
  
  #ir.pca <- prcomp(GeneExpressionMatrixlocal,
  #                center = TRUE,
  #               scale. = TRUE) 
  #GeneExpressionMatrixlocal=GeneExpressionMatrixlocal
  #Groups=decision_Linda2
  #variable="PCA_Adult"
  ir.pca <- prcomp(GeneExpressionMatrixlocal,
                   center = TRUE) 
  
  
  g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
                groups =as.factor(Groups), ellipse = FALSE, 
                circle = FALSE,var.axes = FALSE)
  g <- g + scale_color_discrete(name = variable)
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top',panel.background = element_rect(fill = "lightgray",
                                                                       # colour = "lightgray",
                                                                         size = 0.2, linetype = "solid"),
                 panel.grid.major = element_line(size = 0.2, linetype = 'solid',
                                                 colour = "white"), 
                 panel.grid.minor = element_line(size = 0.2, linetype = 'solid',
                                                 colour = "white")
  )
  #print(g)
  #dev.off()
  plot(g)
  myplots=recordPlot()
  return(list(plots=myplots,pca=ir.pca))
}


plotVenn=function(geneLS,parameters)
{
  
  VENN.LIST <- geneLS
  venn.plot <- venn.diagram(VENN.LIST , NULL, fill=parameters[[1]], alpha=parameters[[2]], cex = as.numeric(parameters[[3]]), cat.fontface=4, category.names=parameters[[5]],cat.pos=5, main="Gene Lists")
  
  # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
  grid.draw(venn.plot)
  
  # To get the list of gene present in each Venn compartment we can use the gplots package
  require("gplots")
  
  a <- venn(VENN.LIST, show.plot=FALSE)
  
  # You can inspect the contents of this object with the str() function
  inters <- attr(a,"intersections")
  return(inters)
  
  
}

plotVenn=function(geneLS,parameters)
{
  
  VENN.LIST <- geneLS
  venn.plot <- venn.diagram(VENN.LIST , NULL, fill=parameters[[1]], alpha=parameters[[2]], cex = as.numeric(parameters[[3]]), cat.fontface=4, category.names=parameters[[5]],cat.pos=5, main="Gene Lists")
  
  # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
  grid.draw(venn.plot)
  
  # To get the list of gene present in each Venn compartment we can use the gplots package
  require("gplots")
  
  a <- venn(VENN.LIST, show.plot=FALSE)
  
  # You can inspect the contents of this object with the str() function
  inters <- attr(a,"intersections")
  return(inters)
  
  
}


plotDistributionCategory=function(GE,titles,variable,filepath,discretizemethod)
{  
  
  print(length(GE))
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    op <- par(mar = c(2,2,2,2) + 0.1)
    end_point = 0.5 + length(GE[[i]]) #+ length(GE[[i]])-1 #this is the line which does the trick (together with barplot "space = 1" parameter)
    
    plt=barplot(table(GE[[i]]),xlab = "",xaxt = "n",space=1,cex.axis = 0.5 ,cex.names = 0.5,cex.main=0.8,main=titles[i])
    par(op)
    text( plt, par("usr")[3]-0.25, 
          srt = 60, adj= 1, xpd = TRUE,
          labels = names(table(GE[[i]])), cex=0.65)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
    
  }
  
  graphics.off()
  # update_geom_default("point", list(size=1))
  # theme_set(theme_grey(base_size=6))
  sc <- c(0.5,0.75,1)
  
  pdf(width=7*sc[1],height=7*sc[1],pointsize=1,paste(filepath,variable,"-",discretizemethod,".pdf",sep=""))
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}

#--------Plotting Functions--------------------
plotDistributionCategory=function(GE,titles,variable,filepath,discretizemethod)
{  
  
  print(length(GE))
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    op <- par(mar = c(16,7,7,4) + 0.1)
    end_point = 0.5 + length(GE[[i]]) #+ length(GE[[i]])-1 #this is the line which does the trick (together with barplot "space = 1" parameter)
    
    plt=barplot(table(GE[[i]]),xlab = "",xaxt = "n",space=1,cex.axis = 0.5 ,cex.names = 0.5,cex.main=0.8,main=titles[i])
    par(op)
    text( plt, par("usr")[3]-0.25, 
          srt = 60, adj= 1, xpd = TRUE,
          labels = names(table(GE[[i]])), cex=0.65)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
    
  }
  
  graphics.off()
  # update_geom_default("point", list(size=1))
  # theme_set(theme_grey(base_size=6))
  sc <- c(0.5,0.75,1)
  
  pdf(width=7*sc[1],height=7*sc[1],pointsize=1,paste(filepath,variable,"-",discretizemethod,".png",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}


#--------Plotting Functions--------------------
plotDistributionCategory=function(GE,titles,variable,filepath,discretizemethod)
{  
  
  print(length(GE))
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    op <- par(mar = c(16,7,7,4) + 0.1)
    end_point = 0.5 + length(GE[[i]]) #+ length(GE[[i]])-1 #this is the line which does the trick (together with barplot "space = 1" parameter)
    
    plt=barplot(table(GE[[i]]),xlab = "",xaxt = "n",space=1,cex.axis = 0.5 ,cex.names = 0.5,cex.main=0.8,main=titles[i])
    par(op)
    text( plt, par("usr")[3]-0.25, 
          srt = 60, adj= 1, xpd = TRUE,
          labels = names(table(GE[[i]])), cex=0.65)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
    
  }
  
  graphics.off()
  # update_geom_default("point", list(size=1))
  # theme_set(theme_grey(base_size=6))
  sc <- c(0.5,0.75,1)
  
  pdf(width=7*sc[1],height=7*sc[1],pointsize=1,paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}
#For TARGET AML project
#Plot clusters Vs Metadata 
plotClustersMetadataPlot=function(ordered_metadata_for_exploratory,clusters,meta,path,date,folder)
{
 # vec = NULL
  myplots <- list()  # new empty list
  for (j in 1:length(meta))
{   for (i in 1:length(clusters))#i in 1:length(meta))
        local({
            i <- i
            theme_set(theme_cowplot(font_size=4)) # reduce default font size
            # The first is for plotting 1 color group for all meta data variables
             # p1 <- qplot(ordered_metadata_for_exploratory[which(ordered_metadata_for_exploratory[,"clusters"]=="red"),meta[i]], as.factor(rownames(ordered_metadata_for_exploratory)[which(ordered_metadata_for_exploratory[,"clusters"]=="red")]),colour = I("red"),xlab = meta[i], ylab = "samples")+theme(plot.background=element_rect(fill="white", colour=NA))#+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
           #The second for plotting all color group for just 1 meta data variable at a time
            
             p1 <- qplot(ordered_metadata_for_exploratory[which(ordered_metadata_for_exploratory[,"clusters"]==as.character(clusters[i])),as.character(meta[j])], as.factor(rownames(ordered_metadata_for_exploratory)[which(ordered_metadata_for_exploratory[,"clusters"]==as.character(clusters[i]))]),colour = I(as.character(clusters[i])),xlab = meta[j], ylab = "samples")+theme(plot.background=element_rect(fill="white", colour=NA))#+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
            
                print(i)
             print(p1)
             myplots[[i]] <<- p1  # add each plot into plot list
         })
   #multiplot(plotlist = myplots ,ncol=2,nrow=6)
  # ggarrange(plotlist =  myplots ,ncol=2,nrow=5)
 # pdf(paste(getwd(),"red",".pdf",sep=""), onefile=TRUE)

  # The solution to the overlaying problem that used to take place in multiple grid was setting the base_aspect_ratio
    dir.create(paste(path,date,"/",folder,"/","MetaData-Clusters",sep = ""))
    
   plot2by2=plot_grid(plotlist = myplots ,labels="auto",nrow=4,ncol=4, label_size=1, rel_widths=c(3,3,3,3),align="hv")
   save_plot(paste(path,date,"/",folder,"/","MetaData-Clusters","/",as.character(meta[j]),"-",date,".png",sep=""), plot2by2,base_aspect_ratio=1.5
           #  ncol = 2, # we're saving a grid plot of 2 columns
           #  nrow = 2, # and 2 rows
             # each individual subplot should have an aspect ratio of 1.3
            # base_aspect_ratio = 1.3
   )
  }
}
plotEnrichment=function(enriched,titleGO,temp)
{
  
print("here")
View(enriched)
  gostplot(enriched, capped = TRUE, interactive = TRUE)
  

#  pp <- publish_gostplot(enriched, width = NA, height = NA, filename = temp )
  
 # my.plots <- vector(length(GE), mode='list')
  
#  graphics.off()
 
 # pdf(temp)
 # dotplot(enriched, showCategory=30,title=titleGO)
 # dev.off()
  
}


plotGWAS=function(GWASTable,variable,filepath,discretizemethod=NULL)
{  
  
  #print(length(GE))
  
  my.plots <- vector(length(unique(GWASTable$CHR)), mode='list')
  
  chromosomes=unique(GWASTable$CHR)
  for(i in 1:length(chromosomes)){
    manhattan(GWASTable[GWASTable$CHR==as.character(chromosomes[i]),], main = paste("chromosome",as.character(chromosomes[i])), annotatePval = 0.01)
    myplots=recordPlot()
    my.plots[[i]]=myplots
    
  }
  
  graphics.off()
  # update_geom_default("point", list(size=1))
  # theme_set(theme_grey(base_size=6))
 # sc <- c(0.5,0.75,1)
  
  #pdf(width=7*sc[1],height=7*sc[1],pointsize=1,paste(filepath,variable,"-",discretizemethod,".pdf",sep=""))
  pdf(file=paste(filepath,"/",variable,"-",discretizemethod,".pdf",sep=""))
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}


plotGeneRulesEnrichment=function(RosettaEnrichment,path,date,folder)
{
  
  #RosettaEnrichment=TARGET_ReleventBP
  temp=getAnnotatedGenes(genes,as.character(RosettaEnrichment$features))
  print(RosettaEnrichment$features)
  print(unlist(RosettaEnrichment$features))
  print(temp)
  temp=annotateInOrder(temp,unlist(RosettaEnrichment$features))
 RosettaEnrichment$features=temp
  BP=RosettaEnrichment$BP
  p1 <-qplot(as.character(RosettaEnrichment$features), as.factor(as.character(BP)),xlab = "Feature", ylab = "Biological Process")+theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1),plot.background=element_rect(fill="white", colour=NA))#+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  save_plot(paste(path,date,"/",folder,"/","Gene-RulesEnrichment","-",date,".png",sep=""), p1,base_aspect_ratio=1.5)
            
}
