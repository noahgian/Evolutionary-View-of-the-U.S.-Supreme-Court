#----Libraries----#

library(MASS)
library(ape)
library(phytools)
library(adegenet)
library(softImpute)
library(prodlim)

#----Read and clean data----#

#Read in main data matrix of cases (update file version and working directory as needed)
bigdat <- read.csv(file="SCDB_2020_01_justiceCentered_Citation.csv",header=TRUE)

#Remove initials from justice names
levels(bigdat$justiceName) <- list("Fortas"="AFortas","Goldberg"="AJGoldberg","Kennedy"="AMKennedy","Scalia"="AScalia","Kavanaugh"="BMKavanaugh","White"="BRWhite","Whittaker"="CEWhittaker","Thomas"="CThomas","Souter"="DHSouter","Kagan"="EKagan","Warren"="EWarren","Frankfurter"="FFrankfurter","Murphy"="FMurphy","Vinson"="FMVinson","Blackmun"="HABlackmun","Burton"="HHBurton","Black"="HLBlack","Roberts"="JGRoberts","Harlan"="JHarlan2","Stevens"="JPStevens","Powell"="LFPowell","Gorsuch"="NMGorsuch","Stewart"="PStewart","Ginsburg"="RBGinsburg","Jackson"="RHJackson","Alito"="SAAlito","O'Connor"="SDOConnor","Reed"="SFReed","Breyer"="SGBreyer","Minton"="SMinton","Sotomayor"="SSotomayor","Clark"="TCClark","Marshall"="TMarshall","Rutledge"="WBRutledge","Burger"="WEBurger","Rehnquist"="WHRehnquist","Brennan"="WJBrennan","Douglas"="WODouglas")

#Remove the cases that are missing values for the "majority" variable
badCases=bigdat[which(is.na(bigdat[,"majority"])),"caseId"]
badRows=which(bigdat[,"caseId"] %in% badCases)
newdat=bigdat[-badRows,] 

#make vector of natural courts that have all nine justices
fullcourtnums <- c()
for(i in 1:length(natcourtnums)){
  if (length(unique(newdat[which(newdat$naturalCourt == natcourtnums[i]),]$justiceName)) == 9){fullcourtnums <- c(fullcourtnums,natcourtnums[i])}
}

#---Misc functions for extracting court/case info----#

#input: natural court number, case matrix
#output: starting and concluding year (as length 2 integer vector)
#example: 
# > natCourtYears(1406,bigdat)
# [1] 1958 1961
natCourtYears <- function(nat,mycases){
  X <- mycases[which(mycases$naturalCourt == nat),]$term
  return(c(min(X),max(X)))
}

#input: ordinal number of case in natural court, natural court number, case matrix
#output: list of case info: justice names, majority/minority votes, case name
#note: assumes 9 justices in the natural court
#example: caseinfo(1,1406,newdat)
caseinfo <- function(m,natcourt,mycases){
  B <- mycases[which(mycases$naturalCourt == natcourt),]
  newcases <- B[(1:9)+9*(m-1),]
  return(list(newcases$justiceName,newcases$majority,unique(newcases$caseName)))
}


#----Functions related to the dissimilarity matrix----#

#input: integer vector of natural court number(s), case matrix
#output: dissimilarity matrix
#example: findSimsLong(1406,newdat)
findSimsLong=function(natCourt,mycases){
  #extract case ID, justice names, and majority votes from cases in specified natural court(s)
  myCourt=mycases[mycases[,"naturalCourt"] %in% natCourt,c("caseId","justiceName","majority")]
  #remove cases with missing majority vote
  badCases=myCourt[which(is.na(myCourt[,"majority"])),"caseId"]
  badRows=which(myCourt[,"caseId"] %in% badCases)
  if (length(badRows) == 0){newCourt <- myCourt}
  if (length(badRows) > 0){newCourt=myCourt[-badRows,]}
  total=nrow(newCourt)
  #get names of all justices in relevant court(s)
  justiceList=unique(as.matrix(newCourt[,"justiceName"]))
  #get case IDs of relevant cases
  caseList=unique(as.matrix(newCourt[,"caseId"]))
  #compute similarity matrix
  voteSims=matrix(0, ncol=length(justiceList), length(justiceList))
  rownames(voteSims)=justiceList
  colnames(voteSims)=justiceList
  voteCount=matrix(0, ncol=length(justiceList), length(justiceList))
  rownames(voteCount)=justiceList
  colnames(voteCount)=justiceList
  for(case in caseList){
    currentCase <- newCourt[newCourt$caseId == case,]
    for(i in currentCase$justiceName){
      for(j in currentCase$justiceName){
        voteCount[matrix(i),matrix(j)] <- voteCount[matrix(i),matrix(j)]+1
        if(currentCase[which(currentCase$justiceName == i),"majority"]==currentCase[which(currentCase$justiceName == j),"majority"]){
          voteSims[i,j]=voteSims[i,j]+1
        }  
      }
    }
  }
  #convert similarity to dissimilarity
  return(1-(voteSims/voteCount))
}

#compute the dissimilarity matrix for the entire data set
natcourtnums <- unique(newdat$naturalCourt)
fulldiss <- findSimsLong(natcourtnums,newdat)

#input: natural court number, full dissimilarity matrix of all justices/cases
#output: sub-dissimilarity matrix given by justices in the specified natural court
#note: assumes 9 justices in the natural court
#example: dissMax(1707,fulldiss)
dissMax <- function(natcourt,fulldiss){
  return(fulldiss[as.matrix(unique(newdat[which(newdat$naturalCourt==natcourt),]$justiceName)),as.matrix(unique(newdat[which(newdat$naturalCourt==natcourt),]$justiceName))])
}

#----Phylogenetic tree functions----#

#input: dissimilarity matrix, optional boolean "flip" to reverse color palette, optional "rotate" for plot, option fontsize, optional "majaxis" bolding of major axis 
#output: plots OLS tree with tips colored by 1-dim MDS
#example: phylotreeplot(dissMax(1707,fulldiss),flip=F,rotate=70,fontsize=0.7)
phylotreeplot <- function(dis, flip=F, rotate=0, fontsize = 0.6,majaxis=F){
  myPal <- colorRampPalette(c("red","blue"))
  mds <- cmdscale(dis,k=1)
  if(flip==T){mds <- -mds}
  tree <- optim.phylo.ls(dis, stree=nj(dis))
  tree$edge.length <- round(tree$edge.length,digits=3)*100
  edgewidth <- rep(1,length(tree$edge.length))
  if (majaxis==TRUE)
  {
    nodes <- which(dis == max(dis),arr.ind=TRUE)[1,]
    path <- nodepath(tree,from=nodes[1],to=nodes[2])
    pathedge <- rep(0,length(path)-1)
    for (i in 1:(length(path)-1)){
      pathedge[i] <- max(row.match(path[(i+1):i],tree$edge),row.match(path[i:(i+1)],tree$edge),na.rm=TRUE)
    }
    edgewidth[pathedge] <- 3
  }
  plot(tree, type="unrooted", edge.width=edgewidth, rotate.tree=rotate, show.tip=FALSE, no.margin=TRUE)
  tiplabels(tree$tip.label, col="white",bg=num2col(mds[tree$tip.label,], col.pal=myPal), cex=fontsize)
  edgelabels(tree$edge.length, bg="white", col="darkgray",cex=0.6,frame="none",adj=c(0.5,-0.55))
}

#---Make plots for paper----#

#1707
pdf("Dropbox/Shared work/mca-1103472/1707.pdf",onefile = TRUE, height=4, width=7.5)
phylotreeplot(dissMax(1707,fulldiss),flip=F,rotate=70,fontsize=0.7, majaxis=T)
dev.off()

#1704
pdf("Dropbox/Shared work/mca-1103472/1704.pdf",onefile = TRUE, height=4, width=7.5)
phylotreeplot(dissMax(1704,fulldiss),flip=F,rotate=290,fontsize=0.7, majaxis=T)
dev.off()

#1606
pdf("Dropbox/Shared work/mca-1103472/1606.pdf",onefile = TRUE, height=4, width=7.5)
phylotreeplot(dissMax(1606,fulldiss),flip=F,rotate=350,fontsize=0.7, majaxis=T)
dev.off()

#1506
pdf("Dropbox/Shared work/mca-1103472/1506.pdf",onefile = TRUE, height=4, width=7.5)
phylotreeplot(dissMax(1506,fulldiss),flip=T,rotate=75,fontsize=0.7, majaxis=T)
dev.off()

#1410
pdf("Dropbox/Shared work/mca-1103472/1410.pdf",onefile = TRUE, height=6, width=6)
phylotreeplot(dissMax(1410,fulldiss),flip=T,rotate=165,fontsize=0.7, majaxis=T)
dev.off()

#1403
pdf("Dropbox/Shared work/mca-1103472/1403.pdf",onefile = TRUE, height=3, width=7)
phylotreeplot(dissMax(1403,fulldiss),flip=T,rotate=220,fontsize=0.7, majaxis=T)
dev.off()

#1301
pdf("Dropbox/Shared work/mca-1103472/1301.pdf",onefile = TRUE, height=4, width=6)
phylotreeplot(dissMax(1301,fulldiss),flip=F,rotate=155,fontsize=0.7, majaxis=T)
dev.off()

#1mds
y <- rep(0,9)
mds <- -cmdscale(dissMax(1707,fulldiss),k=1)
pdf("Dropbox/Shared work/mca-1103472/1mds.pdf",onefile = TRUE, width=8, height=4)
plot(data.frame(mds,y),type='n',xlim=c(-0.25,0.25), ylim=c(-0.1,0.1), xlab="", ylab="", asp=0.2, yaxt='n')
text(mds,y,labels=rownames(mds), cex=0.7, srt=90)
dev.off()

#hist
D <- dissMax(1707,fulldiss)
mds <- cmdscale(D,k=2)
Dmds1 <- as.matrix(dist(mds[,1]))
Dmds2 <- as.matrix(dist(mds))
Dtree <- cophenetic.phylo(optim.phylo.ls(D, stree=nj(D)))
sum((function(x)x*x)(100*(D-Dmds1)))
sum((function(x)x*x)(100*(D-Dmds2)))
sum((function(x)x*x)(100*(D-Dtree)))
pdf("Dropbox/Shared work/mca-1103472/hist.pdf",onefile = TRUE)
hist(100*abs(D-Dmds1), col=rgb(0,1,0,0.6),xlim=c(0,100*0.2), ylim=c(0,45), xlab="Absolute value of residuals")
hist(100*abs(D-Dmds2), col=rgb(1,0,0,0.6), add=T)
hist(100*abs(D-Dtree), col=rgb(0,0,1,0.6), add=T)
box()
dev.off()

#residuals
res1 <- rep(0,length(fullcourtnums))
res2 <- rep(0,length(fullcourtnums))
restree <- rep(0,length(fullcourtnums))
max1 <- rep(0,length(fullcourtnums))
max2 <- rep(0,length(fullcourtnums))
maxtree <- rep(0,length(fullcourtnums))
for (i in 1:length(fullcourtnums)){
  D <- dissMax(fullcourtnums[i],fulldiss)
  mds <- cmdscale(D,k=2)
  Dmds1 <- as.matrix(dist(mds[,1]))
  Dmds2 <- as.matrix(dist(mds))
  Dtree <- cophenetic.phylo(optim.phylo.ls(D, stree=nj(D)))
  res1[i] <- sum((function(x)x*x)(100*(D-Dmds1)))
  res2[i] <- sum((function(x)x*x)(100*(D-Dmds2)))
  restree[i] <- sum((function(x)x*x)(100*(D-Dtree)))
  max1[i] <- max(100*abs(D-Dmds1))
  max2[i] <- max(100*abs(D-Dmds2))
  maxtree[i] <- max(100*abs(D-Dtree))
}
xyrs <- c(1946,1949,1953,1955,1956,1957,1958,1962,1965,1967,1970,1972,1975,1981,1986,1988,1990,1991,1993,1994,2005,2006,2009,2010,2017,2018) 
res1plot <- loess(res1~xyrs)
res2plot <- loess(res2~xyrs)
restreeplot <- loess(restree~xyrs)
max1plot <- loess(max1~xyrs)
max2plot <- loess(max2~xyrs)
maxtreeplot <- loess(maxtree~xyrs)
pdf("Dropbox/Shared work/mca-1103472/residuals1.pdf",onefile = TRUE, height=4, width=5)
plot(xyrs,res1, ylim=c(390,16000),col='green', ylab='Sum of residual squares',xlab='First year of natural court')
points(xyrs,res2,pch=2, col='red')
points(xyrs,restree,pch=3, col='blue')
lines(xyrs,predict(res1plot), col='green', lwd=2)
lines(xyrs,predict(res2plot), col='red', lwd=2)
lines(xyrs,predict(restreeplot), col='blue', lwd=2)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/residuals2.pdf",onefile = TRUE, height=4, width=5)
plot(xyrs,max1, ylim=c(5,32),col='green', ylab='Maximum absolute residual',xlab='First year of natural court')
points(xyrs,max2,pch=2, col='red')
points(xyrs,maxtree,pch=3, col='blue')
lines(xyrs,predict(max1plot), col='green', lwd=2)
lines(xyrs,predict(max2plot), col='red', lwd=2)
lines(xyrs,predict(maxtreeplot), col='blue', lwd=2)
dev.off()
length(which(res2 > restree))/length(fullcourtnums)
length(which(max2 > maxtree))/length(fullcourtnums)

#2-dim MDS plot
mdsPlot <- function(courtmds,marg=0.5){
  plot(courtmds,type="n",xlim=c(0,marg),ylim=c(0,marg),xlab="",ylab="",cex.axis=0.75)
  for (i in 1:9){
    text(courtmds[i,1],courtmds[i,2],labels=rownames(courtmds)[i], cex=0.5)
  }
}

#Break into issue areas
fulldississue <- list()
for (i in 1:14){
  fulldississue[[i]] <- findSimsLong(natcourtnums,newdat[which(newdat$issueArea == i),])
}
pdf("Dropbox/Shared work/mca-1103472/1707-1.pdf",onefile = TRUE, height=3, width=5)
phylotreeplot(dissMax(1707,fulldississue[[1]]),flip=F,rotate=345,fontsize=0.9, majaxis=T)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/1707-2.pdf",onefile = TRUE, height=3, width=5)
phylotreeplot(dissMax(1707,fulldississue[[2]]),flip=F,rotate=0,fontsize=0.9, majaxis=T)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/1707-3.pdf",onefile = TRUE, height=3, width=9)
phylotreeplot(dissMax(1707,fulldississue[[3]]),flip=F,rotate=135,fontsize=0.9, majaxis=T)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/1707-4.pdf",onefile = TRUE, height=3, width=10)
phylotreeplot(dissMax(1707,fulldississue[[4]]),flip=F,rotate=0,fontsize=0.9, majaxis=T)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/1707-5.pdf",onefile = TRUE, height=3, width=10)
phylotreeplot(dissMax(1707,fulldississue[[5]]),flip=T,rotate=200,fontsize=0.9, majaxis=T)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/1707-8.pdf",onefile = TRUE, height=3, width=10)
phylotreeplot(dissMax(1707,fulldississue[[8]]),flip=F,rotate=300,fontsize=0.9, majaxis=T)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/1707-9.pdf",onefile = TRUE, height=3, width=10)
phylotreeplot(dissMax(1707,fulldississue[[9]]),flip=F,rotate=110,fontsize=0.9, majaxis=T)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/1707-10.pdf",onefile = TRUE, height=3, width=10)
phylotreeplot(dissMax(1707,fulldississue[[10]]),flip=F,rotate=250,fontsize=0.9, majaxis=T)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/1707-12.pdf",onefile = TRUE, height=3, width=10)
phylotreeplot(dissMax(1707,fulldississue[[12]]),flip=F,rotate=120,fontsize=0.9, majaxis=T)
dev.off()

#evolution of the Roberts 4 Court
add_label_legend <- function(pos = "bottomleft", label, ...) {
  legend(pos, label, bty = "n", ...)
}
A <- newdat[which(newdat$naturalCourt == 1704),]
n <- floor(dim(A)[1]/9)
treelist <- vector("list", length(n))
treechange <- c(1)
for (i in 2:n){
  treelist[[i]] <- fastme.bal(findSimsLong(1704,A[1:(9*i),]))
  if (all.equal.phylo(treelist[[i]],treelist[[i-1]],use.edge.length=FALSE) == FALSE){
    treechange <- c(treechange,i)
    plot(treelist[[i]],type='unrooted',cex=0.5)
    print(paste(i,n))
  }
}
treechange <- treechange[-1]
for (i in 2:length(treelist)){
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Roberts")] <- "    Roberts    "
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Scalia")] <- "    Scalia    "
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Kennedy")] <- "    Kennedy    "
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Thomas")] <- "    Thomas    "
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Ginsburg")] <- "    Ginsburg    "
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Breyer")] <- "    Breyer    "
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Alito")] <- "    Alito    "
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Sotomayor")] <- "    Sotomayor    "
  treelist[[i]]$tip.label[which(treelist[[i]]$tip.label == "Kagan")] <- "    Kagan    "
}
par(mfrow=c(4,10))
for (i in treechange){
  plot(treelist[[i]],type='unrooted',cex=0.5,no.margin=TRUE,lab4ut='a')
}
treechange2 <- c(1,4,6,8,9,12,13,17,20,21,25,29)
par(mfrow=c(4,3))
for (i in treechange2){
  plot(treelist[[treechange[i]]],type='unrooted',cex=0.6,no.margin=TRUE,lab4ut='a')
  print(treechange[i])
  print(caseinfo(treechange[i],1704))
}
treechange2 <- c(1,4,6,8,9,12,13,17,20,21,25,29)
pdf("Dropbox/Shared work/mca-1103472/evol1704.pdf",onefile = NA)
par(mfrow=c(4,3),oma = c(10, 10, 10, 10))
plot(treelist[[treechange[treechange2[1]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=120,x.lim=c(-0.8,0.8))
add_label_legend("bottomleft", paste0("(", letters[1], ")"))
plot(treelist[[treechange[treechange2[2]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=80,x.lim=c(-0.5,0.5))
add_label_legend("bottomleft", paste0("(", letters[2], ")"))
plot(treelist[[treechange[treechange2[3]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=200,x.lim=c(-0.4,0.4),y.lim=c(-0.4,0.4))
add_label_legend("bottomleft", paste0("(", letters[3], ")"))
plot(treelist[[treechange[treechange2[4]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=60,x.lim=c(-0.25,0.25),y.lim=c(-0.25,0.25))
add_label_legend("bottomleft", paste0("(", letters[4], ")"))
plot(treelist[[treechange[treechange2[5]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=80,y.lim=c(-0.4,0.4))
add_label_legend("bottomleft", paste0("(", letters[5], ")"))
plot(treelist[[treechange[treechange2[6]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=100,x.lim=c(-0.25,0.25),y.lim=c(-0.25,0.25))
add_label_legend("bottomleft", paste0("(", letters[6], ")"))
plot(treelist[[treechange[treechange2[7]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=80,x.lim=c(-0.4,0.4),y.lim=c(-0.4,0.4))
add_label_legend("bottomleft", paste0("(", letters[7], ")"))
plot(treelist[[treechange[treechange2[8]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=220,x.lim=c(-0.4,0.4),y.lim=c(-0.4,0.4))
add_label_legend("bottomleft", paste0("(", letters[8], ")"))
plot(treelist[[treechange[treechange2[9]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=90,x.lim=c(-0.4,0.4),y.lim=c(-0.4,0.4))
add_label_legend("bottomleft", paste0("(", letters[9], ")"))
plot(treelist[[treechange[treechange2[10]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=170,x.lim=c(-0.4,0.4),y.lim=c(-0.4,0.4))
add_label_legend("bottomleft", paste0("(", letters[10], ")"))
plot(treelist[[treechange[treechange2[11]]]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=90,x.lim=c(-0.5,0.5),y.lim=c(-0.5,0.5))
add_label_legend("bottomleft", paste0("(", letters[11], ")"))
plot(treelist[[length(treelist)]],no.margin=TRUE,type='unrooted',cex=0.6,lab4ut='a',rotate.tree=85,x.lim=c(-0.5,0.5),y.lim=c(-0.5,0.5))
add_label_legend("bottomleft", paste0("(", letters[12], ")"))
dev.off()

# -- function for plotting a matrix -- #
myImagePlot <- function(x, ...){
     x[which(x == 'NaN')] <- 1
     x <- x*100
     min <- min(x)
     max <- 95 #max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,0.95,length=256),  # Red
                   seq(0,0.95,length=256),  # Green
                   seq(0,0.95,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 par(mar = c(4,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
 axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, las=3, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las=1, cex.axis=0.7)

 # Color Scale
 par(mar = c(4,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="Dissimilarity rate",ylab="",
      xaxt="n")

 layout(1)
}
# -- END matrix plot function -- #

#matrix completion
comp <- complete(fulldiss, softImpute(fulldiss,type="svd",maxit=100000))
comp[comp < 0] <- 0
pdf("Dropbox/Shared work/mca-1103472/imlong.pdf",onefile = TRUE, height=5, width=7)
myImagePlot(fulldiss)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/imcomplong.pdf",onefile = TRUE, height=5, width=7)
myImagePlot(comp)
dev.off()
pdf("Dropbox/Shared work/mca-1103472/treecomplong.pdf",onefile = TRUE, height=6, width=7)
phylotreeplot(comp, flip=T, rotate=50, fontsize = 0.5, majaxis=T)
dev.off()
