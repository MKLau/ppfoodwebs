### Network modelling and analysis of 
### Hawley and Molly bog data from AME and NJG
### Collected 1999-2000

if (grepl("/src",getwd())){setwd("..")}

## Installing packages
packs <- c("Rgraphviz", "Hmisc", "enaR", "ggplot2", "igraph",
           "intergraph", "reshape", "grid")

if (!("Rgraphviz" %in% installed.packages()[,1])){
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rgraphviz")
}

if (!("Hmisc" %in% installed.packages()[,1])){
    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("sjmgarnier/viridis")
    install.packages("Hmisc", repos = "http://cloud.r-project.org/")
}

sapply(packs[-1:-2], function(x, y) 
    if (!(x %in% y)){install.packages(x, 
                       repos='http://cran.us.r-project.org')}, 
    y = installed.packages()[,1])
sapply(packs, require, character.only = TRUE)

### Load the network models
source('src/load_ppnets.R')

### Check for results directory
dir.create('results',showWarnings = FALSE)
## remove networks with no pitcher present
## Base this on the treatments
if (nrow(info) == 3101){
    rm.nets <- !((info[,"trt"] == 2) & (info[,"date"] < "1999-06-29") | 
                     (info[,"trt"] == 3) & (info[,"date"] < "1999-07-20"))
### trt 4
### remove Volume is zero or NA before the "open" date
    trt4.dates <- sort(unique(info[info["trt"] == 4,"date"]))
    rm.trt4 <- !(info[,"trt"] == 4 & 
                     (((info[,"leaf"] == 2) & (info[,"date"] < trt4.dates[2])) | 
                          ((info[,"leaf"] == 3) & (info[,"date"] < trt4.dates[6])) | 
                              ((info[,"leaf"] == 4) & (info[,"date"] < trt4.dates[9])) | 
                                  ((info[,"leaf"] == 5) & (info[,"date"] < trt4.dates[12])) 
                      )
                 )
    pp.nets <- pp.nets[rm.nets & rm.trt4]
    info <- info[rm.nets & rm.trt4,]
    storage <- out[rm.nets & rm.trt4]
}
### Convert to network model type
if (!"pp.ena" %in% ls()){
    pp.ena <- list()
    for (i in 1:length(pp.nets)){
        x <- pp.nets[[i]]
        inp <- x["ant",]
        res <- x[,"respiration"]
        flo <- x[!rownames(x) %in% c("respiration","ant"),!colnames(x) %in% c("respiration","ant")]
        exp <- rep(0,nrow(flo))
        sto <- rep(0,nrow(flo))
        output <- exp
        pp.ena[[i]] <- pack(flo,inp,res,exp,output,sto,living = c(FALSE,FALSE,TRUE,TRUE,TRUE,TRUE))
    }
}
## Generate a vector of months
month <- substr(info[,"date"],6,7)
## Generate a vector of sampling times
time <- as.Date(info[,'date'],format = '%Y-%m-%d')
### Structural analysis
if (!"pp.ns" %in% ls()){
    pp.ns <- list()
    for (i in 1:length(pp.ena)){
        if (any(is.na(as.matrix(pp.ena[[i]])))){
            pp.ns[[i]] <- rep(NA,13)
        }else{
            pp.ns[[i]] <- enaStructure(pp.ena[[i]])$ns
        }
    }
    pp.ns <- do.call(rbind,pp.ns)
}
### Ascendency
if (!"pp.asc" %in% ls()){
    pp.asc <- lapply(pp.ena,enaAscendency)
    pp.asc <- do.call(rbind,pp.asc)
}
### NOTE: TSTp is based on an adjusted AMI, which has an small value added to 
### avoid dividing by zero.
if (!"TSTp" %in% colnames(pp.asc)){
    pp.asc <- cbind(pp.asc,TSTp = pp.asc[,'ASC'] / (pp.asc[,'AMI'] + 1e-10))
}

### Basic network structural characteristics
pdf("results/pp_structure.pdf",width = 15, height = 15)
pairs(pp.ns[,c("n","L","C","LD","ppr","no.scc","no.scc.big")],pch = 19, cex = 0.5)
dev.off()


### Ascendency through time with variance
haw.asc <- data.frame(info[info[,'site'] == 'hawley',-2],
                date = as.Date(info[info[,'site'] == 'hawley','date'],format="%Y-%m-%d"),
                pp.asc[info[,'site'] == 'hawley',])
mol.asc <- data.frame(info[info[,'site'] == 'Molly',-2],
                date = as.Date(info[info[,'site'] == 'Molly','date'],format="%Y-%m-%d"),
                pp.asc[info[,'site'] == 'Molly',])
## divide overhead by TSTp
haw.asc$rOH <- haw.asc$OH / haw.asc$TSTp
mol.asc$rOH <- mol.asc$OH / mol.asc$TSTp
## zero NA when CAP and AMI are 0
haw.asc[(haw.asc[,"AMI"] == 0 & haw.asc[,"CAP"] == 0),c("ASC.CAP","OH.CAP","rOH","robustness")] <- 0
haw.asc[is.na(haw.asc[,"robustness"]),"robustness"] <- 0
mol.asc[(mol.asc[,"AMI"] == 0 & mol.asc[,"CAP"] == 0),c("ASC.CAP","OH.CAP","rOH","robustness")] <- 0
mol.asc[is.na(mol.asc[,"robustness"]),"robustness"] <- 0

### Combine all network metrics
pp.dat <- cbind(rbind(haw.asc,mol.asc),pp.ns)
if (!(all(pp.dat[,c("trt", "plant", "leaf")] == info[,c("trt", "plant", "leaf")]))){
    warning("Data are corrupted. Make sure to check pp.dat.")
}
### Make sure to deal with both NAs and Inf
sum(na.inf <- apply(pp.dat,1,function(x) any(is.na(x)) | any(x == Inf) | any(x == -Inf)))
### For now, I'm making any NA or Inf or -Inf equal to 0
pp.dat[is.na(pp.dat) | pp.dat == Inf | pp.dat == -Inf] <- 0
sum(na.inf <- apply(pp.dat,1,function(x) any(is.na(x)) | any(x == Inf) | any(x == -Inf)))

### Examining correlations among network stats
pp.cor <- cor(pp.dat[,6:ncol(pp.dat)])
diag(pp.cor) <- NA
pp.cor <- pp.cor[!(apply(pp.cor, 1, function(x) all(is.na(x)))), 
                 !(apply(pp.cor, 2, function(x) all(is.na(x))))]
diag(pp.cor) <- 1

### Figures without treatment 4
     ## ASC: Returns the ascendnecy of a network.  Ascendency is a scaled
     ##      form of AMI relative to the total system throughput
     ##      (Ulanowicz 1997; 2004).  Total system throughput is the sum
     ##      of all activity in a network (Kay et al. 1989).
     ## ASC = TSTp * AMI
gg.asc.h <- ggplot(haw.asc[haw.asc[,"trt"] != 4,],aes(x=date,y=ASC,color = factor(trt)))+geom_point() +
    stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
        guides(color=FALSE)
gg.asc.m <- ggplot(mol.asc[mol.asc[,"trt"] != 4,],aes(x=date,y=ASC, color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
     ## TSTp = Total system throughput is the sum
     ## of all activity in a network (Kay et al. 1989).
gg.TSTp.h <- ggplot(haw.asc[haw.asc[,"trt"] != 4,],aes(x=date,y=TSTp,color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
      guides(color=FALSE)
gg.TSTp.m <- ggplot(mol.asc[mol.asc[,"trt"] != 4,],aes(x=date,y=TSTp, color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
     ## AMI: Returns the Average Mutual Information (AMI) in a network.
     ##      AMI provides a measure of the constraints placed on a given
     ##      peice of energy matter moving through a network (Patricio et
     ##      al. 2006)
gg.AMI.h <- ggplot(haw.asc[haw.asc[,"trt"] != 4,],aes(x=date,y=AMI,color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
      guides(color=FALSE)
gg.AMI.m <- ggplot(mol.asc[mol.asc[,"trt"] != 4,],aes(x=date,y=AMI, color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
      ## OH: Returns the overhead of a network.  Overhead is the
      ##     proportion of the capacity in a network that is not used as
      ##     ascendency (Ulanowicz 2004).
gg.OH.h <- ggplot(haw.asc[haw.asc[,"trt"] != 4,],aes(x=date,y=OH,color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt')
gg.OH.m <- ggplot(mol.asc[mol.asc[,"trt"] != 4,],aes(x=date,y=OH, color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt')

gg.rOH.h <- ggplot(haw.asc[haw.asc[,"trt"] != 4,],
                   aes(x=date,y=rOH,color = factor(trt)))+ geom_point() +
                       stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt')
gg.rOH.m <- ggplot(mol.asc[mol.asc[,"trt"] != 4,],
                   aes(x=date,y=rOH, color = factor(trt)))+ geom_point() +
                       stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt')


     ## ELD: Returns the Effective Link Density of the network(c)
     ##      (Ulanowicz et al. 2014).
gg.ELD.h <- ggplot(haw.asc[haw.asc[,"trt"] != 4,],aes(x=date,y=ELD,color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
      guides(color=FALSE)
gg.ELD.m <- ggplot(mol.asc[mol.asc[,"trt"] != 4,],aes(x=date,y=ELD, color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
      ## TD: Returns the Trophic Depth of the network(r) (Ulanowicz et al.
      ##     2014).
gg.TD.h <- ggplot(haw.asc[haw.asc[,"trt"] != 4,],aes(x=date,y=TD,color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
      guides(color=FALSE)
gg.TD.m <- ggplot(mol.asc[mol.asc[,"trt"] != 4,],aes(x=date,y=TD, color = factor(trt)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
      ## State space of ELD and TD
gg.ss.h <- ggplot(haw.asc[haw.asc[,"trt"] != 4,],aes(x=TD,y=ELD,color = factor(trt)))+geom_point() +
##  stat_summary(fun.data=mean_cl_boot, geom="smooth") + 
      labs(color = 'Hawley trt')
gg.ss.m <- ggplot(mol.asc[mol.asc[,"trt"] != 4,],aes(x=TD,y=ELD, color = factor(trt)))+geom_point() +
##  stat_summary(fun.data=mean_cl_boot, geom="smooth") + 
      labs(color = 'Molly trt')


### Treatment 4: use pitcher number as factor
     ## ASC: Returns the ascendnecy of a network.  Ascendency is a scaled
     ##      form of AMI relative to the total system throughput
     ##      (Ulanowicz 1997; 2004).  Total system throughput is the sum
     ##      of all activity in a network (Kay et al. 1989).
     ## ASC = TSTp * AMI
gg.asc.h.4 <- ggplot(haw.asc[haw.asc[,"trt"] == 4,],aes(x=date,y=ASC,color = factor(leaf)))+geom_point() +
    stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
        guides(color=FALSE)
gg.asc.m.4 <- ggplot(mol.asc[mol.asc[,"trt"] == 4,],aes(x=date,y=ASC, color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
     ## TSTp = Total system throughput is the sum
     ## of all activity in a network (Kay et al. 1989).
gg.TSTp.h.4 <- ggplot(haw.asc[haw.asc[,"trt"] == 4,],aes(x=date,y=TSTp,color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
      guides(color=FALSE)
gg.TSTp.m.4 <- ggplot(mol.asc[mol.asc[,"trt"] == 4,],aes(x=date,y=TSTp, color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
     ## AMI: Returns the Average Mutual Information (AMI) in a network.
     ##      AMI provides a measure of the constraints placed on a given
     ##      peice of energy matter moving through a network (Patricio et
     ##      al. 2006)
gg.AMI.h.4 <- ggplot(haw.asc[haw.asc[,"trt"] == 4,],aes(x=date,y=AMI,color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
      guides(color=FALSE)
gg.AMI.m.4 <- ggplot(mol.asc[mol.asc[,"trt"] == 4,],aes(x=date,y=AMI, color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
      ## OH: Returns the overhead of a network.  Overhead is the
      ##     proportion of the capacity in a network that is not used as
      ##     ascendency (Ulanowicz 2004).
gg.OH.h.4 <- ggplot(haw.asc[haw.asc[,"trt"] == 4,],aes(x=date,y=OH,color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt')
gg.OH.m.4 <- ggplot(mol.asc[mol.asc[,"trt"] == 4,],aes(x=date,y=OH, color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt')

gg.rOH.h.4 <- ggplot(haw.asc[haw.asc[,"trt"] == 4,],
                     aes(x=date,y=rOH,color = factor(leaf)))+ geom_point() +
                         stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt')
gg.rOH.m.4 <- ggplot(mol.asc[mol.asc[,"trt"] == 4,],
                     aes(x=date,y=rOH, color = factor(leaf)))+ geom_point() +
                         stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt')


     ## ELD: Returns the Effective Link Density of the network(c)
     ##      (Ulanowicz et al. 2014).
gg.ELD.h.4 <- ggplot(haw.asc[haw.asc[,"trt"] == 4,],aes(x=date,y=ELD,color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
      guides(color=FALSE)
gg.ELD.m.4 <- ggplot(mol.asc[mol.asc[,"trt"] == 4,],aes(x=date,y=ELD, color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
      ## TD: Returns the Trophic Depth of the network(r) (Ulanowicz et al.
      ##     2014).
gg.TD.h.4 <- ggplot(haw.asc[haw.asc[,"trt"] == 4,],aes(x=date,y=TD,color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Hawley trt') +
      guides(color=FALSE)
gg.TD.m.4 <- ggplot(mol.asc[mol.asc[,"trt"] == 4,],aes(x=date,y=TD, color = factor(leaf)))+geom_point() +
  stat_summary(fun.data=mean_cl_boot, geom="smooth") + labs(color = 'Molly trt') +
      guides(color=FALSE)
      ## State space of ELD and TD
gg.ss.h.4 <- ggplot(haw.asc[haw.asc[,"trt"] == 4,],aes(x=TD,y=ELD,color = factor(leaf)))+geom_point() +
##  stat_summary(fun.data=mean_cl_boot, geom="smooth") + 
      labs(color = 'Hawley trt')
gg.ss.m.4 <- ggplot(mol.asc[mol.asc[,"trt"] == 4,],aes(x=TD,y=ELD, color = factor(leaf)))+geom_point() +
##  stat_summary(fun.data=mean_cl_boot, geom="smooth") + 
      labs(color = 'Molly trt')




## Values vs AMI
Volume <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Volume","value"]))
Volume[is.na(Volume)] <- 0
gg.Volume.ami <- ggplot(data.frame(info, Volume = Volume,pp.asc)[info[,"trt"] != 4,],
                  aes(x = Volume,y = AMI, color = factor(trt))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(trt), color = factor(trt)))  +
                          guides(color=FALSE)
Ff <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Ff","value"]))
Ff[is.na(Ff)] <- 0
gg.Ff.ami <- ggplot(data.frame(info, Ff = Ff,pp.asc)[info[,"trt"] != 4,],
                  aes(x = Ff,y = AMI, color = factor(trt))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(trt), color = factor(trt)))  +
                          guides(color=FALSE)
Hr <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Hr","value"]))
Hr[is.na(Hr)] <- 0
gg.Hr.ami <- ggplot(data.frame(info, Hr = Hr,pp.asc)[info[,"trt"] != 4,],
                  aes(x = Hr,y = AMI, color = factor(trt))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(trt), color = factor(trt))) +
                          guides(color=FALSE)
Mk <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Mk","value"]))
Mk[is.na(Mk)] <- 0
gg.Mk.ami <- ggplot(data.frame(info, Mk = Mk,pp.asc)[info[,"trt"] != 4,],
                  aes(x = Mk,y = AMI, color = factor(trt))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(trt), color = factor(trt))) +
                          guides(color=FALSE)
Sg <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Sg","value"]))
Sg[is.na(Sg)] <- 0
gg.Sg.ami <- ggplot(data.frame(info, Sg = Sg,pp.asc)[info[,"trt"] != 4,],
                  aes(x = Sg,y = AMI, color = factor(trt))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(trt), color = factor(trt))) +
                          guides(color=FALSE)
Prey <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Prey","value"]))
Prey[is.na(Prey)] <- 0
gg.Prey.ami <- ggplot(data.frame(info, Prey = Prey,pp.asc)[info[,"trt"] != 4,],
                  aes(x = Prey,y = AMI, color = factor(trt))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(trt), color = factor(trt))) +
                          guides(color=FALSE)
Ws <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Ws","value"]))
Ws[is.na(Ws)] <- 0
gg.Ws.ami <- ggplot(data.frame(info, Ws = Ws,pp.asc)[info[,"trt"] != 4,],
                  aes(x = Ws,y = AMI, color = factor(trt))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(trt), color = factor(trt))) 
## trt 4
Volume <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Volume","value"]))
Volume[is.na(Volume)] <- 0
gg.Volume.ami4 <- ggplot(data.frame(info, Volume = Volume,pp.asc)[info[,"trt"] == 4,],
                  aes(x = Volume,y = AMI, color = factor(leaf))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(leaf), color = factor(leaf)))  +
                          guides(color=FALSE)
Ff <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Ff","value"]))
Ff[is.na(Ff)] <- 0
gg.Ff.ami4 <- ggplot(data.frame(info, Ff = Ff,pp.asc)[info[,"trt"] == 4,],
                  aes(x = Ff,y = AMI, color = factor(leaf))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(leaf), color = factor(leaf)))  +
                          guides(color=FALSE)
Hr <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Hr","value"]))
Hr[is.na(Hr)] <- 0
gg.Hr.ami4 <- ggplot(data.frame(info, Hr = Hr,pp.asc)[info[,"trt"] == 4,],
                  aes(x = Hr,y = AMI, color = factor(leaf))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(leaf), color = factor(leaf))) +
                          guides(color=FALSE)
Mk <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Mk","value"]))
Mk[is.na(Mk)] <- 0
gg.Mk.ami4 <- ggplot(data.frame(info, Mk = Mk,pp.asc)[info[,"trt"] == 4,],
                  aes(x = Mk,y = AMI, color = factor(leaf))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(leaf), color = factor(leaf))) +
                          guides(color=FALSE)
Sg <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Sg","value"]))
Sg[is.na(Sg)] <- 0
gg.Sg.ami4 <- ggplot(data.frame(info, Sg = Sg,pp.asc)[info[,"trt"] == 4,],
                  aes(x = Sg,y = AMI, color = factor(leaf))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(leaf), color = factor(leaf))) +
                          guides(color=FALSE)
Prey <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Prey","value"]))
Prey[is.na(Prey)] <- 0
gg.Prey.ami4 <- ggplot(data.frame(info, Prey = Prey,pp.asc)[info[,"trt"] == 4,],
                  aes(x = Prey,y = AMI, color = factor(leaf))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(leaf), color = factor(leaf))) +
                          guides(color=FALSE)
Ws <- unlist(lapply(storage,function(x) x[x[,"variable"] == "Ws","value"]))
Ws[is.na(Ws)] <- 0
gg.Ws.ami4 <- ggplot(data.frame(info, Ws = Ws,pp.asc)[info[,"trt"] == 4,],
                  aes(x = Ws,y = AMI, color = factor(leaf))) + geom_point() + 
                      geom_smooth(method = "lm", se=FALSE, 
                                  aes(group = factor(leaf), color = factor(leaf))) 


### 
pdf("results/ami_trt_time.pdf",width = 21, height = 9)
multiplot(gg.TSTp.h,gg.TSTp.m,
          gg.AMI.h,gg.AMI.m,
          gg.rOH.h,gg.rOH.m,
          cols=3)
dev.off()


pdf("results/amiVvalues.pdf",width = 20, height = 15)
multiplot(gg.Volume.ami,
          gg.Prey.ami,
          gg.Hr.ami,
          gg.Mk.ami,
          gg.Sg.ami,
          gg.Ff.ami,
          gg.Ws.ami,
          cols=3)
dev.off()

pdf("results/amiVvalues_trt4.pdf",width = 20, height = 15)
multiplot(gg.Volume.ami4,
          gg.Prey.ami4,
          gg.Hr.ami4,
          gg.Mk.ami4,
          gg.Sg.ami4,
          gg.Ff.ami4,
          gg.Ws.ami4,
          cols=3)
dev.off()

pdf("results/asc_time.pdf",width = 21, height = 9)
multiplot(gg.asc.h,gg.asc.m,
          gg.TSTp.h,gg.TSTp.m,
          gg.AMI.h,gg.AMI.m,
          gg.OH.h,gg.OH.m,
          cols=4)
dev.off()


pdf("results/ami_trt_time_trt4.pdf",width = 21, height = 9)
multiplot(gg.TSTp.h.4,gg.TSTp.m.4,
          gg.AMI.h.4,gg.AMI.m.4,
          gg.rOH.h.4,gg.rOH.m.4,
          cols=3)
dev.off()

pdf("results/eldtd_statespace.pdf",width = 21, height = 9)
multiplot(gg.ELD.h,gg.ELD.m,
          gg.TD.h,gg.TD.m,
          gg.ss.h,gg.ss.m,
          cols = 3)
dev.off()

pdf("results/eldtd_statespace_trt4.pdf",width = 21, height = 9)
multiplot(gg.ELD.h.4,gg.ELD.m.4,
          gg.TD.h.4,gg.TD.m.4,
          gg.ss.h.4,gg.ss.m.4,
          cols = 3)
dev.off()

system("cp results/ami* docs/manuscript")
system("cp results/pp_* docs/manuscript")

### Window of Vitality 
if (!"ena.asc" %in% ls()){
    data(enaModels)
    ena.asc <- do.call(rbind,lapply(enaModels,enaAscendency))
}

wov <- rbind(ena.asc,pp.asc[,-24])
wov <- data.frame(wov, model = c(rep('ena',nrow(ena.asc)),
             rep('pp',nrow(pp.asc[,-24]))))
gg.wov <- ggplot(wov,aes(x=TD,y=ELD,color = factor(model))) + geom_point() + 
      labs(color = 'Model Type')

pdf("results/wov_enaVpp.pdf",height = 9,width = 9)
gg.wov
dev.off()

### Heatmap of correlations among network variables
pdf("results/net_corheatmap.pdf", height = 9, width = 9)
heatmap(pp.cor)
dev.off()

### Multi-layer network analysis

### Base Carbon model plot
c.mod <- mod
c.mod[c.mod != "0"] <- 1
c.mod <- apply(c.mod,2,as.numeric)
rownames(c.mod) <- colnames(c.mod)
rgv <- igraph:::as_graphnel(graph_from_adjacency_matrix(c.mod))

pdf('results/pp_c_model.pdf')
par(mfrow = c(1,1))
plot(rgv,attrs = list(edge = list(arrowsize = 0.5)))
dev.off()

### cp results to manuscript directory
system("cp results/*.pdf docs/manuscript")

