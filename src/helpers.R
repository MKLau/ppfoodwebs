# Construct a supraadjacency matrix from a list

### Populate the Mouquet nitrogen and carbon models
### using biomass from the observed pitcher plants

## Multiplot function
##
## ggplot objects can be passed in ..., 
## or to plotlist (as a list of ggplot objects)
## - cols:   Number of columns in layout
## - layout: A matrix specifying the layout. If present, 'cols' is ignored.
##
## If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
## then plot 1 will go in the upper left, 2 will go in the upper right, and
## 3 will go all the way across the bottom.
##
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    ## Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    ## If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        ## Make the panel
        ## ncol: Number of columns of plots
        ## nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        ## Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        ## Make each plot, in the correct location
        for (i in 1:numPlots) {
            ## Get the i,j matrix positions of the regions 
            ## that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
        }
    }
}

#function to count rotifers, removes "other" and "rotifer volume" and naming
count.and.name <- function(x, zero.na = TRUE){
    var.set <- c("Fletcherimyia",
                 "Habrotrocha",
                 "headcap",
                 "Metriocnemus",
                 "rotvol ml",
                 "Sarraceniopus",
                 "volume ml",
                 "Wyeomyia")
    if (!("variable" %in% colnames(x))){
        warning('Variable column missing')
    }else{
        if (!all(var.set %in% x[,"variable"])){
            missing <- var.set[!(var.set %in% x[,"variable"])]
            add <- x[rep(1,length(missing)),]
            add[,"value"] <- 0
            add[,"variable"] <- missing
            x <- rbind(x,add)
            x <- x[match(var.set,x[,"variable"]),]
            rownames(x) <- 1:nrow(x)
        }else{}
        x <- x[match(var.set,x[,"variable"]),]
        ## replaces rotifer counts per volume with total counts
        x[x[,"variable"] == "Habrotrocha" , "value"] <- 
            x[x[,"variable"] == "Habrotrocha" , "value"] * 
                x[x[,"variable"] == "volume ml" , "value"] / 
                    x[x[,"variable"] == "rotvol ml" , "value"] 
        ## rename variable for easier referencing
        levels(x$variable) <- c("Ff", "Hr", "Prey", "Mk", "rotvol", 
                                "Sg", "Volume", "Ws")
        ## Remove rotvol as it's now redundant with Habrotrocha
        x <- subset(x, variable != "rotvol")
        ## Make NA's zero if present, which could happen
        ## When there is no rotvol
        if (zero.na & is.na(x[x[,"variable"] == "Hr","value"])){
            x[x[,"variable"] == "Hr","value"] <- 0
        }else{}
        return(x) #returns the matrix output
    }
}

#function to match names to interaction matrix order
match.mat <- function(x,interaction.names = 
                          c("Ff", "Ws", "Hr", "Mk", 
                            "Sg", "Prey", "Volume")){
  x <- x[match(interaction.names, x$variable),]
  return(x)
}  

## Make the pitcher plant food web based on the Mouquet 2008 model
mkFoodWeb <- function(x = "variable value matrix",
                      mod = "flow model matrix",
                      mod.par = "param value matrix",
                      Volume.ml = TRUE,
                      zero.na = TRUE) {
    ## Rename variables to conform to flow model
    x[,"variable"] <- c("fle", "hab", "met", "det", 
                        "sar", "Volume", "wye")
    x[is.na(x[,"value"]),"value"] <- 0
    ## re-format the 
    flow <- array(NA,dim = dim(mod))
    rownames(flow) <- colnames(flow) <- colnames(mod)
    ## create objects for each of the parameters
    for (i in 1:nrow(mod.par)){
        assign(as.character(mod.par[i,"name"]),
               as.numeric(as.character(mod.par[i,"value"])))
    }
    ## create objects for the variables
    for (i in 1:nrow(x)){
        assign(as.character(x[i,"variable"]), 
               as.numeric(as.character(x[i,"value"])))
    }
    ## convert ml to L 
    if (Volume.ml){Volume <- Volume / 1000} 
    ## Build the flow matrix based on the model and observed values
    for (i in 1:nrow(mod)){
        for (j in 1:ncol(mod)){
            flow[i,j] <- eval(parse(text = as.character(mod[i,j])))
        }
    }
    if (zero.na){flow[is.na(flow)] <- 0}
    flow
}

### Interplote based on previous timestep
interp.obs <- function(data.split.2, obs){
    ## Fit interpolated values for NAs using ordered time
    out <- list()
    for (i in 1:length(unique(obs))){
        x <- data.split.2[obs == unique(obs)[i]]
        x <- do.call(rbind,x)
        var <- as.character(x[,'variable'])
        x <- split(x,var)
        for (j in 1:length(x)){
            date <- as.Date(x[[j]][,'date'],format = "%Y-%m-%d")    
            x[[j]] <- x[[j]][order(date),]
            if (any(is.na(x[[j]][,'value']))){
                for (k in 1:length(x[[j]][,'value'])){
                    if (is.na(x[[j]][k,'value'])){
                        ## value is zero if initial and NA
                        if (k == 1){x[[j]][k,'value'] <- 0}
                        ## value is the preceding value if NA
                        x[[j]][k,'value'] <- mean(x[[j]][k-1,'value'])
                    }else{}
                }
            }else{}
        }
        x <- do.call(rbind,x)
        rownames(x) <- 1:nrow(x)
        x <- as.data.frame(x)
        x[,'value'] <- as.numeric(as.character(x[,'value']))
        out[[i]] <- x
    }
    out
}
