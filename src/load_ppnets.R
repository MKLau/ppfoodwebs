### Makes the models of the pitcher plant food webs
### based on the data collected by AM Ellison and NJ 
### Gotelli.

### Models generated from empirical abundances of inquilines 
### and non-living compartments (e.g. water level and prey inputs)

### Data from the Harvard Forest Data Archive hf193

packs <- c("igraph", "sna")
sapply(packs[-1], function(x, y) 
    if (!(x %in% y)){install.packages(x, 
         repos = 'http://cran.us.r-project.org')}, 
    y = installed.packages()[, 1])
sapply(packs, require, character.only = TRUE)

source('helpers.R')

## 1999 dataset is HF193-01
## hf193-01-hm99-rs.csv
rawdata <- read.csv("../data/hf193-01-hm99-rs.csv")

## remove value of "." for "value" and reclass from factor to numeric
rawdata[which(rawdata$value=="."), "value"] <- NA
rawdata$value <- as.numeric(as.character(rawdata$value))
rawdata[which(rawdata$date=="1999-70-20"), "date"] <-"1999-07-20"

## Fix issue with Molly_1999-09-10_1_ 5_1, rotvol ml is equal to 0
## Aaron determined that this should be 0.5
rawdata[apply(rawdata[, 1:5], 1, paste, collapse="_") == "Molly_1999-09-10_1_ 5_1" 
       & rawdata[, "variable"] == "rotvol ml", "value"] <- 0.5
## 
rawdata<-rawdata[rawdata$variable != "other", ]
rawdata <- droplevels(rawdata)

## Pitcher plant food web model based on Mouquet et al. 2008
## The current model differs mainly in terms of:
## 1. Detritus and water volume are known
## 2. Bacteria are unknown and combined into deteritus
## 3. Several additional fly and one mite species have been introduced
## All flows are units of mg C
mod.par <- na.omit(read.csv("../data/mouquet_model_2008.csv"))
mod <- na.omit(read.csv("../data/pp_C_flow.csv")[, -1])
mod <- as.matrix(mod)
colnames(mod)[colnames(mod) == "prey"] <- "det"
rownames(mod) <- colnames(mod)

## split into a list of observations by individual pitchers in time
data.split <- split(rawdata, (apply(rawdata[, 1:5] , 1, paste, collapse="_")))

### counting rotifers, removing "other" and "rotifer volume" and naming
data.split.1 <- lapply(data.split, count.and.name)

### match names to interaction matrix order
data.split.2 <-lapply(data.split.1, match.mat)

### Create the information data.frame
info <- do.call(rbind, strsplit(names(data.split.2), "_"))
colnames(info) <- c('site', 'date', 'trt', 'plant', 'leaf')
obs <- apply(info[, c('site', 'date', 'trt', 'plant', 'leaf')], 1, paste, collapse="_")

### Interplote missing values based on the assumption
### that there is 100% survival and/or 100% replacement
out <- interp.obs(data.split.2, obs)

### Re-order by variables
out <- do.call(rbind, out)
out <- split(out, apply(out[, 1:5], 1, paste, collapse="_"))

### Make array conformable
for (i in 1:length(out)){
    if (nrow(out[[i]]) != nrow(data.split.2[[i]])){
        x <- as.matrix(out[[i]])
        na.mat <- matrix(NA, 
                         nrow=abs(nrow(x) - nrow(data.split.2[[i]])), 
                         ncol = ncol(x))
        out[[i]] <- rbind(na.mat, x)
        rownames(out[[i]]) <- 1:nrow(out[[i]])
    }else{}
}
out <- lapply(out, as.data.frame)

### Get interaction matrix for foodwebs
pp.nets <- lapply(out, mkFoodWeb, 
                  mod = mod, mod.par = mod.par)

### final output info 
pp.info <- do.call(rbind, strsplit(names(pp.nets), "_"))
pp.info <- as.data.frame(pp.info)
colnames(pp.info) <- c('site', 'date', 'trt', 'plant', 'leaf')
obs <- apply(pp.info[, c('site', 'plant', 'leaf')], 1, paste, collapse="_")
pp.info[, "date"] <- as.Date(pp.info[, "date"], format = "%Y-%m-%d")

