### Mean and variance

if (grepl("/src",getwd())){setwd("..")}
source('src/pp_nets.R')

## color scaling
c_scale <- colorRamp(grey(seq(0.01,0.75,length = 10)))

### group by site, date and trt
grp.sdt <- paste(info[,"site"],info[,"date"],info[,"trt"],sep = "_")
ppn.sdt <- split(pp.nets,grp.sdt)
ppn.prb <- list()
ppn.avg <- list()
ppn.std <- list()
ppn.info <- do.call(rbind,strsplit(names(ppn.sdt),split = "_"))
colnames(ppn.info) <- c("site","date","trt")

for (i in 1:length(ppn.sdt)){
    x <- sign(ppn.sdt[[i]][[1]])
    for (j in 2:length(ppn.sdt[[i]])){
        x <- x + sign(ppn.sdt[[i]][[j]])
    }
    ppn.prb[[i]] <- x / length(ppn.sdt[[i]])
}
for (i in 1:length(ppn.sdt)){
    x <- ppn.sdt[[i]][[1]]
    for (j in 2:length(ppn.sdt[[i]])){
        x <- x + ppn.sdt[[i]][[j]]
    }
    ppn.avg[[i]] <- x / length(ppn.sdt[[i]])
}
for (i in 1:length(ppn.sdt)){
    ss <- list()
    for (j in 1:length(ppn.sdt[[i]])){
        ss[[j]] <- (ppn.sdt[[i]][[j]] - ppn.avg[[i]])^2
    }
    for (j in 2:length(ss)){
        ss[[1]] <- ss[[1]] + ss[[j]]
    }
    ppn.std[[i]] <- (sqrt(ss[[1]])) / sqrt(length(ppn.sdt[[i]]) - 1)
}

site <- "hawley"
igp <- ppn.prb
ppn.igp.prb <- list(trt1 = list(
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-07-20" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-07-27" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-08-10" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-08-24" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-09-07" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-09-14" & ppn.info[,"site"] == site]),
                trt2 = list(
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-07-20" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-07-27" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-08-10" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-08-24" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-09-07" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-09-14" & ppn.info[,"site"] == site]),
                trt3 = list(
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-07-20" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-07-27" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-08-10" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-08-24" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-09-07" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-09-14" & ppn.info[,"site"] == site])
                )
igp <- ppn.avg
ppn.igp.avg <- list(trt1 = list(
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-07-20" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-07-27" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-08-10" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-08-24" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-09-07" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == "1999-09-14" & ppn.info[,"site"] == site]),
                trt2 = list(
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-07-20" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-07-27" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-08-10" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-08-24" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-09-07" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == "1999-09-14" & ppn.info[,"site"] == site]),
                trt3 = list(
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-07-20" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-07-27" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-08-10" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-08-24" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-09-07" & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == "1999-09-14" & ppn.info[,"site"] == site])
                )
ppn.dat <- c("1999-07-20","1999-07-27","1999-08-10","1999-08-24","1999-09-07","1999-09-14")

tmp.lay <- as.matrix(read.csv('data/tmp_lay.csv'))

pdf("results/pp_prob_time_Hawley.pdf",width = 26, height = 15)
par(mfrow = c(3,length(ppn.dat)),mar = c(0,0,0,3))
for (i in 1:length(ppn.igp.prb)){
    for (j in 1:length(ppn.igp.prb[[i]])){
        g <- ppn.igp.prb[[i]][[j]][[1]]
        vs.lab <- sign((apply(g,1,sum) + apply(g,2,sum)))
        vs <- vs.lab * 20
        vs.lab <- vs.lab * 2
        g <- g / max(g)
        rownames(g) <- colnames(g) <- substr(rownames(g),1,2)
        tmp <- graph_from_adjacency_matrix(g,weighted = TRUE)
        tmp.prb <- ppn.igp.avg[[i]][[j]][[1]] / max(ppn.igp.avg[[i]][[j]][[1]])
        tmp.prb <- graph_from_adjacency_matrix(tmp.prb,weighted = TRUE)
        E(tmp)$color = apply(c_scale(E(tmp)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
        V(tmp)$color = grey(0.5,alpha = 0.1)
        plot.igraph(tmp,layout = tmp.lay, 
                    vertex.label.color = 'black', 
                    vertex.label.cex = vs.lab + 0.000001,
                    vertex.label.font = 1,
                    vertex.size = vs + 0.000001,
                    edge.arrow.size = 1.25,
                    edge.width = (10 * E(tmp.prb)$weight))
        if (i == 3){legend('bottom', legend = ppn.dat[j], bty = 'n', cex = 2.5)}
        if (j == 1){legend('topleft', legend = paste0("Trt ",i), bty = 'n', cex = 3)}
    }
}
dev.off()

site <- "Molly"
igp <- ppn.prb
ppn.dat <- c("1999-07-20","1999-07-27","1999-08-03","1999-08-10","1999-08-17","1999-08-25")
ppn.igp.prb <- list(trt1 = list(
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[1] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[2] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[3] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[4] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[5] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[6] & ppn.info[,"site"] == site]),
                trt2 = list(
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[1] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[2] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[3] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[4] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[5] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[6] & ppn.info[,"site"] == site]),
                trt3 = list(
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[1] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[2] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[3] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[4] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[5] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[6] & ppn.info[,"site"] == site])
                )
igp <- ppn.avg
ppn.igp.avg <- list(trt1 = list(
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[1] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[2] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[3] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[4] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[5] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 1 & ppn.info[,"date"] == ppn.dat[6] & ppn.info[,"site"] == site]),
                trt2 = list(
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[1] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[2] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[3] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[4] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[5] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 2 & ppn.info[,"date"] == ppn.dat[6] & ppn.info[,"site"] == site]),
                trt3 = list(
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[1] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[2] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[3] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[4] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[5] & ppn.info[,"site"] == site],
                    igp[ppn.info[,"trt"] == 3 & ppn.info[,"date"] == ppn.dat[6] & ppn.info[,"site"] == site])
                )
tmp.lay <- as.matrix(read.csv('data/tmp_lay.csv'))

pdf("results/pp_prob_time_Molly.pdf",width = 26, height = 15)
par(mfrow = c(3,length(ppn.dat)),mar = c(0,0,0,3))
for (i in 1:length(ppn.igp.prb)){
    for (j in 1:length(ppn.igp.prb[[i]])){
        g <- ppn.igp.prb[[i]][[j]][[1]]
        vs.lab <- sign((apply(g,1,sum) + apply(g,2,sum)))
        vs <- vs.lab * 20
        vs.lab <- vs.lab * 2
        g <- g / max(g)
        rownames(g) <- colnames(g) <- substr(rownames(g),1,2)
        tmp <- graph_from_adjacency_matrix(g,weighted = TRUE)
        tmp.prb <- ppn.igp.avg[[i]][[j]][[1]] / max(ppn.igp.avg[[i]][[j]][[1]])
        tmp.prb <- graph_from_adjacency_matrix(tmp.prb,weighted = TRUE)
        E(tmp)$color = apply(c_scale(E(tmp)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
        V(tmp)$color = grey(0.5,alpha = 0.1)
        plot.igraph(tmp,layout = tmp.lay, 
                    vertex.label.color = 'black', 
                    vertex.label.cex = vs.lab + 0.000001,
                    vertex.label.font = 1,
                    vertex.size = vs + 0.000001,
                    edge.arrow.size = 1.25,
                    edge.width = (10 * E(tmp.prb)$weight))
        if (i == 3){legend('bottom', legend = ppn.dat[j], bty = 'n', cex = 2.5)}
        if (j == 1){legend('topleft', legend = paste0("Trt ",i), bty = 'n', cex = 3)}
    }
}
dev.off()
