library("HMMcopy")

args <- commandArgs(TRUE)

tumour_wig <- args[1]
map <- args[2]
gc <- args[3]
output <- args[4]
output_pdf <- args[5]
mappability <- args[6]
sampleid <- args[7]

pdf(output_pdf)


x <- wigsToRangedData(tumour_wig, gc, map)

x$valid <- TRUE
x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
x$ideal <- TRUE
routlier <- 0.01
range <- quantile(x$reads[x$valid], prob = c(0, 1 - routlier), na.rm = TRUE)
doutlier <- 0.001
domain <- quantile(x$gc[x$valid], prob = c(doutlier, 1 - doutlier), na.rm = TRUE)
x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] | x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE

ideal <- subset(x, x$ideal)

rough = loess(ideal$reads ~ ideal$gc, span = 0.03)




i <- seq(0, 1, by = 0.001)
final <- loess(predict(rough, i) ~ i, span = 0.3)
#corrected
cor_reads <- predict(loess(predict(rough, i) ~ i, span = 0.3))

#plot start
reads <- predict(rough, i)
gc <- i
df <- data.frame(reads, gc)
df <- df[order(gc),]
df <- df[complete.cases(df), ]

df$cor_reads <-  cor_reads

plot(ideal$reads, ideal$gc, type='p', col=rgb(1, 0, 0, 0.5), main=paste(sampleid, 'loess', sep=' '))
lines(df$cor_reads, df$gc)
#plot stop


x$cor.gc <- x$reads / predict(final, x$gc)


coutlier <- 0.01
range <- quantile(x$cor.gc[which(x$valid)], prob = c(0, 1 - coutlier), na.rm = TRUE)

x$map_keep <- (x$cor.gc < range[2])

y <- subset(x, x$map_keep)

final <- approxfun(lowess(y$map, y$cor.gc))


print (head(y$map))
print (head(y$cor.gc))
print (head(final(y$map)))
 
x$cor.map <- x$cor.gc / final(x$map)
x$copy <- x$cor.map
x$copy[x$copy <= 0] = NA
x$copy <- log(x$copy, 2)

y <- as.data.frame(x)
names(y)[names(y) == 'space'] <- 'chromosome'
y$map_keep <- NULL

write.table(y, file=output, sep=',', quote=FALSE, row.names=FALSE)


dev.off()
