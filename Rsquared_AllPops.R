rm(list=ls())

version 2

#Starting fresh to make a script to calculate LD (R^2) for SNP genotype data from populations.

setwd("~/Dropbox/ChromosomeOmy5/Omy5_survey2014/LD_86SNPs")

# Start with 'toolkit' formatted file (.csv), first column called 'Sample_ID'.

x <- read.csv(file="~/Dropbox/ChromosomeOmy5/Omy5_survey2014/Omy5survey2014_86SNPs.csv", na.strings=c("0"))


# Remove the Omy28 markers (R000200 and R018277)

drops <- c("R000200", "R000200.1", "R018277", "R018277.1")
x <- x[,!(names(x) %in% drops)]


######                                                  #############
###### NOTE SUPER IMPORTANT TO DO, REMOVE:              #############
######   INDIVIDUALS WHITH MISSING DATA (HIGH MISSERS)  #############
######                                                  #############

#head(x),  str(x)

# now, make the first column be a Population column, by stripping the numbers off the Individuals (use appropriate y for your infile):

y <- cbind(Pop=factor(gsub("_.*", "", x$Sample_ID)),x)
#head(y)
#levels(y$Pop)

#y <- cbind(Pop=factor(gsub("[0-9]","",x$TK)),x)
#y <- cbind(Pop=factor(gsub("[0-9]", "", x$Individual)),x)

#Convert numeric SNP alleles into bases. Seems like there should be a better way to do this, but I (ECA) couldn't figure it out:

f <- y
f[f=="1"] <- "A"  
f[f=="2"] <- "C"
f[f=="3"] <- "G"
f[f=="4"] <- "T"

str(f)
#Now remove loci with ALL missing data (this is probably unnecessary if the infile is already clean):

hh <- f[!sapply(f, function(x) all(is.na(x)))]

# make a list of data frames, one for each population:

pop.idxs <- split(1:nrow(hh), hh$Pop)
g.list<-list(); for(i in seq_along(pop.idxs)) {g.list[[names(pop.idxs)[i]]]=hh[pop.idxs[[i]],] }

# Drop loci that are missing all observations on a population-by-population basis:

g.list.dropped <- list(); for(i in seq_along(g.list)) {
  g.list.dropped[[i]]=g.list[[i]][,!sapply(g.list[[i]],
                                           function(x) all(is.na(x)))]}
names(g.list.dropped) <- names(g.list) # put the names back there

# make Genotypes from those, using a unique convert list for each of them... Sort of messy:

library(genetics)
g.list.genotypes <- lapply(g.list.dropped, function(x) {
  ll<-list();
  for(i in seq(4,ncol(x),by=2)) {
    ll[[i/2-1]] = c(i-1,i)
  }
  makeGenotypes(x, convert=ll)
})

# Finally fix the names of the loci in all of them (From R24370/R24370.1 to R24370):

for(i in seq_along(g.list.genotypes)) {names(g.list.genotypes[[i]])<-gsub("/.*", "", names(g.list.genotypes[[i]])) }

#### Everything to this point is only data manipulation and transformation to meet the Library "genetics" input data format which is a list of 43 populations with respective genotypes called g.list.genotypes #####

# and use this data.frame to get the LD results:

LD.by.pop <- lapply(g.list.genotypes, LD)

##DEP note: This takes a couple of minutes for large datasets and gives a bunch or warning messages indicating that there are monomorphic loci in each of the individual populations. No worries-- ignore them!
#It takes a long time with a RAD data set!

# Save this file in case you want to come back to it later and you don't want to run the LD again:
# save(LD.by.pop, file = "LD.by.pop.Rdata")

# Make list of locus names (which can be reordered for figures later on):

all.locus.names<-names(hh)[-c(1,2)]

# Extract all R^2 matrices and put them in a single list:

Rsquared.by.pop <- lapply(LD.by.pop, "[[", "R^2")

# Copy the upper triangle to the lower triangle of all matrices stored in the list.

Fu <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
} # First, work up a function that copies the upper triangle to the lower triangle to a single list element (a single upper-triangular matrix)

Rsquared.M <- lapply(Rsquared.by.pop, Fu) # Then apply that funtion to all the list elements

# Define a function to expand the matrix for all loci so all pops are comparable:

FLM <- all.locus.names[c(T,F)]

full.locus.matrix <- function(m, aln=FLM) {
  missers <- setdiff(FLM, dimnames(m)[[1]])
  nm <- length(missers)
  if(nm==0) {
    return(m)
  }
  else {
    m1 <- cbind(m,matrix(NA,nrow=nrow(m),ncol=nm))
    m2 <- rbind(m1, matrix(NA, nrow=nm, ncol=ncol(m1)))
    dimnames(m2) <- list(c(dimnames(m)[[1]],missers), c(dimnames(m)[[2]],missers)) 
    return(m2)
  }}

# Make a full.locus.matrix of Rsquared.M to export:
Rsquared.MM <- lapply(Rsquared.M, full.locus.matrix)

#write.csv(Rsquared.MM$MFEel,file='MMtest_MFEel.csv') to extract matrices by population.


## Calculate the mean R-square by locus (by row):
Rsquared.mean <- lapply(Rsquared.MM, function(m) {rowMeans(m, na.rm = TRUE)}) # note the use of na.rm = TRUE to take NA into acount, if this is FALSE then the result is NA.

# sort loci to be in the same order for all populations:
Rsquared.mean.names <- names(Rsquared.mean[[1]])
Mean.by.pop <- do.call(what = cbind, lapply(Rsquared.mean, function(y) y[Rsquared.mean.names]))

write.csv(Mean.by.pop,file='Rsquared_mean.csv')


## Count the number of pairwise R-square grater than 0.95 by locus (by row):
Rsquared.sum95 <- lapply(Rsquared.MM, function(m) {rowSums(m >=0.95, na.rm = TRUE)})
Rsquared.sum95.names <- names(Rsquared.sum95[[1]])
sum95.by.pop <- do.call(what = cbind, lapply(Rsquared.sum95, function(y) y[Rsquared.sum95.names]))

write.csv(sum95.by.pop,file='Rsquared_sum95.csv')






## There are multiple loci for which only a few pairwise R^2 values were obtained, therefore, we want to get a weighted estimate of linkage that takes into account the number of R^2 values for each loci.
## One approach could be to divide the total number of R^2 values by the number of values grater than 0.95 (this value can change to 0.80, 0.50, etc.)

#Count the number of values grater than 0 (non-NAs):
Rsquared.NaN <- lapply(Rsquared.MM, function(m) {rowSums(m >=0, na.rm = TRUE)})
Rsquared.NaN.names <- names(Rsquared.NaN[[1]])
NaN.by.pop <- do.call(what = cbind, lapply(Rsquared.NaN, function(y) y[Rsquared.NaN.names]))

write.csv(NaN.by.pop,file='Rsquared_NaN.csv')

## Divide the number of values (non-NAs) by the number of pairwise R-square grater than 0.95:

Rsquared.wt95 <- mapply('/', Rsquared.sum95, Rsquared.NaN, SIMPLIFY=FALSE)

Rsquared.wt.names <- names(Rsquared.wt95[[1]])
for(i in seq_along(Rsquared.wt.names)) {names(Rsquared.wt.names[[i]])<-gsub("*.", "", names(Rsquared.wt.names[[i]])) }

write.csv(Rsquared.wt95,file='Rsquared_wt95.csv')

# Sort Rsquared.wt95: #### working on this. I've tried a few things but nothing worked ####






#Alicia, figure how to loop this shit! First we need to decide how we want to define the linkage block. Then, it would be great if the end product was a matrix of pops by loci with some metric of LD 
#(mean, or counts exceeding a threshold, or ???). Or, conversely, for each pair of loci we could count the number or proportion of populations in which they are highly linked.	

#Make LD heatmap--I have ditched 'feel.the.heat2' function here and just run LDheatmap directly:
#For unordered heatmap:
#manord <- FLM
#To modify to sort by R^2, etc.
#manord<-FLM[c(1, 17, 24, 29, 31, 33, 46, 51, 52, 56, 59, 61, 71, 72, 77, 79, 81, 86, 88, 11, 20, 41, 43, 49, 55, 62, 66, 67, 74, 75, 3, 9, 14, 16, 18, 19, 44, 58, 68, 80, 83, 87, 4, 40, 65, 85, 78, 6, 28, 37, 69, 70, 5, 45, 76, 63, 21, 38, 36, 10, 25, 60, 2, 27, 23, 82, 53, 64, 12, 15, 35, 54, 32, 42, 57, 73, 7, 22, 34, 30, 47, 84, 26, 50, 48, 13, 8, 39)]
#For RAD data sorted by R^2:
manord <- FLM[c(12, 29, 30, 32, 37, 42, 43, 59, 76, 82, 104, 107, 118, 123, 125, 126, 130, 135, 137, 142, 147, 152, 174, 175, 201, 204, 210, 218, 226, 239, 268, 270, 287, 309, 311, 321, 333, 213, 228, 232, 246, 247, 196, 203, 211, 220, 323, 339, 176, 184, 187, 189, 190, 191, 192, 194, 68, 14, 40, 41, 48, 49, 96, 99, 112, 116, 117, 119, 120, 62, 124, 136, 149, 74, 75, 78, 153, 159, 160, 170, 282, 285, 295, 300, 303, 307, 315, 316, 86, 164, 172, 115, 271, 33, 291, 161, 110, 183, 64, 83, 44, 343, 215, 263, 22, 36, 145, 60, 199, 167, 69, 293, 38, 148, 165, 332, 105, 77, 81, 163, 255, 3, 233, 290, 216, 224, 98, 206, 273, 334, 221, 70, 15, 294, 235, 251, 236, 185, 13, 241, 94, 229, 92, 322, 154, 17, 93, 103, 193, 79, 39, 205, 178, 237, 180, 217, 7, 274, 80, 58, 259, 299, 143, 336, 8, 28, 208, 240, 265, 46, 134, 254, 296, 138, 283, 222, 338, 281, 87, 158, 177, 212, 209, 121, 344, 51, 335, 318, 26, 168, 288, 129, 312, 289, 326, 238, 173, 1, 267, 166, 106, 186, 258, 277, 146, 279, 155, 66, 19, 308, 261, 101, 280, 306, 162, 84, 340, 156, 5, 248, 314, 301, 297, 55, 16, 327, 131, 305, 319, 97, 231, 35, 292, 276, 132, 34, 18, 54, 157, 298, 207, 6, 56, 329, 253, 179, 53, 202, 9, 256, 182, 302, 65, 100, 328, 286, 198, 330, 25, 269, 234, 57, 214, 133, 266, 169, 197, 23, 31, 72, 264, 242, 243, 10, 91, 71, 262, 67, 249, 181, 108, 128, 252, 140, 188, 200, 227, 45, 225, 260, 275, 313, 63, 111, 223, 73, 150, 90, 272, 278, 250, 122, 244, 127, 245, 139, 284, 50, 317, 102, 257, 109, 171, 230, 21, 2, 331, 47, 52, 113, 95, 151, 89, 341, 4, 61, 27, 304, 320, 325, 195, 20, 337, 11, 24, 324, 114, 342, 141, 219, 88, 85, 144, 310)]
LDheatmap(full.locus.matrix(M)[manord,manord], SNP.name=FLM, title="UpqaWSH_RAD", newpage=F)

