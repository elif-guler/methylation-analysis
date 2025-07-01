# Load necessary libraries
library(minfi)
library(IlluminaHumanMethylation450kmanifest)

# Step 1: Load raw data
read.csv('Input_Data/SampleSheet_Report_II.csv')

baseDir <- ('Input_data')
targets <- read.metharray.sheet(baseDir)# replace with actual path
targets
RGset <- read.metharray.exp(targets = targets)

# Step 2: Create Red and Green channel matrices
Red <- getRed(RGset)
Green <- getGreen(RGset)

dim(Red)
head(Red)
dim(Green)
head(Green)

# Step 3: Extract signal intensities for your assigned address (40796508)
address <- "40796508"
probe_red <- Red[address, ]
probe_green <- Green[address, ]
probe_red
probe_green

Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID =='40796508','Infinium_Design_Type']
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID =='40796508','Color_Channel']

#Get probe type and color
manifest <- getManifest(RGset)
probeInfo <- getProbeInfo(manifest)
type_info <- probeInfo[probeInfo$AddressA == address | probeInfo$AddressB == address, ]
print(type_info)
head(type_info)

intensity_table <- data.frame(
  Sample = colnames(Red),
  Red = as.numeric(probe_red),
  Green = as.numeric(probe_green)
)

print(intensity_table)

# Step 4: Create MSet.raw and QC plot
MSet.raw <- preprocessRaw(RGset)
qc <- getQC(MSet.raw)
qc
plotQC(qc)


# Calculate detection p-values
detP <- detectionP(RGset)
failed_positions <- detP > 0.05
print(failed_positions)
summary(failed_positions)


# Step 6: Raw beta and M values
beta_df <- data.frame(getBeta(MSet.raw))
M_df <- data.frame(getM(MSet.raw))

pheno <- read.csv('Input_Data/SampleSheet_Report_II.csv')
pheno[,'Group']
betaWT <- beta_df[,pheno[,'Group']=='CTRL']
betaMUT <- beta_df[, pheno[,'Group'] == 'DIS']
MWT = M_df[, pheno[,'Group']== 'CTRL']
MMUT = M_df[, pheno[,'Group']== 'DIS']

mean_of_betaWT = apply(betaWT, 1, mean, na.rm = T)
mean_of_betaMUT = apply(betaMUT, 1, mean, na.rm =T)
mean_of_MWT = apply(MWT, 1, mean, na.rm = T)
mean_of_MMUT = apply(MMUT, 1, mean, na.rm = T)

head(mean_of_MMUT)

d_mean_of_betaWT = density(mean_of_betaWT, na.rm = T)
d_mean_of_betaMUT = density(mean_of_betaMUT, na.rm = T)
d_mean_of_MWT = density(mean_of_MWT, na.rm= T)
d_mean_of_MMUT = density(mean_of_MMUT, na.rm = T)

par(mfrow = c(1,2))
plot(d_mean_of_betaWT, main = 'Density of WT Beta Values', col = 'orange')
lines(d_mean_of_betaMUT, main = 'Density of MUT Beta Values', col= 'green')

## density of beta wt m values plot
plot(d_mean_of_MWT , main= 'Density of WT M Values', col='orange')
lines(d_mean_of_MMUT , main= 'Density of MUT M values', col= 'green')


# Step 7: Normalize using preprocessQuantile
MSet.norm <- preprocessQuantile(RGset)
beta.norm <- getBeta(MSet.norm)

# Step 8: 6-panel plot (density + sd + boxplot, raw and norm)
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dim(dfI)
str(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)


# Now we subset the beta matrix in order to retain only the rows whose name is in the first column of dfI...
beta_I <- beta_df[rownames(beta_df) %in% dfI$IlmnID,]
beta_II <- beta_df[rownames(beta_df) %in% dfII$IlmnID, ]
dim(beta_I)

# For each probe in the mean_of_beta_I and mean_of_beta_II matrices, we calculate the mean of beta values across the 8 samples...
mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)

# ... and then we calculate the density distribution of the 2 vectors of mean values:
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)

#standard deviations
sd_of_beta_I <- apply(beta_I, 1, sd, na.rm =T)
sd_of_beta_II <- apply(beta_II, 1, sd, na.rm= T)
d_sd_of_beta_I <- density(sd_of_beta_I, na.rm= T)
d_sd_of_beta_II <- density(sd_of_beta_II, na.rm = T)

#we normalise the data with our assigned function
preprocessQuantile_results <- preprocessQuantile(RGset)
beta_preprocessQuantile <- getBeta(preprocessQuantile_results)


## We calculate the mean and the standartd deviation for each probe across the 8 samples and calculate the density distributions
beta_preprocessQuantile_I <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfI$IlmnID,]
beta_preprocessQuantile_II <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfII$IlmnID,]

mean_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,mean)
mean_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,mean)
d_mean_of_beta_preprocessQuantile_I <- density(mean_of_beta_preprocessQuantile_I,na.rm=T)
d_mean_of_beta_preprocessQuantile_II <- density(mean_of_beta_preprocessQuantile_II,na.rm=T)

sd_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,sd)
sd_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,sd)
d_sd_of_beta_preprocessQuantile_I <- density(sd_of_beta_preprocessQuantile_I,na.rm=T)
d_sd_of_beta_preprocessQuantile_II <- density(sd_of_beta_preprocessQuantile_II,na.rm=T)

# Finally, we can plot the densities.
pheno$Group
palette(c('orange', 'green'))
par(mfrow = c(2,3))
plot(d_mean_of_beta_I,col="blue", main= 'raw beta')
lines(d_mean_of_beta_II,col="red")


# But in this case it is more meaningful to overlay the two density distributions. To this aim, we first open a plot using the plot() function, and then we add a new line, coloured with a different color, using the lines() function
plot(d_mean_of_beta_I,col="blue" , main= 'raw sd')
lines(d_mean_of_beta_II,col="red")

palette(c("orange", "green"))

beta_df 

# Boxplots

boxplot(beta_df, ylab= 'probes mean', xlab= 'samples' , main='boxplot raw beta', col = as.numeric(as.factor(pheno$Group)), names= pheno$SampleID )
plot(d_mean_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile beta",xlim=c(0,1),ylim=c(0,5))
lines(d_mean_of_beta_preprocessQuantile_II,col="red")

#density plot normalizes SD values
plot(d_sd_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile sd")
lines(d_sd_of_beta_preprocessQuantile_II,col="red")

#box plot of normalized b-values.
boxplot(beta_preprocessQuantile, ylab="probes mean",xlab="probes mean",main="boxplot normalized beta",col = as.numeric(as.factor(pheno$Group)), names= pheno$SampleID)


##Step 8: PCA

pca_result = prcomp(t(beta_preprocessQuantile), scale =T)
plot(pca_result, main= 'pca results')
pca_result$x

plot(pca_result$x[,1], pca_result$x[,2],cex=1,pch=2, col = as.numeric(as.factor(pheno$Group)),xlab="PC1",ylab="PC2", main= 'PCA' ,xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_result$x[,1], pca_result$x[,2],labels= as.character(pheno$SampleID),cex=0.8,pos=1)


plot(pca_result$x[,1], pca_result$x[,2],cex=2,pch=17, col = as.numeric(as.factor(pheno$Group)),xlab="PC1",ylab="PC2", main= 'PCA by group' ,xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_result$x[,1], pca_result$x[,2],labels= as.character(pheno$SampleID),cex=0.8,pos=3, col="black")
legend("bottomright",legend = levels(factor(pheno$Group)),col = 1:nlevels(factor(pheno$Group)),pch = 17)
pheno$Sex
palette(c('pink', 'lightblue'))
plot(pca_result$x[,1], pca_result$x[,2],cex=2,pch=17, col = as.numeric(as.factor(pheno$Sex)),xlab="PC1",ylab="PC2", main= 'PCA by sex' ,xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_result$x[,1], pca_result$x[,2],labels= as.character(pheno$SampleID),cex=0.8,pos=1)

pheno$Sentrix_ID = factor(pheno$Sentrix_ID)
palette(c('black', 'gray', 'gold', 'magenta', 'purple'))
plot(pca_result$x[,1], pca_result$x[,2], main= 'pca by batch', xlab = 'PC1', ylab= 'PC2', cex= 2, pch=17, col= as.numeric(as.factor(pheno$Sentrix_ID)), xlim=c(-1000,1000),ylim=c(-1000,1000) )
text(pca_result$x[,1], pca_result$x[,2],labels= as.character(pheno$SampleID),cex=0.8,pos=1)

## Step 9: t-test

t_test <- t.test(beta_preprocessQuantile[1,] ~ pheno$Group)
t_test

str(t_test)

t_test$p.value

my_t_test <- function(x) {
  t_test <- t.test(x~ pheno$Group)
  return(t_test$p.value)
} 


# We apply the function to the rows of beta_preprocessQuantile
pValues_ttest <- apply(beta_preprocessQuantile,1, my_t_test)
pValues_ttest
length(pValues_ttest)

final_ttest <- data.frame(beta_preprocessQuantile, pValues_ttest)
head(final_ttest)

final_ttest$pValues_ttest

final_ttest_0.05 <- final_ttest[final_ttest$pValues_ttest<=0.05,]
head(final_ttest_0.05)
dim(final_ttest_0.05)



##Apply Multiple Correction Method

corrected_pValues_BH <- p.adjust(final_ttest$pValues_ttest,"BH")
corrected_pValues_Bonf <- p.adjust(final_ttest$pValues_ttest,"bonferroni")
final_ttest_corrected <- data.frame(final_ttest, corrected_pValues_BH, corrected_pValues_Bonf)
head(final_ttest_corrected)


# We can visualize the distributions of the p-values and of the corrected p-values by boxplots:
colnames(final_ttest_corrected)
par(mfrow = c(1,1))
boxplot(final_ttest_corrected[,9:11], main= 'distribution of p-values')

# How many probes survive the multiple test correction?
dim(final_ttest_corrected[final_ttest_corrected$pValues_ttest<=0.05,])
dim(final_ttest_corrected[final_ttest_corrected$corrected_pValues_BH<=0.05,])
dim(final_ttest_corrected[final_ttest_corrected$corrected_pValues_Bonf<=0.05,])


## CHUNK 2: Volcano plots
# First of all, we have to calculate the difference between the average of group A values and the average of group B values. To this aim, we will first create two matrixes containing the beta-values of group A and group B samples, and then we will calculate the mean within each group for each row.
pheno$Group

beta_corrected <- final_ttest_corrected[,1:8]

beta_corrected_group_CTRL <- beta_corrected[,pheno$Group=="CTRL"]
mean_beta_corrected_group_CTRL <- apply(beta_corrected_group_CTRL,1,mean)
beta_corrected_group_DIS <- beta_corrected[,pheno$Group=="DIS"]
mean_beta_corrected_group_DIS <- apply(beta_corrected_group_DIS,1,mean)

head(beta_corrected_group_CTRL)

# Now we can calculate the difference between average values:
delta_corrected <- mean_beta_corrected_group_DIS-mean_beta_corrected_group_CTRL
head(delta_corrected)

# Now we create a dataframe with two columns, one containing the delta values and the other with the -log10 of p-values
toVolcPlot <- data.frame(delta_corrected, -log10(final_ttest_corrected$pValues_ttest))
head(toVolcPlot)
plot(toVolcPlot[,1], toVolcPlot[,2], main= 'volcano plot') 


plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5, main= 'volcano plot with treshold')
-log10(0.01)
abline(h=-log10(0.01),col="red")

# Now we want to highlight the probes, that have a nominal pValue<0.01 and a delta > 0.1)
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5, main= 'pvalue<0.01 and delta > 0.1')
toHighlight <- toVolcPlot[toVolcPlot[,1]>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="yellow")

# If we want to highlight the points with an absolute delta > 0.01
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5, main= 'delta > 0.01')
toHighlight <- toVolcPlot[abs(toVolcPlot[,1])>0.1 & toVolcPlot[,2]>(-log10(0.01)),]
head(toHighlight)
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="red")

install.packages("qqman")
library(qqman)

load("final_ttest_reduced.RData")
# First we have to annotate our dataframe, that is add genome annotation information for each cpg probe. We will use the Illumina450Manifest_clean object that we previously created:
load('~/Illumina450Manifest_clean.RData')

# We will use the merge() function to merge the final_ttest_corrected with the Illumina450Manifest_clean object
?merge
# The merge function performs the merging by using a column which is common to two dataframes and which has the same name in the two dataframes
head(Illumina450Manifest_clean)
head(final_ttest_reduced)
# We want to merge on the basis of the CpG probes, but unfortunately in the final_ttest_corrected object the CpG probes are stored in the rownames, not in a column. We can overcome this issue as follows:
final_ttest_reduced <- data.frame(rownames(final_ttest_reduced),final_ttest_reduced)
head(final_ttest_reduced)
colnames(final_ttest_reduced)
colnames(final_ttest_reduced)[1] <- "IlmnID"
colnames(final_ttest_reduced)

final_ttest_reduced_annotated <- merge(final_ttest_reduced, Illumina450Manifest_clean,by="IlmnID")
dim(final_ttest_reduced)
dim(Illumina450Manifest_clean)
dim(final_ttest_reduced_annotated)
str(final_ttest_reduced_annotated)


input_Manhattan <- final_ttest_reduced_annotated[colnames(final_ttest_reduced_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_ttest")]
dim(input_Manhattan)
head(input_Manhattan)
str(input_Manhattan$CHR)
levels(input_Manhattan$CHR)
# It is better to reorder the levels of the CHR
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )
levels(input_Manhattan$CHR)

input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)

# and finally we can produce our Manhattan plot
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest" , main= 'manhattan plot')
# what is the blue line?
-log10(0.00001)
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest",annotatePval = 0.00001,col=rainbow(24), main= 'manhattan plot' )

install.packages('gplots')
library(gplots)

input_heatmap=as.matrix(final_ttest[1:100,1:8])

# Finally, we create a bar of colors for Group A or Group B membership
pheno$Group
colorbar <- c("green","green","orange","orange","green","green","orange","orange")

# In the following lines we will compare the results of hierchical clustering using different methods.

# Complete (default options)
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage")

# You can see that hierarchical clustering divides well group A and group B samples; in addition, you see that some probes are hypermethylated in Group A compared to group B, others are hypomethylated.

# Single
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Single linkage")

# Average
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Average linkage")

col2=colorRampPalette(c("green","black","red"))(100)
heatmap.2(input_heatmap,col=col2,Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F)
