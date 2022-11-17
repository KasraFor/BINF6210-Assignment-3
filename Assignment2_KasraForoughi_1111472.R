#BINF6210 Assignment 2---- 

#Name: Kasra Foroughi
#Student ID: 1111472

#Collaborator name: Omar Amin
#Student ID: 1251813

#November 17, 20229


#Part 1: Loading Packages----
#install.packages('tidyverse')
library(tidyverse)
library(ggplot2)
library(cluster)
library(grid)

#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install(c("Biostrings", "DECIPHER", "muscle"))
library(Biostrings)
library(muscle)
library(DECIPHER)

#install.packages("rentrez")
library(rentrez)
#install.packages("seqinr")
library(seqinr)

#install.packages("stringi")
library(stringi)
#install.packages("ape")
library(ape)
#install.packages("RSQLite") #a required package for DECIPHER
library(RSQLite)

library(factoextra)

#install.packages('dendextend')
library(dendextend) 


#Part 2: Load Data sets----
#Step 1 Load Nuccore database of Rattus Genus for Genes COI & CYTB for analysis

#Search for Rattus Species with Gene CYTB from NUCCORE Database--> 2085 HITS
#Checking to see if CYTB vs CYT B spelling was differnt. CYT B had only 6 hits, so CYTB will be the one that is used 
#Rattus_CYTB_Raw2 <- entrez_search(db = 'nuccore', term = "Rattus[ORGN] AND CYT B[GENE]", retmax = 2085, use_history = TRUE)


#Use Entrez_search to extract data from nuccore database with Rattus as the organism and CYTB Gene
#Need to use web history in this instance as the amount of HITS is too large
Rattus_CYTB_Raw <- entrez_search(db = 'nuccore', term = "Rattus[ORGN] AND CYTB[GENE]", retmax = 2085, use_history = TRUE)
#Fetch the data and save it as a Fasta to be used by DNAStringSet 
Rattus_CYTB_Raw_fetch <- entrez_fetch(db = 'nuccore', web_history = Rattus_CYTB_Raw$web_history, rettype = 'fasta')
#Save the data 
write(Rattus_CYTB_Raw_fetch, "Rattus_CYTB_Raw_fetch.fasta", sep = "\n")
Rattus_CYTB_Raw_stringSet <- readDNAStringSet("Rattus_CYTB_Raw_fetch.fasta")

# Create a DataFrame of the Rattus CYTB Genes
Rattus_DF_CYTB_Raw <- data.frame(CytB_Title = names(Rattus_CYTB_Raw_stringSet), CYTB_Sequence = paste(Rattus_CYTB_Raw_stringSet))
#Clean up Datasets
# Make a new column, which extracts the Species Name from the Title 
Rattus_DF_CYTB_Raw$Species_Name <- word(Rattus_DF_CYTB_Raw$CytB_Title, 2L, 3L)
# Rearrange the columns.
Rattus_DF_CYTB_Raw <- Rattus_DF_CYTB_Raw[, c("CytB_Title", "Species_Name", "CYTB_Sequence")]

#The process above was done again for Rattus COI gene.
#Search for Rattus with GENE COI --> 1184 HITS
Rattus_COI_Raw <- entrez_search(db = 'nuccore', term = "Rattus[ORGN] AND COI[GENE]", retmax = 1184, use_history = TRUE)
Rattus_COI_Raw_fetch <- entrez_fetch(db = 'nuccore', web_history = Rattus_COI_Raw$web_history, rettype = 'fasta')
write(Rattus_COI_Raw_fetch, "Rattus_COI_Raw_fetch.fasta", sep = "\n")
Rattus_COI_Raw_stringSet <- readDNAStringSet("Rattus_COI_Raw_fetch.fasta")
# Create a Dataframe of Rattus COI Genes
Rattus_DF_COI_Raw <- data.frame(CytB_Title = names(Rattus_COI_Raw_stringSet), COI_Sequence = paste(Rattus_COI_Raw_stringSet))
#Clean up the Database
Rattus_DF_COI_Raw$Species_Name <- word(Rattus_DF_COI_Raw$CytB_Title, 2L, 3L)
# Rearrange the columns.
Rattus_DF_COI_Raw <- Rattus_DF_COI_Raw[, c("CytB_Title", "Species_Name", "COI_Sequence")]
#View(Rattus_DF_COI_Raw)

#Part 3: Data Exploration and Cleaning----

#For clustering, wanted to make sure there is an intersect between both gene pools
#Determine the unique Species name in CYTB and COI gene information
unique_Rattus_COI_Raw <- unique(Rattus_DF_COI_Raw$Species_Name)
unique_Rattus_CYTB_Raw <- unique(Rattus_DF_CYTB_Raw$Species_Name)
#Determine which Species names are common amongst both datasets
Rattus_Species <- intersect(unique_Rattus_COI_Raw, unique_Rattus_CYTB_Raw)
#remove Rattus Sp. and Rattus novaeguineae
#Rattus Sp. is an unkown species. Rattus novaeguineae was removed due to issues with alignment due to lack of data and huge discrepency between DNA length sequences
Rattus_Species <- Rattus_Species[!Rattus_Species %in% c("Rattus sp.", "Rattus novaeguineae")]
#Rattus_Species

#Summary of the COI and CYTB sequences
#This is to determine if there are outliers
summary(nchar(Rattus_DF_COI_Raw$COI_Sequence))
summary(nchar(Rattus_DF_CYTB_Raw$CYTB_Sequence))
#We can see that there are sequences which exceed the average bp length of both genes and ones which are too little

#Quick Histogram to determine if outliers are present
hist(nchar(Rattus_DF_COI_Raw$COI_Sequence))
hist(nchar(Rattus_DF_CYTB_Raw$CYTB_Sequence))
#Histograms show that there are outliers which need to be handled

#In both sequences, there seems to be outliers present
#For COI Gene, the common BP length is 650 - so I have chosen to only use data points that are within +/- 200 bp length of 650 bp
#Remove outliers in COI Gene sequences
#Only use Rattus Species which are present in both genes
Rattus_DF_COI <- Rattus_DF_COI_Raw %>% 
  filter(Species_Name %in% Rattus_Species) %>% 
  filter(str_length(COI_Sequence) >= 450 & str_length(COI_Sequence) <850) %>% 
  filter(!is.na(COI_Sequence)) 
#Determine how many data points were removed
#123 Data points were removed post filter of COI BP sequences
nrow(Rattus_DF_COI_Raw) - nrow(Rattus_DF_COI)
#new summary and histogram to check if outliers are removed
summary(nchar(Rattus_DF_COI$COI_Sequence))

#EDIT 1: Just edit the boxplot
#Create a box plot to visualize the spread of the COI Sequence lengths
COI_His <- nchar(Rattus_DF_COI$COI_Sequence)

#####################my edit######################
COI_df <- as.data.frame(COI_His)

findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

df <- COI_df %>%
  mutate(outlier = ifelse(findoutlier(COI_His), COI_His, NA))



ggplot(data.frame(df), aes(x="", y = COI_His)) +
  geom_boxplot() +
  geom_text(aes(label=outlier), na.rm=TRUE, hjust=-.5)

#################################################
#boxplot(COI_His, horizontal =TRUE,
#        main = 'Distribution of COI Sequence Lengths',
#        xlab = "Length of COI Sequence (Bp)")

#For Cyt B gene, the common BP length is 1140, however, there is higher variance within the range, thus I have chose to use
#data points that are within 500+/- BP 
Rattus_DF_CYTB <- Rattus_DF_CYTB_Raw %>% 
  filter(Species_Name %in% Rattus_Species) %>% 
  filter(str_length(CYTB_Sequence) >= 600 & str_length(CYTB_Sequence) <1400) %>% 
  filter(!is.na(CYTB_Sequence))

#Determine how many data points were removed
#423 Data points were removed post filter of COI BP sequences
nrow(Rattus_DF_CYTB_Raw) - nrow(Rattus_DF_CYTB)
#summary to check 
summary(nchar(Rattus_DF_CYTB$CYTB_Sequence))

#Create a box plot to visualize the spread of the CYTB Sequence lengths
CYTB_His <- nchar(Rattus_DF_CYTB$CYTB_Sequence)
boxplot(CYTB_His, horizontal =TRUE,
        main = 'Distribution of CYTB Sequence Lengths',
        xlab = "Length of CYTB Sequence (Bp)")


#Part 4: Data Analysis: Alignment----
#This section will involve alignment of the respective genes 
#The basis is to first align each sequence grouped by species name, create a consensus sequence for each species name before an alignment

############################## Omar's EDIT (1 & 2) ##############################
#Edit 1
#In part 4 of Kasra's code, he performs a very smart method of separating the large df into a list of dfs, each df has the sequences of a specific  species in the group, and then he does a function that generates a consensus sequence after alignment in each species, after this function is done, he performs two lapply functions on the list of dfs, he changes the DNAStringSet character and then unlists all the dfs into a single df (containing al the species' consensus sequence). then, he converts the rownames into a single column and he changes the column names as "Species_Name" and "Consensus_COI_Sequence". 
#I think it might be more efficient to do all these steps inside one function, the input for the function will be the large df without splitting it into a list of dfs, thus making it more efficient and easier to use with fewer lines and it is reproducible with any data frame that has many species with many sequences if we want to use this function for future projects.
#Also, he created two functions for both COI & CYTB, because they have a different column name for the sequence (COI_Sequence, CYTB_Sequence), I think it will be more efficient to use one function that takes the name of the sequence column of any input df and deals with it throughout the function (as long as the column containing the sequences is in the third column).

#Edit 2
#Inside the function, there was a filtering process in which it removes the N if it passes 1 percent of sequence, I added a very basic regular expression inside (pattern = [^ATCG]) that looks for anything other than the ATGC nucleotides, while also keeping the 1 percent threshold.

omar_cons_func <- function(df){
  missing.data <- 0.01
  length.var <- 50
  seq_col_name <- paste0("Consensus_", colnames(df)[3])
  df_final <- lapply(split(df, f = df$Species_Name), function(df1){
    df1 <-df1 %>%
      mutate(nucleotides2 = str_remove_all(df1[,3], "^N+|N+$|-")) %>%
      filter(str_count(nucleotides2, "[^ATCG]") <= (missing.data * str_count(df1[,3]))) %>%
      filter(str_count(nucleotides2) >= median(str_count(nucleotides2)) - length.var & str_count(nucleotides2) <= median(str_count(nucleotides2)) + length.var)
    df1$nucleotides2 <- DNAStringSet(df1$nucleotides2)
    names(df1$nucleotides2) <- df1$Species_Name
    testCOI_alignment <- DNAStringSet(muscle::muscle(df1$nucleotides2))
    cons <- as.character(ConsensusSequence(testCOI_alignment))})
  df_final <- as.data.frame(unlist(df_final)) %>% 
    tibble::rownames_to_column("Species_Name") %>%
    `colnames<-`(c("Species_Name", seq_col_name))
  return(df_final)}

#And just like that, we use the function on the two large dfs 
COI_seq <- omar_cons_func(Rattus_DF_COI)
CYTB_seq <- omar_cons_func(Rattus_DF_CYTB)

####################################End of Edit 1 & 2 #####################################

#Alignment of Consensus Sequences of COI
#First create another column known as nucleotides 2 and convert it to a DNAStringSet for Bioconductor alignment
COI_All <- COI_seq %>% mutate(nucleotides2 = Consensus_COI_Sequence)
COI_All$nucleotides2 <- DNAStringSet(COI_All$nucleotides2)
names(COI_All$nucleotides2) <- COI_All$Species_Name
#Align all consensus sequence based on Muscle Alignment
COI_alignment <- DNAStringSet(muscle::muscle(COI_All$nucleotides2))

#Alignment of Consensus Sequences of CYTB
CYTB_All <- CYTB_seq %>% mutate(nucleotides2 = Consensus_CYTB_Sequence)
CYTB_All$nucleotides2 <- DNAStringSet(CYTB_All$nucleotides2)
names(CYTB_All$nucleotides2) <- CYTB_All$Species_Name
CYTB_alignment <- DNAStringSet(muscle::muscle(CYTB_All$nucleotides2))
#BrowseSeqs(COI_alignment)


#Part 5: Clustering & Silhouette---- 
#Before clustering, we need to create a distance Matrix for hierarchical clustering
#JC model was chosen because for equal base frequencies, all substitutions are likely 
chosen.model <- "JC69"

#COI distance matrix
#To create a distance matrix, need to convert the alignment into a DNAbin class
test_BIN_COI <- as.DNAbin(COI_alignment)
#Perform distance matrix for clustreing based on JC69 model
distanceMatrix_COI <- dist.dna(test_BIN_COI, model = chosen.model, as.matrix = TRUE, pairwise.deletion = TRUE)
distanceMatrix_COI=as.dist(distanceMatrix_COI)

#CYTB Distance Matrix
test_BIN_CYTB <- as.DNAbin(CYTB_alignment)
distanceMatrix_CYTB <- dist.dna(test_BIN_CYTB, model = chosen.model, as.matrix = TRUE, pairwise.deletion = TRUE)
distanceMatrix_CYTB = as.dist(distanceMatrix_CYTB)

#Source: https://bradleyboehmke.github.io/HOML/hierarchical.html
#In this section, we first figure out which hierarchical clustering method is best suited for the genes. 
#The 4 different hierarchical clustering method
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

#This function uses agnes() which is similar to hclust, but also returns a coefficient (AC) that measures the strength of the respective clustering methods
#COI AC - determine which clustering method is highest for either 
ac_COI <- function(x) {
  agnes(distanceMatrix_COI, method = x)$ac
}
ac_COI_values <- purrr::map_dbl(m, ac_COI) 

#CYTB AC
ac_CYTB <- function(x) {
  agnes(distanceMatrix_CYTB, method = x)$ac
}
ac_CYTB_values<- purrr::map_dbl(m, ac_CYTB) 

which.max(ac_COI_values)
which.max(ac_CYTB_values)
#Ward Clustering method has the highest AC value for both COI and CYTB genes, thus the ward method will be used for clustering 

#Clustering COI
#Using a Ward clustering method 
clustering.method <- "ward.D2"

#Create COI cluster 
COI_Cluster <-hclust(distanceMatrix_COI,method=clustering.method)
#Plot the COI cluster 
plot(COI_Cluster, main = "Rattus Species COI Cluster", xlab = "Rattus Species")
#Clustering CYTB
CYTB_Cluster <-hclust(distanceMatrix_CYTB,method=clustering.method)
#Plot the CYTB Cluster
plot(CYTB_Cluster, main = "Rattus Species CYTB Cluster", xlab = "Rattus Species")

#Source: https://rdrr.io/cran/factoextra/man/fviz_silhouette.html 
#Create a Silhoutte index of both COI and CYTB to determine which is a stronger clustering method 
hc_cut_COI <- hcut(distanceMatrix_COI, k = 3, hc_method = "ward.D2")
# Visualize silhouhette information
fviz_silhouette(hc_cut_COI, main = "COI Silhouette")

hc_cut_CYTB <- hcut(distanceMatrix_CYTB, k = 3, hc_method = "ward.D2")
# Visualize silhouhette information
fviz_silhouette(hc_cut_CYTB, main = "CYTB Silhouette")

#######################################Omar Edit 3########################################
#In Kasra's code above, he used the silhouette index to determine which clustering method of CYTB and COI is the best based on the silhouette width
#In my edit, I wanted to compare the two clusters facing each other using a tanglegram, this method aligns the two clusters together while untangling the dendrograms (while maintaining the structure of the dendrograms accurately). I used the step2side method for untangling as it showed the most efficient untangling. The results were very promising and agrees with Kasra's discussion and conclusion. The figure shows only 3 terminal nodes aligning between both genes (Rattus tiomanicus - Rattus rattus), (Rattus brunneusculus, Rattus andamanensis) and (Rattus norvegicus - Rattus nitidus), and no larger clades aligned between both genes.
COI_dend <- as.dendrogram (COI_Cluster)
CYTB_dend <- as.dendrogram(CYTB_Cluster)

# Create a list to hold dendrograms
cluster_dend <- dendlist(COI_dend, CYTB_dend) %>%
  untangle() %>%
  set("rank_branches") %>%
  tanglegram(margin_inner=5)

pdf("trial.pdf")
dendlist(COI_dend, CYTB_dend) %>%
  untangle(method="step2side") %>%
  set("rank_branches") %>%
  tanglegram(margin_inner=9)
grid.text("COI vs CYTB Clusters", x = (0.5), y = (0.95))

dev.off()
#######################################Edit End#####################################
