
library(Biostrings)
library(tidyverse)
library(ggtree)

###########################################################################

### Visualize Can1/Can2 Phylogenetic Tree
# Tanner Wiegand
setwd("data/")

rm(list=ls())

###########################################################################

# Align sequences using MAFFT (LINSI option)
  # linsi can1_whole_vs_can2_whole_msa_v1.fasta > can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.afa

# Import alignment
x <- readAAStringSet("can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.afa")

# Rename sequences and export FASTA
x.df <- data.frame(name = names(x), label = paste0(word(names(x)), "--", seq(1, length(x))), stringsAsFactors = FALSE)
names(x) <- x.df$label
#writeXStringSet(x, "can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.renamed.afa") # Run this command only once

# Trim columns composed of >70% gaps from MSA using Trimal
  # trimal -in can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.renamed.afa -out can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.renamed.70trm.afa -gt 0.7

# Convert alignment to phylip format
  # seqret -osf phylip -sequence can1_trimmed_vs_can2_whole_msa.70trm.renamed.afa -outseq can1_trimmed_vs_can2_whole_msa.70trm.renamed.phy
  # Above didn't work, instead, output from AliView: can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.renamed.70trm.phy

# Run ProtTest3
  # java -jar ~/Bioinformatics/prottest-3.4.2/prottest-3.4.2.jar -i can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.renamed.70trm.phy -all-distributions -F -AIC -BIC -tc 0.5 > Best_model_for_phyML.txt
    # Best model selected = LG+G+F
      # LG = Generalized model from Le and Gascuel (2008)
      # G = Discrete Gamma model, with default  4 rate categories
      # F = Empirical AA frequencies from the data. AliSim will randomly generate the base frequencies from Uniform distribution.

# Run IQ tree to build phylogenetic tree
  # iqtree -s can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.renamed.70trm.phy -m LG+G+F -bb 1000 -nt 80
    # Output tree = can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.renamed.70trm.phy.tre

# Import tree
tree <- read.tree("can1_whole_vs_can2_whole_msa_v1.CARF-Nuc.renamed.70trm.phy.tre")
tree.df <- as.tibble(tree)
tree.df <- tree.df %>% left_join(x.df, by = "label")

# Store sequences that we want to highlight on the tree
these <- c("WP_011229147.1--1", "WP_012638362.1--208", "WP_013702306.1--207", "WP_053959480.1--206", "637048371--1945", "1XMX--212", "2730024700--2893", "WP_195972101.1--2949", "WP_090010930.1--2730", "WP_025817671.1--3062", "WP_152640288.1--3108")
these.2 <- c("TthCan1", "TsuCan2", "TsuCARD1", "SthCan2", "VC1889", "VC1889-PDB", "AaCan2", "CthCan2", "MbaCan2_a", "MbaCan2_b", "TthCan2")
these.df <- data.frame(label = these, Protein = these.2, stringsAsFactors = FALSE)

# Visualize tree
p <- ggtree(tree, layout = "equal_angle") %<+% these.df
p + geom_tippoint(aes(color = Protein)) + theme(legend.position = "right")



### Visualize tree with originating organism info
# Add genus name to tree
tree.df$organism <- word(tree.df$name, start = -2, end = -1)
y <- grepl("\\[", tree.df$name)
tree.df$organism[y] <- gsub("\\]", "", word(tree.df$name[y], start = 2, end = -1, sep = "\\["))
tree.df$genus <- word(tree.df$organism)

# Import temperature data from Murat
temp <- read_tsv("temperature_data.tsv")
temp$class <- word(temp$lineage_text, start = 4, sep = "\\; ")
temp$genus <- word(temp$organism, sep = "\\_")
temp$genus <- str_to_title(temp$genus)
#temp <- temp %>% group_by(genus) %>% summarize(class = first(class))
temp <- temp %>% group_by(genus) %>% summarize(temp = mean(temperature))
temp <- temp %>% ungroup()

# Add notation
temp$kind <- NA
temp$kind[temp$temp <= 18] <- "Psychrophile"
temp$kind[temp$temp > 18 & temp$temp <= 25] <- "Psycrotroph"
temp$kind[temp$temp > 25 & temp$temp <= 45] <- "Mesophile"
temp$kind[temp$temp > 45 & temp$temp <= 70] <- "Thermophile"
temp$kind[temp$temp > 70] <- "Hyperthermophile"
temp <- temp %>% select(genus, kind)

# Join class info to tree data
tree.df <- tree.df %>% left_join(temp, by = "genus")

# Make list of lists for OTU mapping of genus name
class.list <- list()
#y <- unique(tree.df$class)
y <- unique(tree.df$kind)
y <- y[!is.na(y)]

#for(i in 1:length(y)){ class.list[[paste0(y[i])]] <- tree.df$node[tree.df$class == y[i] & !is.na(tree.df$class)] }
for(i in 1:length(y)){ class.list[[paste0(y[i])]] <- tree.df$node[tree.df$kind == y[i] & !is.na(tree.df$kind)] }
class.list

# Visualize tree
groupOTU(p, class.list, 'Organism') + aes(color = Organism) + theme(legend.position = "right")



# Put Can1 vs Can2 on tree as color
key <- read_tsv("label_can1-can2.txt")
colnames(key)[1] <- "acc"
tree.df$acc <- word(tree.df$label, sep = "--")
tree.df2 <- tree.df %>% left_join(key, by = "acc")

key.list <- list()
y2 <- unique(tree.df2$Protein)
y2 <- y2[!is.na(y2)]

for(i in 1:length(y2)){ key.list[[paste0(y2[i])]] <- tree.df2$node[tree.df2$Protein == y2[i] & !is.na(tree.df2$Protein)] }


# Visualize tree
groupOTU(p, key.list, 'Protein') + aes(color = Protein) + theme(legend.position = "right")






