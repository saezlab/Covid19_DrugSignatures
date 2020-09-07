CARNIVAL on A549 and CALU-3 cell lines after SARS-CoV-2 infection
================
Alberto Valdeolivas: <alberto.valdeolivas@bioquant.uni-heidelberg.de>;
Date:
27/08/2020

### License Info

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

Please check <http://www.gnu.org/licenses/>.

## Introduction

The present script takes pathway and tf activity on A549 and CALU-3 cell
lines after SARS-CoV-2 infection and runs CARNIVAL

### Getting Ready

We first load the required libraries and we define a function to export
CARNIVAL results to Cytoscape and another one to extract CARNIVAL
results.

``` r
library(tidyverse)
library(OmnipathR)
library(CARNIVAL)
library(gprofiler2)
library(plotly)

OutputCyto <- function(CarnivalResults, outputFile) {
    CarnivalNetwork <- 
        as.data.frame(CarnivalResults$weightedSIF, stringsAsFactors = FALSE) %>%
        dplyr::mutate(Sign = as.numeric(Sign), Weight = as.numeric(Weight)) %>% 
        dplyr::mutate(Weight = Sign * Weight) %>%
        dplyr::select(Node1, Weight, Node2)
        
    CarnivalNetworkNodes <- 
        unique(c(CarnivalNetwork$Node1,CarnivalNetwork$Node2))
    
    CarnivalAttributes <- CarnivalResults$nodesAttributes %>% 
        as.data.frame() %>%
        dplyr::filter(Node %in% CarnivalNetworkNodes) %>%
        dplyr::mutate(NodeType = as.character(NodeType)) %>%
        dplyr::mutate(NodeType=if_else(NodeType =="", "I", NodeType))
            
    nameOutputNetwork <- paste0(outputFile, "Network.sif")
    nameOutputAttributes <-  paste0(outputFile, "Attributes.txt")    
    
    write.table(CarnivalNetwork, file = nameOutputNetwork,
        quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
    
    write.table(CarnivalAttributes, file = nameOutputAttributes,
        quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

extractCARNIVALnodes <- function(CarnivalResults){

    CarnivalNetwork <- 
        as.data.frame(CarnivalResults$weightedSIF, stringsAsFactors = FALSE)
    
    colnames(CarnivalNetwork) <- c("source", "sign", "target", "Weight")

    ## We define the set of nodes interesting for our condition
    sucesses <- unique(c(gsub("_.*","",CarnivalNetwork$source), 
        gsub("_.*","",CarnivalNetwork$target)))

    CarnivalAttributes <- as.data.frame(CarnivalResults$nodesAttributes, 
        stringsAsFactors = FALSE)

    ## We define the background as all the genes in our prior knowledge network.
    bg <- unique(gsub("_.*","",CarnivalAttributes$Node))     
    
    return(list(sucesses = sucesses, bg= bg))
}
```

### Generating the prior knowledge network from Omnipath

``` r
### Prior knowledge network from Omnipath. 
ia_omnipath <- import_omnipath_interactions() %>% as_tibble()
ia_pwextra <- import_pathwayextra_interactions() %>% as_tibble()
ia_kinaseextra <- import_kinaseextra_interactions() %>% as_tibble()

## We bind the datasets
interactions <- as_tibble(
  bind_rows(
    ia_omnipath %>% mutate(type = 'ppi'),
    ia_pwextra %>% mutate(type = 'ppi'),
    ia_kinaseextra %>% mutate(type = 'ppi')))

signed_directed_interactions <- 
  dplyr::filter(interactions, consensus_direction==1) %>%
  filter(consensus_stimulation == 1 | consensus_inhibition == 1) %>% 
  dplyr::mutate(sign = if_else(consensus_stimulation==1,1,-1))  %>%
  dplyr::select(source_genesymbol, sign,  target_genesymbol) %>%
  dplyr::rename(source ="source_genesymbol", target ="target_genesymbol")

carnival_pkn <- signed_directed_interactions %>%
  dplyr::distinct(source, target, .keep_all = TRUE)
saveRDS(carnival_pkn, file="Intermediate_FIles/carnival_pkn.rds")
```

``` r
carnival_pkn<- readRDS("Intermediate_FIles/carnival_pkn.rds")
all_source_nodes <- unique(carnival_pkn$source)
all_target_nodes <- unique(carnival_pkn$target)
all_nodes_network <- unique(c(all_source_nodes,all_target_nodes))
```

### Reading TF and pathway activity and perturbed nodes.

``` r
### Dorothea + Progeny Results 
dorothea_results <- 
  read_csv(file = "CARNIVAL_input/virus_signatures_dorothea.csv") %>% 
  dplyr::rename(Tf = "X1")

progeny_results <- 
  read_csv(file = "CARNIVAL_input/virus_signatures_progeny.csv") %>% 
  dplyr::rename(Pathway = "X1")

### Perturbation nodes
perturbations_nodes <- 
  read_csv(file = "CARNIVAL_input/virus_signatures_pert.csv") %>% 
  dplyr::select(Gene, Mechanism)
```

## Results

### Early phase: RIG-I-Like receptors as perturbation

For the early response to SARS-CoV-2 infection, we set the RIG-I-Like
receptors as perturbed.

``` r
RIG_perturbation <- perturbations_nodes %>%
  dplyr::filter(Mechanism == "RIG") %>% 
  dplyr::mutate(sign = 1) %>% 
  dplyr::select(Gene, sign) %>% 
  tibble::column_to_rownames(var = "Gene") %>%
  t() %>% as.data.frame()
```

We also check if these nodes are present in our prior knowledge network
as source of interactions (otherwise it does not make sense to include
its perturbation).

``` r
colnames(RIG_perturbation) 
```

    ## [1] "DDX58" "IFIH1" "DHX58"

``` r
colnames(RIG_perturbation) %in% all_source_nodes
```

    ## [1] TRUE TRUE TRUE

#### GSE147507.S5\_A549\_SARS.CoV.2

We select this condition and prepared the inputs to run CARNIVAL

``` r
dorothea_results_GSE147507_A549 <- dorothea_results %>% 
  dplyr::select(Tf, GSE147507.S5_A549_SARS.CoV.2) %>% 
  dplyr::filter(Tf %in% all_nodes_network) %>%
  dplyr::top_n(25, wt = abs(GSE147507.S5_A549_SARS.CoV.2)) %>%
  dplyr::arrange(GSE147507.S5_A549_SARS.CoV.2) %>%
  tibble::column_to_rownames(var = "Tf")  %>%
  t() %>% as.data.frame()

progeny_results_GSE147507_A549 <- progeny_results %>% 
  dplyr::select(Pathway,GSE147507.S5_A549_SARS.CoV.2) %>%
  tibble::column_to_rownames(var = "Pathway")  %>%
  t()
```

And we finally run CARNIVAL:

``` r
carnival_results_GSE147507_A549 <-runCARNIVAL(
    solverPath="/home/alvaldeolias/Downloads/cplex",
    netObj=carnival_pkn,
    measObj=dorothea_results_GSE147507_A549,
    inputObj = RIG_perturbation,
    # dir_name="Carnival_Results",
    weightObj=progeny_results_GSE147507_A549,
    # nodeID = 'gene',
    timelimit = 7200,
    solver = "cplex")
saveRDS(carnival_results_GSE147507_A549, 
    file = "CARNIVAL_results/viral_signatures/carnival_results_GSE147507_A549.rds")
OutputCyto(carnival_results_GSE147507_A549, 
    outputFile="CARNIVAL_results/viral_signatures/carnival_results_GSE147507_A549")
```

We now run the enrichment analysis:

``` r
carnival_results_GSE147507_A549 <- 
    readRDS("CARNIVAL_results/viral_signatures/carnival_results_GSE147507_A549.rds")
nodes_carnival_GSE147507_A549 <- 
  extractCARNIVALnodes(carnival_results_GSE147507_A549)

enrichment_results_GSE147507_A549 <- 
  gost(nodes_carnival_GSE147507_A549$sucesses, user_threshold = 0.01, 
  correction_method = "gSCS", custom_bg = nodes_carnival_GSE147507_A549$bg,
  sources = c("KEGG","REAC","WP","GO:BP","MIRNA")) 
```

    ## Detected custom background input, domain scope is set to 'custom'

``` r
gostplot(enrichment_results_GSE147507_A549, capped = FALSE, interactive = TRUE)
```

<img src="Carnival_virus_signatures_files/figure-gfm/plot_GSE147507_A549-1.png" width="800" height="700" />

#### GSE147507.S7\_Calu3\_SARS.CoV.2

We select this condition and prepared the inputs to run CARNIVAL

``` r
dorothea_results_GSE147507_CALU3 <- dorothea_results %>% 
  dplyr::select(Tf, GSE147507.S7_Calu3_SARS.CoV.2) %>% 
  dplyr::filter(Tf %in% all_nodes_network) %>%
  dplyr::top_n(25, wt = abs(GSE147507.S7_Calu3_SARS.CoV.2)) %>%
  dplyr::arrange(GSE147507.S7_Calu3_SARS.CoV.2) %>%
  tibble::column_to_rownames(var = "Tf")  %>%
  t() %>% as.data.frame()

progeny_results_GSE147507_CALU3 <- progeny_results %>% 
  dplyr::select(Pathway,GSE147507.S7_Calu3_SARS.CoV.2) %>%
  tibble::column_to_rownames(var = "Pathway")  %>%
  t()
```

And we finally run CARNIVAL:

``` r
carnival_results_GSE147507_CALU3 <-runCARNIVAL(
    solverPath="/home/alvaldeolias/Downloads/cplex",
    netObj=carnival_pkn,
    measObj=dorothea_results_GSE147507_CALU3,
    inputObj = RIG_perturbation,
    # dir_name="Carnival_Results",
    weightObj=progeny_results_GSE147507_CALU3,
    # nodeID = 'gene',
    timelimit = 7200,
    solver = "cplex")
saveRDS(carnival_results_GSE147507_CALU3, 
    file = "CARNIVAL_results/viral_signatures/carnival_results_GSE147507_CALU3.rds")
OutputCyto(carnival_results_GSE147507_CALU3, 
    outputFile="CARNIVAL_results/viral_signatures/carnival_results_GSE147507_CALU3")
```

We now run the enrichment analysis:

``` r
carnival_results_GSE147507_CALU3 <- 
    readRDS("CARNIVAL_results/viral_signatures/carnival_results_GSE147507_CALU3.rds")
nodes_carnival_GSE147507_CALU3 <- 
  extractCARNIVALnodes(carnival_results_GSE147507_CALU3)

enrichment_results_GSE147507_CALU3 <- 
  gost(nodes_carnival_GSE147507_CALU3$sucesses, user_threshold = 0.01, 
  correction_method = "gSCS", custom_bg = nodes_carnival_GSE147507_CALU3$bg,
  sources = c("KEGG","REAC","WP","GO:BP","MIRNA")) 
```

    ## Detected custom background input, domain scope is set to 'custom'

``` r
gostplot(enrichment_results_GSE147507_CALU3, capped = FALSE, interactive = TRUE)
```

<img src="Carnival_virus_signatures_files/figure-gfm/plot_GSE147507_CALU3-1.png" width="800" height="700" />

#### GSE148729\_Calu3\_SARS.CoV.2\_24H

We select this condition and prepared the inputs to run CARNIVAL

``` r
dorothea_results_GSE148729_CALU3 <- dorothea_results %>% 
  dplyr::select(Tf, GSE148729_Calu3_SARS.CoV.2_24H) %>% 
  dplyr::top_n(25, wt = abs(GSE148729_Calu3_SARS.CoV.2_24H)) %>%
  dplyr::arrange(GSE148729_Calu3_SARS.CoV.2_24H) %>%
  dplyr::filter(Tf %in% all_nodes_network) %>%
  tibble::column_to_rownames(var = "Tf")  %>%
  t() %>% as.data.frame()

progeny_results_GSE148729_CALU3 <- progeny_results %>% 
  dplyr::select(Pathway,GSE148729_Calu3_SARS.CoV.2_24H) %>%
  tibble::column_to_rownames(var = "Pathway")  %>%
  t()
```

And we finally run CARNIVAL:

``` r
carnival_results_GSE148729_CALU3 <-runCARNIVAL(
    solverPath="/home/alvaldeolias/Downloads/cplex",
    netObj=carnival_pkn,
    measObj=dorothea_results_GSE148729_CALU3,
    inputObj = RIG_perturbation,
    # dir_name="Carnival_Results",
    weightObj=progeny_results_GSE148729_CALU3,
    # nodeID = 'gene',
    timelimit = 7200,
    solver = "cplex")
saveRDS(carnival_results_GSE148729_CALU3, 
    file = "CARNIVAL_results/viral_signatures/carnival_results_GSE148729_CALU3.rds")
OutputCyto(carnival_results_GSE148729_CALU3, 
    outputFile="CARNIVAL_results/viral_signatures/carnival_results_GSE148729_CALU3")
```

We now run the enrichment analysis:

``` r
carnival_results_GSE148729_CALU3 <- 
    readRDS("CARNIVAL_results/viral_signatures/carnival_results_GSE148729_CALU3.rds")
nodes_carnival_GSE148729_CALU3 <- 
  extractCARNIVALnodes(carnival_results_GSE148729_CALU3)

enrichment_results_GSE148729_CALU3 <- 
  gost(nodes_carnival_GSE148729_CALU3$sucesses, user_threshold = 0.01, 
  correction_method = "gSCS", custom_bg = nodes_carnival_GSE148729_CALU3$bg,
  sources = c("KEGG","REAC","WP","GO:BP","MIRNA")) 
```

    ## Detected custom background input, domain scope is set to 'custom'

``` r
gostplot(enrichment_results_GSE148729_CALU3, capped = FALSE, interactive = TRUE)
```

<img src="Carnival_virus_signatures_files/figure-gfm/plot_GSE148729_CALU3-1.png" width="800" height="700" />

### Later phase: RIG-I-Like and interferon receptors as perturbation

For a later response to SARS-CoV-2 infection, we set the RIG-I-Like and
the interferon receptors as perturbed.

``` r
RIG_IFN_perturbation <- perturbations_nodes %>%
  dplyr::filter(Mechanism %in% c("RIG","IFN")) %>% 
  dplyr::mutate(sign = 1) %>% 
  dplyr::select(Gene, sign) %>% 
  tibble::column_to_rownames(var = "Gene") %>%
  t() %>% as.data.frame()
```

We also check if these nodes are present in our prior knowledge network
as source of interactions (otherwise it does not make sense to include
its perturbation).

``` r
colnames(RIG_IFN_perturbation) 
```

    ## [1] "DDX58"  "IFIH1"  "DHX58"  "IFNAR1" "IFNAR2" "IFNGR1" "IFNGR2"

``` r
colnames(RIG_IFN_perturbation) %in% all_source_nodes
```

    ## [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE

#### GSE147507.S5\_A549\_SARS.CoV.2

We select this condition and prepared the inputs to run CARNIVAL

``` r
dorothea_results_GSE147507_A549 <- dorothea_results %>% 
  dplyr::select(Tf, GSE147507.S5_A549_SARS.CoV.2) %>% 
  dplyr::filter(Tf %in% all_nodes_network) %>%
  dplyr::top_n(25, wt = abs(GSE147507.S5_A549_SARS.CoV.2)) %>%
  dplyr::arrange(GSE147507.S5_A549_SARS.CoV.2) %>%
  tibble::column_to_rownames(var = "Tf")  %>%
  t() %>% as.data.frame()

progeny_results_GSE147507_A549 <- progeny_results %>% 
  dplyr::select(Pathway,GSE147507.S5_A549_SARS.CoV.2) %>%
  tibble::column_to_rownames(var = "Pathway")  %>%
  t()
```

And we finally run CARNIVAL:

``` r
carnival_results_GSE147507_A549_IFN <-runCARNIVAL(
    solverPath="/home/alvaldeolias/Downloads/cplex",
    netObj=carnival_pkn,
    measObj=dorothea_results_GSE147507_A549,
    inputObj = RIG_IFN_perturbation,
    # dir_name="Carnival_Results",
    weightObj=progeny_results_GSE147507_A549,
    # nodeID = 'gene',
    timelimit = 7200,
    solver = "cplex")
saveRDS(carnival_results_GSE147507_A549_IFN, 
    file = "CARNIVAL_results/viral_signatures/carnival_results_GSE147507_A549_IFN.rds")
OutputCyto(carnival_results_GSE147507_A549_IFN, 
    outputFile="CARNIVAL_results/viral_signatures/carnival_results_GSE147507_A549_IFN")
```

We now run the enrichment analysis:

``` r
carnival_results_GSE147507_A549_IFN <- 
    readRDS("CARNIVAL_results/viral_signatures/carnival_results_GSE147507_A549_IFN.rds")
nodes_carnival_GSE147507_A549_IFN <- 
  extractCARNIVALnodes(carnival_results_GSE147507_A549_IFN)

enrichment_results_GSE147507_A549_IFN <- 
  gost(nodes_carnival_GSE147507_A549_IFN$sucesses, user_threshold = 0.01, 
  correction_method = "gSCS", custom_bg = nodes_carnival_GSE147507_A549_IFN$bg,
  sources = c("KEGG","REAC","WP","GO:BP","MIRNA")) 
```

    ## Detected custom background input, domain scope is set to 'custom'

``` r
gostplot(enrichment_results_GSE147507_A549_IFN, capped = FALSE, interactive = TRUE)
```

<img src="Carnival_virus_signatures_files/figure-gfm/plot_GSE147507_A549_IFN-1.png" width="800" height="700" />

#### GSE147507.S7\_Calu3\_SARS.CoV.2

We select this condition and prepared the inputs to run CARNIVAL

``` r
dorothea_results_GSE147507_CALU3 <- dorothea_results %>% 
  dplyr::select(Tf, GSE147507.S7_Calu3_SARS.CoV.2) %>% 
  dplyr::filter(Tf %in% all_nodes_network) %>%
  dplyr::top_n(25, wt = abs(GSE147507.S7_Calu3_SARS.CoV.2)) %>%
  dplyr::arrange(GSE147507.S7_Calu3_SARS.CoV.2) %>%
  tibble::column_to_rownames(var = "Tf")  %>%
  t() %>% as.data.frame()

progeny_results_GSE147507_CALU3 <- progeny_results %>% 
  dplyr::select(Pathway,GSE147507.S7_Calu3_SARS.CoV.2) %>%
  tibble::column_to_rownames(var = "Pathway")  %>%
  t()
```

And we finally run CARNIVAL:

``` r
carnival_results_GSE147507_CALU3_IFN <-runCARNIVAL(
    solverPath="/home/alvaldeolias/Downloads/cplex",
    netObj=carnival_pkn,
    measObj=dorothea_results_GSE147507_CALU3,
    inputObj = RIG_IFN_perturbation,
    # dir_name="Carnival_Results",
    weightObj=progeny_results_GSE147507_CALU3,
    # nodeID = 'gene',
    timelimit = 7200,
    solver = "cplex")
saveRDS(carnival_results_GSE147507_CALU3_IFN, 
    file = "CARNIVAL_results/viral_signatures/carnival_results_GSE147507_CALU3_IFN.rds")
OutputCyto(carnival_results_GSE147507_CALU3_IFN, 
    outputFile="CARNIVAL_results/viral_signatures/carnival_results_GSE147507_CALU3_IFN")
```

We now run the enrichment analysis:

``` r
carnival_results_GSE147507_CALU3_IFN <- 
    readRDS("CARNIVAL_results/viral_signatures/carnival_results_GSE147507_CALU3_IFN.rds")
nodes_carnival_GSE147507_CALU3_IFN <- 
  extractCARNIVALnodes(carnival_results_GSE147507_CALU3_IFN)

enrichment_results_GSE147507_CALU3_IFN <- 
  gost(nodes_carnival_GSE147507_CALU3_IFN$sucesses, user_threshold = 0.01, 
  correction_method = "gSCS", custom_bg = nodes_carnival_GSE147507_CALU3_IFN$bg,
  sources = c("KEGG","REAC","WP","GO:BP","MIRNA")) 
```

    ## Detected custom background input, domain scope is set to 'custom'

``` r
gostplot(enrichment_results_GSE147507_CALU3_IFN, capped = FALSE, interactive = TRUE)
```

<img src="Carnival_virus_signatures_files/figure-gfm/plot_GSE147507_CALU3_IFN-1.png" width="800" height="700" />

#### GSE148729\_Calu3\_SARS.CoV.2\_24H

We select this condition and prepared the inputs to run CARNIVAL

``` r
dorothea_results_GSE148729_CALU3 <- dorothea_results %>% 
  dplyr::select(Tf, GSE148729_Calu3_SARS.CoV.2_24H) %>% 
  dplyr::top_n(25, wt = abs(GSE148729_Calu3_SARS.CoV.2_24H)) %>%
  dplyr::arrange(GSE148729_Calu3_SARS.CoV.2_24H) %>%
  dplyr::filter(Tf %in% all_nodes_network) %>%
  tibble::column_to_rownames(var = "Tf")  %>%
  t() %>% as.data.frame()

progeny_results_GSE148729_CALU3 <- progeny_results %>% 
  dplyr::select(Pathway,GSE148729_Calu3_SARS.CoV.2_24H) %>%
  tibble::column_to_rownames(var = "Pathway")  %>%
  t()
```

And we finally run CARNIVAL:

``` r
carnival_results_GSE148729_CALU3_IFN <-runCARNIVAL(
    solverPath="/home/alvaldeolias/Downloads/cplex",
    netObj=carnival_pkn,
    measObj=dorothea_results_GSE148729_CALU3,
    inputObj = RIG_IFN_perturbation,
    # dir_name="Carnival_Results",
    weightObj=progeny_results_GSE148729_CALU3,
    # nodeID = 'gene',
    timelimit = 7200,
    solver = "cplex")
saveRDS(carnival_results_GSE148729_CALU3_IFN, 
    file = "CARNIVAL_results/viral_signatures/carnival_results_GSE148729_CALU3_IFN.rds")
OutputCyto(carnival_results_GSE148729_CALU3_IFN, 
    outputFile="CARNIVAL_results/viral_signatures/carnival_results_GSE148729_CALU3_IFN")
```

We now run the enrichment analysis:

``` r
carnival_results_GSE148729_CALU3_IFN <- 
    readRDS("CARNIVAL_results/viral_signatures/carnival_results_GSE148729_CALU3_IFN.rds")
nodes_carnival_GSE148729_CALU3_IFN <- 
  extractCARNIVALnodes(carnival_results_GSE148729_CALU3_IFN)

enrichment_results_GSE148729_CALU3_IFN <- 
  gost(nodes_carnival_GSE148729_CALU3_IFN$sucesses, user_threshold = 0.01, 
  correction_method = "gSCS", custom_bg = nodes_carnival_GSE148729_CALU3_IFN$bg,
  sources = c("KEGG","REAC","WP","GO:BP","MIRNA")) 
```

    ## Detected custom background input, domain scope is set to 'custom'

``` r
gostplot(enrichment_results_GSE148729_CALU3_IFN, capped = FALSE, interactive = TRUE)
```

<img src="Carnival_virus_signatures_files/figure-gfm/plot_GSE148729_CALU3_IFN-1.png" width="800" height="700" />

## Session Info Details

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] plotly_4.9.2.1   gprofiler2_0.1.9 CARNIVAL_1.0.1   OmnipathR_1.3.5 
    ##  [5] jsonlite_1.7.0   igraph_1.2.5     forcats_0.5.0    stringr_1.4.0   
    ##  [9] dplyr_1.0.1      purrr_0.3.4      readr_1.3.1      tidyr_1.1.1     
    ## [13] tibble_3.0.3     ggplot2_3.3.2    tidyverse_1.3.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_1.4-1     ellipsis_0.3.1       class_7.3-17        
    ##   [4] fs_1.5.0             rstudioapi_0.11      bit64_4.0.2         
    ##   [7] AnnotationDbi_1.50.3 fansi_0.4.1          lubridate_1.7.9     
    ##  [10] xml2_1.3.2           codetools_0.2-16     splines_4.0.2       
    ##  [13] doParallel_1.0.15    knitr_1.29           broom_0.7.0         
    ##  [16] annotate_1.66.0      kernlab_0.9-29       dbplyr_1.4.4        
    ##  [19] graph_1.66.0         shiny_1.5.0          compiler_4.0.2      
    ##  [22] httr_1.4.2           backports_1.1.8      fastmap_1.0.1       
    ##  [25] assertthat_0.2.1     Matrix_1.2-18        lazyeval_0.2.2      
    ##  [28] cli_2.0.2            later_1.1.0.1        htmltools_0.5.0     
    ##  [31] tools_4.0.2          gtable_0.3.0         glue_1.4.1          
    ##  [34] Category_2.54.0      rappdirs_0.3.1       Rcpp_1.0.5          
    ##  [37] Biobase_2.48.0       cellranger_1.1.0     vctrs_0.3.2         
    ##  [40] iterators_1.0.12     crosstalk_1.1.0.1    xfun_0.16           
    ##  [43] ps_1.3.4             rvest_0.3.6          mime_0.9            
    ##  [46] lpSolve_5.6.15       lifecycle_0.2.0      XML_3.99-0.5        
    ##  [49] MASS_7.3-52          scales_1.1.1         promises_1.1.1      
    ##  [52] hms_0.5.3            parallel_4.0.2       RBGL_1.64.0         
    ##  [55] yaml_2.2.1           curl_4.3             memoise_1.1.0       
    ##  [58] viper_1.22.0         segmented_1.2-0      stringi_1.4.6       
    ##  [61] RSQLite_2.2.0        genefilter_1.70.0    S4Vectors_0.26.1    
    ##  [64] foreach_1.5.0        e1071_1.7-3          BiocGenerics_0.34.0 
    ##  [67] rlang_0.4.7          pkgconfig_2.0.3      bitops_1.0-6        
    ##  [70] evaluate_0.14        lattice_0.20-41      labeling_0.3        
    ##  [73] htmlwidgets_1.5.1    processx_3.4.3       bit_4.0.4           
    ##  [76] tidyselect_1.1.0     GSEABase_1.50.1      magrittr_1.5        
    ##  [79] R6_2.4.1             IRanges_2.22.2       UniProt.ws_2.28.0   
    ##  [82] generics_0.0.2       DBI_1.1.0            pillar_1.4.6        
    ##  [85] haven_2.3.1          withr_2.2.0          mixtools_1.2.0      
    ##  [88] survival_3.1-12      RCurl_1.98-1.2       modelr_0.1.8        
    ##  [91] crayon_1.3.4         KernSmooth_2.23-17   BiocFileCache_1.12.1
    ##  [94] rmarkdown_2.3        grid_4.0.2           readxl_1.3.1        
    ##  [97] data.table_1.13.0    callr_3.4.3          blob_1.2.1          
    ## [100] webshot_0.5.2        reprex_0.3.0         digest_0.6.25       
    ## [103] xtable_1.8-4         httpuv_1.5.4         stats4_4.0.2        
    ## [106] munsell_0.5.0        viridisLite_0.3.0
