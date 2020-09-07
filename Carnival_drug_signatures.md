CARNIVAL on activities for selected drugs with anti SARS-CoV-2 effect
================
Alberto Valdeolivas: <alberto.valdeolivas@bioquant.uni-heidelberg.de>;
Date:
31/08/2020

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

The present script takes pathway and tf activity on selected drugs with
an anti SARS-CoV-2 effect.

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

We now read the TF and pathway activity extracted from cell lines where
the drugs with anti SARS-CoV-2 effect was used. We also read a file
containing the target genes of those drugs.

``` r
### Dorothea + Progeny Results 
dorothea_results <- 
  read_csv(file = "CARNIVAL_input/drug_signatures_dorothea.csv") %>% 
  dplyr::rename(Tf = "X1")

progeny_results <- 
  read_csv(file = "CARNIVAL_input/drug_signatures_progeny.csv") %>% 
  dplyr::rename(Pathway = "X1")

### Perturbation nodes
perturbations_nodes <- 
  read_csv(file = "CARNIVAL_input/drug_signatures_pert.csv") 
```

## Results

We now are going to run CARNIVAL for the different drugs.

### Alprostadil

We set the gene targets of **Alprostadil** as perturbed.

``` r
alprostadil_targets_perturbation <- perturbations_nodes %>%
  dplyr::filter(pert_iname == "alprostadil") %>% 
  dplyr::pull(target) %>% 
  stringr::str_split(pattern = "\\|",simplify=FALSE) %>% 
  unlist()

alprostadil_sign_perturbation <- perturbations_nodes %>%
  dplyr::filter(pert_iname == "alprostadil") %>% 
  dplyr::pull(act_or_inh) %>% 
  as.integer() %>% 
  rep(length(alprostadil_targets_perturbation)) %>% 
  as.data.frame() %>% t()

alprostadil_perturbation <- alprostadil_sign_perturbation
colnames(alprostadil_perturbation) <- alprostadil_targets_perturbation
```

We also check if these nodes are present in our prior knowledge network
as source of interactions (otherwise it does not make sense to include
its perturbation).

``` r
colnames(alprostadil_perturbation) 
```

    ## [1] "CATSPER1" "CATSPER2" "CATSPER3" "CATSPER4" "PTGDR"    "PTGER1"   "PTGER2"  
    ## [8] "PTGER4"   "PTGIR"

``` r
colnames(alprostadil_perturbation) %in% all_source_nodes
```

    ## [1] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE

We select the dorothea and progeny results containing the TF and pathway
activity effect of this drug.

``` r
dorothea_results_alprostadil <- dorothea_results %>% 
  dplyr::select(Tf, alprostadil) %>% 
  dplyr::filter(Tf %in% all_nodes_network) %>%
  dplyr::top_n(25, wt = abs(alprostadil)) %>%
  dplyr::arrange(alprostadil) %>%
  tibble::column_to_rownames(var = "Tf")  %>%
  t() %>% as.data.frame()

progeny_results_alprostadil <- progeny_results %>% 
  dplyr::select(Pathway,alprostadil) %>%
  tibble::column_to_rownames(var = "Pathway")  %>%
  t()
```

And we finally run CARNIVAL:

``` r
carnival_results_alprostadil <-runCARNIVAL(
    solverPath="/home/alvaldeolias/Downloads/cplex",
    netObj=carnival_pkn,
    measObj=dorothea_results_alprostadil,
    inputObj = alprostadil_perturbation,
    # dir_name="Carnival_Results",
    weightObj=progeny_results_alprostadil,
    # nodeID = 'gene',
    timelimit = 100,
    solver = "cplex")
saveRDS(carnival_results_alprostadil, 
    file = "CARNIVAL_results/drug_signatures/carnival_results__alprostadil")
OutputCyto(carnival_results_GSE147507_A549, 
    outputFile="CARNIVAL_results/drug_signatures/carnival_results__alprostadil")
```

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
    ##  [1] segmented_1.2-0      Category_2.54.0      bitops_1.0-6        
    ##  [4] fs_1.5.0             lubridate_1.7.9      bit64_4.0.2         
    ##  [7] doParallel_1.0.15    httr_1.4.2           tools_4.0.2         
    ## [10] backports_1.1.8      R6_2.4.1             KernSmooth_2.23-17  
    ## [13] lazyeval_0.2.2       DBI_1.1.0            BiocGenerics_0.34.0 
    ## [16] colorspace_1.4-1     withr_2.2.0          tidyselect_1.1.0    
    ## [19] curl_4.3             bit_4.0.4            compiler_4.0.2      
    ## [22] graph_1.66.0         cli_2.0.2            rvest_0.3.6         
    ## [25] Biobase_2.48.0       xml2_1.3.2           scales_1.1.1        
    ## [28] genefilter_1.70.0    RBGL_1.64.0          rappdirs_0.3.1      
    ## [31] digest_0.6.25        mixtools_1.2.0       rmarkdown_2.3       
    ## [34] pkgconfig_2.0.3      htmltools_0.5.0      dbplyr_1.4.4        
    ## [37] htmlwidgets_1.5.1    rlang_0.4.7          readxl_1.3.1        
    ## [40] rstudioapi_0.11      RSQLite_2.2.0        generics_0.0.2      
    ## [43] viper_1.22.0         RCurl_1.98-1.2       magrittr_1.5        
    ## [46] Matrix_1.2-18        Rcpp_1.0.5           munsell_0.5.0       
    ## [49] S4Vectors_0.26.1     fansi_0.4.1          lifecycle_0.2.0     
    ## [52] UniProt.ws_2.28.0    stringi_1.4.6        yaml_2.2.1          
    ## [55] MASS_7.3-52          BiocFileCache_1.12.1 grid_4.0.2          
    ## [58] blob_1.2.1           parallel_4.0.2       crayon_1.3.4        
    ## [61] lattice_0.20-41      haven_2.3.1          splines_4.0.2       
    ## [64] annotate_1.66.0      hms_0.5.3            knitr_1.29          
    ## [67] pillar_1.4.6         lpSolve_5.6.15       codetools_0.2-16    
    ## [70] stats4_4.0.2         XML_3.99-0.5         reprex_0.3.0        
    ## [73] glue_1.4.1           evaluate_0.14        data.table_1.13.0   
    ## [76] modelr_0.1.8         vctrs_0.3.2          foreach_1.5.0       
    ## [79] cellranger_1.1.0     gtable_0.3.0         kernlab_0.9-29      
    ## [82] assertthat_0.2.1     xfun_0.16            xtable_1.8-4        
    ## [85] broom_0.7.0          e1071_1.7-3          viridisLite_0.3.0   
    ## [88] class_7.3-17         survival_3.1-12      iterators_1.0.12    
    ## [91] AnnotationDbi_1.50.3 memoise_1.1.0        IRanges_2.22.2      
    ## [94] ellipsis_0.3.1       GSEABase_1.50.1
