## Considering top candidate among "hypothesis generating" gene sets
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )

synapseLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/" )@filePath }

## Quiet version of read_csv
readq_csv <- function( ... ) { suppressMessages( read_csv(...) ) }

## Loads a file from Synapse and annotates it with the pre-defined region and method
synLoad <- function( id, region, method )
{
    syn( id ) %>% readq_csv %>% select( Drug = id, AUC, p_value ) %>%
        mutate( Region = !!region, Method = !!method ) %>%
        arrange( p_value ) %>% mutate( Rank = 1:nrow(.) )
}

fTest <- function( id, region, method )
{
    paste( id, region, method )
}

main <- function()
{
    ## Define the slice of data relevant to the analysis
    S <- tribble(
        ~Region, ~Method, ~SynID,
        "BM22", "MC-lin", "syn11770897",
        "BM22", "MC-nln", "syn11770904",
        "BM22", "Ordinal", "syn11770910",
        "BM36", "MC-lin", "syn11770921",
        "BM36", "MC-nln", "syn11770936",
        "BM36", "Ordinal", "syn11770959" )

    ## Load all the data
    XX <- S %>% rowwise %>% do( X1 = synLoad( .$SynID, .$Region, .$Method ) ) %>% unnest

    ## Measure consistency of metrics across regions and across methods
    R1 <- XX %>% group_by( Drug, Region ) %>%
        summarize( AveAUC = mean(AUC), MedianRank = median(Rank) ) %>% ungroup %>%
        arrange( MedianRank )

##    R1 %>% group_by( Region ) %>% top_n( -15, MedianRank ) %>% ungroup %>% as.data.frame %>%
##        arrange( Region ) %>% write_tsv( "Top15-by-region.tsv" )

    ## Look at Ordinal regression only
    ORD <- XX %>% filter( Method == "Ordinal" ) %>% group_by( Region ) %>% top_n( -20, Rank ) %>%
        ungroup %>% select( Drug, Rank, Region, everything(), -Method ) %>% as.data.frame

    ## Look at Multi-class Nonlinear only
    NLN <- XX %>% filter( Method == "MC-nln" ) %>% group_by( Region ) %>% top_n( -15, Rank ) %>%
        ungroup %>% select( Drug, Rank, Region, everything(), -Method ) %>% as.data.frame

    ## Join against the meta-data (drug name, nominal target, etc.)
    M <- syn( "syn11801537" ) %>% readq_csv %>%
        rename( URL = link, LINCS_ID = lincs_id, Drug = name, Target = target_name )
    X <- ORD %>% rename( URL = Drug ) %>% inner_join( M, by = "URL" ) %>%
        mutate( Target = ifelse( is.na(Target), "", Target ) )

    ## Join against Steve's annotations
##    S <- read_tsv( "~/Downloads/Top15-old.tsv" ) %>% select( URL = Drug, TargetS = target ) %>%
##        filter( !duplicated( URL ) )
##    X %>% left_join( S ) %>% mutate( TargetS = ifelse( is.na(TargetS), "", TargetS ) )

    ## Fix several target entries by hand
    X <- X %>% mutate( Target = ifelse( LINCS_ID == "HMSL10226", "JAK1,JAK3", Target ) ) %>%
        mutate( Target = ifelse( LINCS_ID == "HMSL10510", "PI3K", Target ) ) %>%
        mutate( Target = ifelse( LINCS_ID == "HMSL10209", "IGF1R", Target ) ) %>%
        mutate( Target = ifelse( LINCS_ID %in% c( "HMSL10470", "HMSL10475" ), "DAGK", Target ) )
    
    write_tsv( X, "Top20-by-region.tsv" )
}

## Extracts a gene set associated with dsRNA
geneset.dsRNA <- function()
{
    ## Retrieve the raw data
    X <- syn( "syn11807753" ) %>% readq_csv %>%
        select( Gene=`Gene Symbol`, Control1=`Control Replicate 1`, Control2=`Control Replicate 2`,
               dsRNAmi1 = `DsRNAmi Replicate 1`, dsRNAmi2 = `DsRNAmi Replicate 2` ) %>%
        rowwise %>% mutate( Control = mean( c(Control1, Control2), na.rm=TRUE ),
                           dsRNAmi = mean( c(dsRNAmi1, dsRNAmi2), na.rm=TRUE ) ) %>% ungroup %>%
        mutate( FoldChange = dsRNAmi / Control )

    ## Corner case: protein not expressed in control leading infinite fold change
    v1 <- X %>% filter( Control == 0, dsRNAmi > 10 ) %>% magrittr::extract2( "Gene" )

    ## Consider remaining data
    v2 <- X %>% filter( Control != 0 ) %>% select( Gene, FoldChange ) %>%
        mutate( lFC = abs(log2( FoldChange )) ) %>% arrange( desc(lFC) ) %>% filter( lFC > 1 ) %>%
        magrittr::extract2( "Gene" )

    ## Combine the two lists and retrieve unique entries
    ## Write the result to a file and store the file to Synapse
    v <- c( v1, v2 ) %>% unique %>% cat( file = "dsRNA.txt", sep="\n" )
    f <- File( "dsRNA.txt", parentId = "syn11629934" )
    annotations(f)$Type <- "Gene Set"
    annotations(f)$Category <- "Interferome"
    annotations(f)$Reference <- "LSP Experiment"
    synStore(f)
}
