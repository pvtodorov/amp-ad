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
