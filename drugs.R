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
    syn( id ) %>% readq_csv %>% t %>% as.data.frame %>% rownames_to_column( "Drug" ) %>%
        rename( AUC = V1 ) %>% mutate( Region = !!region, Method = !!method ) %>%
        arrange( desc(AUC) ) %>% mutate( Rank = 1:nrow(.) )
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
        "BM22", "MC-lin", "syn11766888",
        "BM22", "MC-nln", "syn11766890",
        "BM22", "Ordinal", "syn11766891",
        "BM36", "MC-lin", "syn11766896",
        "BM36", "MC-nln", "syn11766894",
        "BM36", "Ordinal", "syn11766893" )

    ## Load all the data
    XX <- S %>% rowwise %>% do( X1 = synLoad( .$SynID, .$Region, .$Method ) ) %>% unnest

    ## Measure consistency of metrics across regions and across methods
    R1 <- XX %>% group_by( Drug, Region ) %>%
        summarize( AveAUC = mean(AUC), MedianRank = median(Rank) ) %>% ungroup %>%
        arrange( MedianRank )

    R1 %>% group_by( Region ) %>% top_n( -15, MedianRank ) %>% ungroup %>% as.data.frame %>%
        arrange( Region ) %>% write_tsv( "Top15-by-region.tsv" )
}
