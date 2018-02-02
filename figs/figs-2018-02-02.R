## Figures for proposal planning meeting
##
## by Artem Sokolov

library( tidyverse )
library( stringr )
library( synapseClient )
library( grid )
library( gridExtra )

synapseLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/" )@filePath }

## Quiet version of read_csv
readq_csv <- function( ... ) { suppressMessages( read_csv(...) ) }

## Loads all the relevant background matrices
loadBK <- function()
{
    ff <- list()
    ff$BM22.lin <- syn( "syn11727060" )
    ff$BM22.nln <- syn( "syn11727045" )
    ff$BM22.ord <- syn( "syn11726985" )
    ff$BM36.lin <- syn( "syn11727053" )
    ff$BM36.nln <- syn( "syn11727051" )
    ff$BM36.ord <- syn( "syn11736123" )

    lapply( ff, readq_csv )
}

mainBK <- function()
{
    ## Combine all background sets into a single matrix for comparison across brain region and
    ##   across prediction method
    f <- function( X ) { X %>% mutate( Run = 1:nrow(X) ) %>% gather( Size, AUC, -Run ) }
    BB <- loadBK() %>% lapply( f )
    BB$BM22.lin <- BB$BM22.lin %>% mutate( Region = "BM22", Method = "MC-lin" )
    BB$BM22.nln <- BB$BM22.nln %>% mutate( Region = "BM22", Method = "MC-nln" )
    BB$BM22.ord <- BB$BM22.ord %>% mutate( Region = "BM22", Method = "Ordinal" )
    BB$BM36.lin <- BB$BM36.lin %>% mutate( Region = "BM36", Method = "MC-lin" )
    BB$BM36.nln <- BB$BM36.nln %>% mutate( Region = "BM36", Method = "MC-nln" )
    BB$BM36.ord <- BB$BM36.ord %>% mutate( Region = "BM36", Method = "Ordinal" )
    RR <- bind_rows(BB) %>% mutate( Size = as.integer(Size) )

    ## Plot a comparison across methods and across brain regions
    etxt <- function( s ) {element_text( size = s, face="bold" )}
    thm <- theme( axis.text = etxt(10), axis.title = etxt(12),
              legend.text = etxt(10), legend.title = etxt(12),
              strip.text = etxt(12) )
    gg1 <- ggplot( RR, aes( x = Size, y = AUC, color = Method ) ) +
        geom_smooth() + theme_bw() + facet_wrap( ~Region, ncol=2 ) + thm
    gg2 <- ggplot( RR, aes( x = Size, y = AUC, color = Region ) ) +
        geom_smooth() + theme_bw() + facet_wrap( ~Method, ncol = 3 ) + thm

    gg <- arrangeGrob( gg1, gg2, ncol = 1, heights=c(3,2) )
    ggsave( "background_compare.png", gg, width = 8, height = 7 )
}

## Retrieves background distribution for gene sets of requested size
getBK <- function( B, n )
{
    ## Identify the closest column
    j <- as.character( round( n / 10 ) * 10 )
    B[[j]]
}

## Generates a figure for a given AUC / Background / gene set size
make.fig <- function( AUC, B, n, lbl, xlbl = "AUC" )
{
    ## Retrieve the background
    v <- getBK( B, n )
    RR <- data_frame( Bk = v )
    pval <- sum( AUC < v ) / length(v)
    txt <- paste0( "p-value=", round(pval,3) )

    ## Plot everything
    ggplot( RR, aes(x=Bk) ) + theme_bw() +
        geom_density( size = 1, fill = "steelblue", alpha = 0.2 ) + ylab( "Density" ) +
        xlab( xlbl ) +
        geom_vline( xintercept = AUC, color = "red", size=1 ) +
        theme( axis.text = element_text( face="bold", size = 12 ),
              axis.title = element_text( face="bold", size = 12 ) ) +
        geom_text( data = RR[1,], label = txt, x = -Inf, y=Inf, size = 6,
                  vjust = 1.2, hjust = -0.2, fontface="bold" ) +
        ggtitle( lbl )
}

## Plots for Metformin performance
metformin <- function( fg, bg )
{
    ## Name of set -> n_genes -> intersect
    ## SPNetwork -> 477 -> 404
    ## upstream -> 65 -> 54
    ## downstream -> 66 -> 60

    ## Load the foreground and background values
    BB <- syn( bg ) %>% readq_csv
    FF <- syn( fg ) %>% readq_csv

    ## Compose individual figures
    gg1 <- make.fig( FF$`Metformin-upstream.txt`, BB, 54, "Upstream", "" )
    gg2 <- make.fig( FF$`Metformin-SPNetwork.txt`, BB, 404, "SPNetwork" )
    gg3 <- make.fig( FF$`Metformin-downstream.txt`, BB, 60, "Downstream", "" )
    arrangeGrob( gg1, gg2, gg3, nrow=1 )
}

main <- function()
{
    ## BM22, multi-class non-linear
    BM22.nln <- metformin( "syn11727046", "syn11727045" )
    ggsave( "metformin-rf.png", BM22.nln, width=12, height=4 )

    ## BM22, ordinal
    BM22.ord <- metformin( "syn11726984", "syn11726985" )
    ggsave( "metformin.png", BM22.ord, width=12, height=4 )
}
