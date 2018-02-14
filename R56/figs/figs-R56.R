## Figures for R56
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

## Loads a single background matrix
synLoad <- function( id, region, method )
{
    syn(id) %>% readq_csv %>% mutate( Run = 1:nrow(.) ) %>% gather( Size, AUC, -Run ) %>%
        mutate( Region = !!region, Method = !!method, Size = as.integer(Size) )
}

## Loads all the relevant background matrices
loadBK <- function()
{
    tribble( ~Region, ~Method, ~SynID,
            "BM22", "Multiclass", "syn11727060",
            "BM22", "Non-linear", "syn11727045",
            "BM22", "Ordinal", "syn11726985",
            "BM36", "Multiclass", "syn11727053",
            "BM36", "Non-linear", "syn11727051",
            "BM36", "Ordinal", "syn11736123" ) %>%
        rowwise %>% do( X1 = synLoad( .$SynID, .$Region, .$Method ) ) %>% unnest
}

etxt <- function( s ) {element_text( size = s, face="bold" )}

fig2a <- function()
{
    RR <- loadBK() %>% filter( Size <= 500, Method != "Non-linear", Region == "BM22" )

    ## Plot a comparison across methods and across brain regions
    thm <- theme( axis.text = etxt(12), axis.title.y = etxt(14), axis.title.x = etxt(12),
              legend.text = etxt(12), legend.title = etxt(14),
              strip.text = etxt(14), legend.position = c(0.7,0.8),
              legend.background = element_rect( color="black") )
    ggplot( RR, aes( x = Size, y = AUC, color = Method ) ) +
        xlab( "Number of randomly-selected genes" ) +
        geom_smooth( size = 2, level = 0.99 ) + theme_bw() + thm
}

fig2b <- function()
{
    ## Match against set sizes
    S <- syn( "syn11770910" ) %>% readq_csv() %>% select( URL = id, Size = intersect )
    X <- read_tsv( "Top20-by-region.tsv" ) %>% filter( Region == "BM22", p_value < 0.06 ) %>%
        select( URL, Drug, AUC ) %>% arrange( AUC ) %>% mutate( Drug = factor(Drug, Drug) ) %>%
        inner_join( S )

    ggplot( X, aes( y = AUC, x = Drug )  ) + coord_flip( ylim=c(0.5,0.7) ) + theme_bw() +
        geom_bar( stat = "identity", fill = "tomato", color="black", width = .5 ) +
        geom_text( aes( label = Size ), hjust=-0.5, fontface = "bold" ) +
        theme( axis.title = etxt(14), axis.text = etxt(12),
              plot.margin = unit( c(5.5,5.5,5.5,30), "pt" ) )
}

fig2 <- function()
{
    gga <- fig2a()
    ggb <- fig2b()

    gg <- arrangeGrob( gga, ggb, widths = c(2,2), nrow=1 )
    grid.draw(gg)
    ggsave( "fig2.png", gg, width = 8, height = 3 )
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
