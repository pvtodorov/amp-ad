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

fig2c <- function()
{
    ## Load dsRNA and Metformin results
    load( "res-SPN-EvI.RData" )
    load( "res-dsRNA-EvI.RData" )

    ## Compute summary statistics across the different metaparameter values
    vSPN <- sapply( rSPN, max )
    vdsRNA <- sapply( rdsRNA, max )

    ## Put together foreground and background data frames
    FG <- data_frame( AUC = c( vSPN[1], vdsRNA[1] ),
                     GeneSet = c( "Metformin", "IFN-I" ),
                     SetSize = str_c( "#Genes: ", c(472, 292) ),
                     pLbl = c( "p<0.01", "p=0.03" ) ) %>%
        mutate( Lbl = str_c( SetSize, "\n  ", pLbl ) )
    BK <- bind_rows( data_frame( Value = vSPN[-1], Bk = str_c( "Background", 1:(length(vSPN)-1) ),
                                GeneSet = "Metformin" ),
                    data_frame( Value = vdsRNA[-1], Bk = str_c( "Background", 1:(length(vdsRNA)-1) ),
                               GeneSet = "IFN-I" ) )
                                
    ## Plot everything
    ggplot( BK, aes(x=Value) ) + theme_bw() +
        geom_density( size = 1, fill = "steelblue", alpha = 0.2 ) + ylab( "Density" ) +
        xlab( "AUC" ) + xlim( c(0.44, 0.66) ) + ylim( 0, 14 ) +
        facet_wrap( ~GeneSet, ncol=1 ) +
        theme( strip.background = element_blank(), strip.text = element_blank(),
              axis.text.x = element_text( face="bold", size=12),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title = element_text( face="bold", size=14 ),
              plot.margin = unit( c(5.5,5.5,5.5,30), "pt" ) ) +
        geom_vline( data = FG, aes(xintercept = AUC), color = "red", size=1 ) +
        geom_text( data = FG, aes(x=AUC, label=GeneSet), y=Inf, color="red",
                  hjust = 1.1, vjust = 1.5, fontface="bold" ) +
        geom_text( data = FG, aes(label=Lbl), x=-Inf, y=Inf, hjust=-0.1, vjust=1.25,
                  fontface="bold" )
}

fig2 <- function()
{
    gga <- fig2a()
    ggb <- fig2b()
    ggc <- fig2c()

    gg <- arrangeGrob( gga, ggb, ggc, widths = c(2,2,1.75), nrow=1 )
##    grid.draw(gg)
    ggsave( "fig2.png", gg, width = 11, height = 3 )
}

