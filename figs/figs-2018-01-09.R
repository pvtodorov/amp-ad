## Figures for R56 call
##
## by Artem Sokolov

library( tidyverse )
library( stringr )
library( synapseClient )

synapseLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/" )@filePath }

## Quiet version of read_csv
readq_csv <- function( ... ) { suppressMessages( read_csv(...) ) }

## Retrieves background distribution for gene sets of requested size
getBK <- function( n )
{
    ## Identify the closest column
    j <- as.character( round( n / 10 ) * 10 )
    Xbk <- syn( "syn11645640" ) %>% readq_csv
    Xbk[[j]]
}

## Generates a figure for a given R^2 / gene set size
make.fig <- function( R2, n )
{
    ## Retrieve the background
    v <- getBK( n )
    RR <- data_frame( Bk = v )
    pval <- sum( R2 < v ) / length(v)
    txt <- paste0( "p-value=", round(pval,3) )

    ## Plot everything
    ggplot( RR, aes(x=Bk) ) + theme_bw() +
        geom_density( size = 1, fill = "steelblue", alpha = 0.2 ) + ylab( "Density" ) +
        xlab( "Performance estimate for predicting Braak score from mRNA (R^2)" ) +
        geom_vline( xintercept = R2, color = "red", size=1 ) +
        theme( axis.text = element_text( face="bold", size = 12 ),
              axis.title = element_text( face="bold", size = 12 ) ) +
        geom_text( data = RR[1,], label = txt, x = -Inf, y=Inf, size = 6,
                  vjust = 1.2, hjust = -0.2, fontface="bold" )
}

## Microglia
figs.microglia <- function()
{
    X <- syn( "syn11664081" ) %>% read_csv

    gg6 <- make.fig( X$R2[1], X$intersect[1] )
    ggsave( "clus6.png", gg6, width=6, height=4 )

    gg3 <- make.fig( X$R2[2], X$intersect[2] )
    ggsave( "clus3.png", gg3, width=6, height=4 )

    gg7 <- make.fig( X$R2[3], X$intersect[3] )
    ggsave( "clus7.png", gg7, width=6, height=4 )
}

## Metformin
figs.metformin <- function()
{
    X <- syn( "syn11664078" ) %>% read_csv

    fns <- c( "met-spn.png", "met-up.png", "met-down.png" )
    for( i in 1:nrow(X) )
    {
        gg <- make.fig( X$R2[i], X$intersect[i] )
        ggsave( fns[i], gg, width=6, height=4 )
    }
}

## Interferome
figs.isg <- function()
{
    X <- syn( "syn11664080" ) %>% read_csv
    gg <- make.fig( X$R2, X$intersect )
    ggsave( "isg.png", gg, width=6, height=4 )
}

## Barplot: small molecular pathways database
barp.smpdb <- function()
{
    ## Load the data and extract pathway names
    X <- syn( "syn11658930" ) %>% readq_csv
    f <- function( v ) { str_split( v, ";", simplify=TRUE )[,1] %>% str_sub( 7 ) }
    X <- X %>% mutate( name = f(description) ) %>% select( name, R2, p_value, intersect, adjusted_p ) %>%
        group_by( name ) %>% slice( 1 ) %>% ungroup %>% arrange( desc(R2) )

    ## Make the barplot
    Y <- head(X,20) %>% mutate( name = factor(name, rev(name)) )
    gg <- ggplot( Y, aes( x = name, y = R2 ) ) + theme_bw() +
        geom_bar( stat = "identity" ) + coord_flip() + xlab("") +
        theme( axis.title = element_text( size = 14, face = "bold" ),
              axis.text = element_text( size = 12, face = "bold" ) )
    ggsave( "smpdb.png", gg, width = 8, height = 6 )

    ## Create a background plot for the top hit
    gg1 <- make.fig( Y$R2[1], Y$intersect[1] )
    ggsave( "smpdb1.png", gg1, width = 6, height = 4 )
}

## Barplot: Nienke's sets
barp.nienke <- function()
{
    ## Load the data and extract pathway names
    X <- syn( "syn11665709" ) %>% readq_csv
    gg <- make.fig( X$R2[1], X$intersect[1] )
    ggsave( "BMS345541.png", gg, width = 6, height = 4 )
}

## Barplot: Hallmarks
barp.hallmarks <- function()
{
    ## Load the data and extract pathway names
    X <- syn( "syn11658944" ) %>% readq_csv %>% mutate( name = str_sub( description, 10 ) )
    Y <- X %>% arrange( desc(R2) ) %>% head( 20 ) %>% mutate( name = factor( name, rev(name) ) )
    gg <- ggplot( Y, aes( x = name, y = R2 ) ) + theme_bw() +
        geom_bar( stat = "identity" ) + coord_flip() + xlab("") +
        theme( axis.title = element_text( size = 14, face = "bold" ),
              axis.text = element_text( size = 12, face = "bold" ) )
    ggsave( "hallmarks.png", gg, width = 8, height = 6 )

    ## Make a background plot for the top hit
    gg1 <- make.fig( Y$R2[1], Y$intersect[1] )
    ggsave( "il2-stat5.png", gg1, width = 6, height = 4 )
}
