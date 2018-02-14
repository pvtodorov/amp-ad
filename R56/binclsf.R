## Binary classification: early-vs-late and early-vs-intermediate
##
## by Artem Sokolov

library( tidyverse )
library( synapseClient )
library( caret )

synapseLogin()

syn <- function( id )
{ synGet( id, downloadLocation = "~/data/AMP-AD/" )@filePath }

## Quiet versions of read_csv and read_tsv
readq_csv <- function( ... ) { suppressMessages( read_csv(...) ) }
readq_tsv <- function( ... ) { suppressMessages( read_tsv(...) ) }

## Maps Braak score to disease stage (early: 0-2, intermediate: 3-4, late: 5-6 )
braak2stage <- function( v )
{
    ## Define the mapping
    m <- c( rep("early",3), rep("intermediate",2), rep("late",2) )

    ## Use Braak+1 as 1-based indexing into the mapping
    m[v+1]
}

## Evaluates a single subset of data in cross-validation
eval1 <- function( X1 )
{
    ## Set up a cross-validation scheme
    fc <- trainControl( method="lgocv", classProbs = TRUE,
                       summaryFunction = twoClassSummary )
    
    ## Train a predictor
    mm <- train( Stage ~ ., data = X1, method="glmnet", metric="ROC", trControl = fc )
    mm$results$ROC
}

## Evaluates a single feature set in cross-validation
## All features in fSet must be present in data
evalSet <- function( fSet, X )
{
    stopifnot( all( fSet %in% colnames(X) ) )
    X %>% select( Stage, one_of(fSet) ) %>% eval1
}

## Evaluates a given gene set against nBk background sets of the same cardinality
evalGeneSet <- function( X, gSet, nBk = 30 )
{
    ## Determine which genes are present in the data
    vTrue <- intersect( gSet, colnames(X) )
    nTrue <- length(vTrue)

    ## Compose background sets
    lBk <- list()
    for( i in 1:nBk )
        lBk[[i]] <- sample( setdiff(colnames(X), "Stage"), nTrue )
    
    ## Prepend the true set and assign names
    lSet <- c( list(vTrue), lBk ) %>% setNames( c("True", paste0( "Background", 1:nBk ) ) )

    ## Evaluate each set
    res <- list()
    for( i in 1:length(lSet) )
    {
        cat( "Running on", names(lSet)[i], "\n" )
        res[[i]] <- evalSet( lSet[[i]], X )
    }
    setNames( res, names(lSet) )
}

## Summarizes a results list
sRes <- function( ll )
{
    v1 <- lapply( ll, mean ) %>% unlist
    v2 <- lapply( ll, max ) %>% unlist
    f <- function( v ) { sum(v[1] < v[-1]) / (length(v)-1)}
    c( "mean" = v1[1], "max" = v2[1],
        "p-mean" = f(v1), "p-max" = f(v2) )
}

## Evaluates a collection of gene sets
evalSets <- function( X )
{
    ## Load relevant gene sets
    vSPN <- syn( "syn11617550" ) %>% scan( what=character() )
    vUps <- syn( "syn11617551" ) %>% scan( what=character() )
    vDns <- syn( "syn11617549" ) %>% scan( what=character() )
    vISG <- syn( "syn11629935" ) %>% scan( what=character() )
    vdsRNA <- syn( "syn11808176" ) %>% scan( what=character() )

    ## Evaluate each set
    rUps <- evalGeneSet( X, vUps )
    rDns <- evalGeneSet( X, vDns )
    rSPN <- evalGeneSet( X, vSPN )
    rISG <- evalGeneSet( X, vISG )
    rdsRNA <- evalGeneSet( X, vdsRNA )

    ## Combine all the results
    bind_rows( sRes(rUps), sRes(rDns), sRes(rSPN), sRes(rISG), sRes(rdsRNA) ) %>%
        mutate( "Set" = c("Upstream","Downstream","SPN", "ISG", "dsRNA") ) %>%
        select( Set, everything() )
}

## Limit the analysis to BM22 area of MSBB data
main.BM22 <- function()
{
    ## Load MSBB data
    ## Map Braak score to disease stage (Early, Intermediate, Late)
    Xraw <- syn( "syn11615752" ) %>% readq_tsv
    XX <- Xraw %>% na.omit %>%
        mutate( Stage = braak2stage( Braak ) ) %>% filter( BrodmannArea == "BM22" ) %>%
        select( -(PMI:Barcode) ) %>% select( ID, Stage, everything() )

    ## Split up the data into early-vs-intermediate and early-vs-late tasks
    cEvI <- c( "early", "intermediate" )
    cEvL <- c( "early", "late" )
    X.EvI <- XX %>% filter( Stage %in% cEvI ) %>% mutate( Stage = factor( Stage, cEvI ) )
    X.EvL <- XX %>% filter( Stage %in% cEvL ) %>% mutate( Stage = factor( Stage, cEvL ) )

    res1 <- evalSets( X.EvL )
    res2 <- evalSets( X.EvI )

    ## Example dataset for prototyping
    vSPN <- syn( "syn11617550" ) %>% scan( what=character() )
    X1 <- X.EvL %>% select( Stage, one_of( intersect(vSPN, colnames(.)) ) )
}

## Limit the analysis to BM22 area of MSBB data
main.BM36 <- function()
{
    ## Load MSBB data
    ## Map Braak score to disease stage (Early, Intermediate, Late)
    Xraw <- syn( "syn11615752" ) %>% readq_tsv
    XX <- Xraw %>% na.omit %>%
        mutate( Stage = braak2stage( Braak ) ) %>% filter( BrodmannArea == "BM36" ) %>%
        select( -(PMI:Barcode) ) %>% select( ID, Stage, everything() )

    ## Split up the data into early-vs-intermediate and early-vs-late tasks
    cEvI <- c( "early", "intermediate" )
    cEvL <- c( "early", "late" )
    X.EvI <- XX %>% filter( Stage %in% cEvI ) %>% mutate( Stage = factor( Stage, cEvI ) )
    X.EvL <- XX %>% filter( Stage %in% cEvL ) %>% mutate( Stage = factor( Stage, cEvL ) )

    res1 <- evalSets( X.EvI )
    res2 <- evalSets( X.EvL )

    vISG <- syn( "syn11629935" ) %>% scan( what=character() )
    rISG <- evalGeneSet( X.EvL, vISG, 100 )
}

main.ROSMAP <- function()
{
    ## Load ROSMAP data
    ## Map Braak score to disease stage (Early, Intermediate, Late)
    Xraw <- syn( "syn11738012" ) %>% readq_tsv
    XX <- Xraw %>% mutate( Stage = braak2stage( Braak ) ) %>%
        select( -(PMI:Barcode) ) %>% select( ID, Stage, everything() )

    ## Split up the data into early-vs-intermediate and early-vs-late tasks
    cEvI <- c( "early", "intermediate" )
    cEvL <- c( "early", "late" )
    X.EvI <- XX %>% filter( Stage %in% cEvI ) %>% mutate( Stage = factor( Stage, cEvI ) )
    X.EvL <- XX %>% filter( Stage %in% cEvL ) %>% mutate( Stage = factor( Stage, cEvL ) )

    res1 <- evalSets( X.EvI )
    res2 <- evalSets( X.EvL )

    vSPN <- syn( "syn11617550" ) %>% scan( what=character() )
    vdsRNA <- syn( "syn11808176" ) %>% scan( what=character() )

    ## In-depth early-vs-intermediate experiments
    rSPN <- evalGeneSet( X.EvI, vSPN, 100 )
    rdsRNA <- evalGeneSet( X.EvI, vdsRNA, 100 )

    save( rSPN, file="figs/res-SPN-EvI.RData" )
    save( rdsRNA, file="figs/res-dsRNA-EvI.RData" )    

    ## In-depth early-vs-late experiments
    rSPN <- evalGeneSet( X.EvL, vSPN, 100 )
    rdsRNA <- evalGeneSet( X.EvL, vdsRNA, 100 )

    save( rSPN, file="figs/res-SPN-EvL.RData" )
    save( rdsRNA, file="figs/res-dsRNA-EvL.RData" )

    ## Dataset slice for prototyping / variance assessment
    X1 <- X.EvL %>% select( Stage, one_of( intersect(vSPN, colnames(.)) ) )
    eval1( X1 )
}
