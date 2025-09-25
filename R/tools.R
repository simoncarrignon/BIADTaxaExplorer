split_plus <- function(listtax) trimws(strsplit(listtax,"\\+")[[1]])

groupTaxons <- function(dataset,filename="groupings/periods/group_IV.csv"){
    gdata <- read.csv(filename,skip=1,header=F)
    gdata  <- gdata[!is.na(gdata[,1]),]
    txl <- trimws(gsub("=","",gdata[,2]))
    lot  <-  sapply(gdata[,3],split_plus)
    names(lot)=txl
    grouping=lot
    correspTx = unlist(lapply(names(grouping),function(e){names(grouping[[e]])=grouping[[e]];grouping[[e]][]=e;grouping[[e]]}))
    dataset$new_txgroups = correspTx[dataset$TaxonCode]
    print(paste0("grouping ",length(unique(dataset$TaxonCode))," Taxons in  ",length(grouping)," larger categories"))
    return(dataset)
}


#' @param filename, a csv file with list of period and group like:
#'  Period,PeriodID
#'  Neolithic,ELN+EMN+EN+FN+LN+MLN+MN+PPN+PPNB+UN
#'  Eneolithic,EE+ELE+EME+FE+LE+LFE+ME+MLE+PE+PEEE+PEEEME+UE
#'  Bronze Age,EBA+EMBA+FBA+LBA+LFBA+MBA+MLBA+UBA
#'  Early Iron Age,EIA

groupPeriod <- function(originalperiod,filename="groupings/periods/periods.csv"){
    periods <- read.csv(filename)
    grouping=strsplit(periods$PeriodID,"\\+")
    names(grouping) <-  periods$Period
    corresp  <-  unlist(lapply(names(grouping),function(e){names(grouping[[e]])=grouping[[e]];grouping[[e]][]=e;grouping[[e]]}))
    print(paste0("grouping ",length(unique(originalperiod))," periods in  ",length(grouping)," larger categories"))
    return(corresp[originalperiod])
}
