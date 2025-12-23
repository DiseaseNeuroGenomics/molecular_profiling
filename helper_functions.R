library(parallel)
library(GeneOverlap)
library(reshape2)
library(fdrtool)
library(parallel)

##################################################################
# Config
##################################################################

ROOT = "~/Desktop/molecular_profiling/molecular_profiling/" # !!! FIXME: SET TO YOUR CUSTOM DIRECTORY !!
QC_ATACSEQ = file.path(ROOT, "inputs", "")  # Pre-calculated QC metrics for ATAC-seq samples from processing computational pipeline

GENOME_SIZE_hg38 = 2923715862  # non-gap non blacklisted
chrSizesFile_hg38 = file.path(ROOT, "inputs", "hg38-filtered_alt_order.chrom.size")
MY_BLACKLIST_hg38 = file.path(ROOT, "inputs", "hg38.blacklist.bed")
ENSEMBL_INFO_hg38 = file.path(ROOT, "inputs", "muchEnsemblInfo_hg38.tsv.gz")
myGenesFallback_hg38 = file.path(ROOT, "inputs", "hg38genes")

mtsv=function(x,outDir,myHeader=T,gz=F,fileBaseName=deparse(substitute(x))){
  myFile=paste0(outDir, "/",  fileBaseName,".tsv",ifelse(gz,".gz",""))
  if(gz) myFile=gzfile(myFile, "w")
  write.table(x,file=myFile,sep="\t",quote=F,row.names=F,col.names=myHeader)
  if(gz) close(myFile)
}

##################################################################
# Do tests (background version)
##################################################################
# 
# Input types:
# - peaks + background(s)
# - genes + background(s)
# - lists of peaks + background
# - list of genes + background(s)
# - full output of topTable for peaks (must contain the "myGRangesCols" as shown below as well as "adj.P.Val")
# - full output of topTable for genes (must contain the "myGRangesCols" as shown below as well as "adj.P.Val", though the position isn't used and "PeakID" must contain the ensembl ID, which is contraintuitive)
# - lists of full outputs of topTable for peaks with or without background(s)
# - lists of full outputs of topTable for genes with or without background(s)
#

universalGsea=function(
    testMethod,      # one of cameraPR, or fisher #or binominal?
    inputForTest,    # multiple different input types are possible see text above
    myDataType,      # "genes" or "peaks" with the latter being any kind of epigenetic or other positional data
    GENOME_VERSION,  # here we don't provide a default, but currently hg19/hg38/mm10 are supported
    inputForTestMetadata = NULL, # optional metadata, must contain collumns "Set","SetFullName" or "peakSets","peaksFullName"
    background=NULL, # "gene"-data ==> ensemblIDs or lists thereof. "peak"-data ==>  data.frame coersible to GRanges or GRanges or a list of either (lists currently not fully implemented). must have col "PeakID"
    testSetToBackgroundMapping=NULL,
    geneRawMetaSets=NULL,
    gseaExcludeChr=c("chrM","chrY"),
    myGenesFallback=NA, #if NA then a default hg19/hg38/mm10 is picked below
    myGeneMetaSets=NULL,
    furtherArgs__geneSetToPeakSet=NULL, # optional further args to geneSetToPeakSet. Only for "peak"-data
    furtherArgs__geneRegDomains=NULL,   # optional further args to geneRegDomains. Only for "peak"-data
    furtherArgs__cameraPR=NULL,         # optional further args to cameraPR CURRENTLY NOT IMPLEMENTED!
    qValueForFisher=0.05,               # if testMethod is "fisher" and input is full output of topTable, then include peaks with p-val less than this in test set
    outDir=NULL,                        # for plots etc.
    forceNoOutDirCheck=F,
    combineAcrossMultipleBackgrounds=T, # if the results are from multiple backgrounds, combine them
    useRankingForCameraPR=F,            # use rankings instead of the statistic. if TRUE then either a "t" or "F" col must be in the topTable must be present
    useAbsTinRankingForCameraPR=F,      # use absolute t-stat instead of just t-stat if col 
    shrinkOutput=F
){
  cleanNull=function(x)Filter(Negate(is.null),x) #discard items of named list that are simply NULL
  testArg=function(x,acceptableInput=c(T,F)){if(length(x)!=1 | any(!x %in% acceptableInput)) stop(paste0("unexpected input to ", deparse(substitute(x)), ". It must be one of the following:", paste(acceptableInput,collapse=", ")))} #test if input argument is one in length and a prespecified type
  myGRangesCols=c("seqnames", "start", "end", "strand", "PeakID")
  testArg(testMethod,c("fisher","cameraPR"))
  testArg(myDataType,c("genes","peaks"))
  testArg(GENOME_VERSION,c("hg19","hg38","mm10"))
  testArg(forceNoOutDirCheck)
  testArg(combineAcrossMultipleBackgrounds)
  testArg(useRankingForCameraPR)
  testArg(useAbsTinRankingForCameraPR)
  testArg(shrinkOutput)
  if(is.null(background) & !is.null(testSetToBackgroundMapping)) stop("you cannot provide a testSetToBackgroundMapping if no background(s) are provided")
  if(is.null(getOption("mc.cores"))) message("Not using parallel processing as option >>mc.cores<< is not specified!")
  
  w=list()
  w$testMethod=testMethod; rm(testMethod)
  w$myDataType=myDataType; rm(myDataType)
  w$GENOME_VERSION=GENOME_VERSION; rm(GENOME_VERSION)
  w$background=background; rm(background)
  w$testSetToBackgroundMapping=testSetToBackgroundMapping; rm(testSetToBackgroundMapping)
  w$myGeneMetaSets=myGeneMetaSets; rm(myGeneMetaSets)
  w$combineAcrossMultipleBackgrounds=combineAcrossMultipleBackgrounds; rm(combineAcrossMultipleBackgrounds)
  w$useRankingForCameraPR=useRankingForCameraPR; rm(useRankingForCameraPR)
  w$useAbsTinRankingForCameraPR=useAbsTinRankingForCameraPR; rm(useAbsTinRankingForCameraPR)
  w$inputForTestMetadata=inputForTestMetadata;rm(inputForTestMetadata)
  
  #address metadata
  if(!is.null(w$inputForTestMetadata)){
    #These two substitutions are for backwards compatability with step4 metadata
    colnames(w$inputForTestMetadata)=gsub("^peakSets$","Set",colnames(w$inputForTestMetadata))
    colnames(w$inputForTestMetadata)=gsub("^peaksFullName$","SetFullName",colnames(w$inputForTestMetadata))
    if(!all(c("Set","SetFullName") %in% colnames(w$inputForTestMetadata))) stop("required column in inputForTestMetadata is missing")
    #NOTE: one could implement better testing between provided inputForTest and inputForTestMetadata
  }
  
  #grab gene list with positions
  if(is.na(myGenesFallback)){ myGenesFallback = myGenesFallback_hg38 } #NOTE: global vars
  w$myGenesTable=read.delim(myGenesFallback,stringsAsFactors=F,col.names=c("chr","tss","strand","ensembl"))
  w$myGenesTable=w$myGenesTable[!w$myGenesTable$chr %in% gseaExcludeChr,]
  w$myGenes=w$myGenesTable$ensembl
  message("grabbed gene list with positions")
  
  #determine input type
  if(!is.list(inputForTest) & !is.data.frame(inputForTest)){
    message(paste("Input is interpreted simply as a vector of",w$myDataType))
    w$inputForTest=list(myPeaks=inputForTest)
    w$inputType="vectorOfItems"
    w$inputCount="single"
  }else if(is.list(inputForTest)){
    if(is.data.frame(inputForTest)){
      message(paste("Input is interpreted to be the output of limma's topTable containing",w$myDataType))
      w$inputForTest=list(myPeaks=inputForTest)
      w$inputType="topTable"
      w$inputCount="single"
    }else if(  all(!sapply(inputForTest,is.list))  &  all(sapply(inputForTest,is.vector)) ){
      message(paste("Input interpreted to be a list of vectors of",w$myDataType))
      w$inputForTest=inputForTest
      w$inputType="vectorOfItems"
      w$inputCount="multiple"
    }else if(all(sapply(inputForTest,is.data.frame))){
      message(paste("Input interpreted to be a list of topTables of",w$myDataType))
      w$inputForTest=inputForTest
      w$inputType="topTable"
      w$inputCount="multiple"
    }else{
      stop("unexpected input to universalGsea")
    }
  }else{
    stop("unexpected input to universalGsea")
  }
  rm(inputForTest)
  
  
  #test parameters for vectorOfItems
  if(w$inputType=="vectorOfItems"){
    if(!w$testMethod %in% c("fisher")) stop("only test method supported for vector(s) of peaks/genes is >>fisher<<")
    if(any(sapply(w$inputForTest,anyDuplicated)>0)) stop("duplicated peaks/genes in input for test")
    if(w$inputCount=="multiple")
      if(anyDuplicated(names(w$inputForTest)))stop("duplicated gene set names in input for test")
    #background
    if(is.null(w$background) & w$myDataType=="peaks") stop("if input is vector(s) of peaks, a background set must be provided")
    if(is.null(w$background) & w$myDataType=="genes"){
      message("used >>myGenesFallback<< to generate background. Please make sure this is what you want, as an impropriate background can bias results")
      w$background=w$myGenes
    }
  }
  
  #test and convert background if necessary #NOTE: one could also test that all inputForTest items are found in the background, when applicable
  if(!is.null(w$background)){
    if(
      (is.data.frame(w$background) & w$myDataType=="peaks") | #data frame of peaks
      (typeof(w$background)=="S4"  & w$myDataType=="peaks") | #genomicRanges of peaks
      (!is.list(w$background)     & w$myDataType=="genes")   #vector of genes
    ){
      w$backgroundCount="single"
      w$background=list(universalBackground=w$background)
      message("seems like there is provided only one background, which will be used for all tests")
      if(!is.null(w$testSetToBackgroundMapping)) stop("In this case a testSetToBackgroundMapping cannot be provided")
    }else if(
      (is.list(w$background) & w$myDataType=="peaks") | #list of sets of peaks (not data.frame as per above)
      (is.list(w$background) & !is.data.frame(w$background) & w$myDataType=="genes")   #list of sets of genes
    ){ #NOTE: the parsing here could be better
      w$backgroundCount="multiple"
      message("seems like multiple backgrounds were provided")
    }else{
      stop("something is wrong with the provided background")
    }
  }
  
  #test parameters for topTable
  if(w$inputType=="topTable"){
    checkTopTableCols=function(x){all(c(myGRangesCols, "adj.P.Val") %in% colnames(x))}
    if(!all(sapply(w$inputForTest,checkTopTableCols))) stop("one or more of the required cols when the input is from topTable")
    #test if the PeakID to genomic position is unambiguous
    uniquePosAndIdCount=nrow(unique(do.call(rbind,lapply(w$inputForTest,function(x)x[,myGRangesCols]))))
    uniqueIdCount=length(unique(unlist(lapply(w$inputForTest,function(x)x$PeakID))))
    if(uniquePosAndIdCount!=uniqueIdCount) stop("ambiguity in topTable PeakIDs")
    rm(uniquePosAndIdCount)
    #test for duplicates
    if(any(sapply(w$inputForTest,function(x)anyDuplicated(x$PeakID))>0)) stop("duplicated peaks in input for test")
    #handle background
    if(is.null(w$background)){
      #background was not provided
      extractSortedBg=function(x){x[order(x$PeakID),myGRangesCols]}
      if(length(unique(sapply(w$inputForTest,nrow)))==1){ #if all topTable data frames have the same length
        if(nrow(w$inputForTest[[1]])==uniqueIdCount){ #since they all are equally long and have the same number of peaks as the number of unique peaks they must be identical
          w$backgroundCount="single"
          if(w$myDataType=="peaks"){
            w$background=list(universalBackground=extractSortedBg(w$inputForTest[[1]]))
          }else if(w$myDataType=="genes"){
            w$background=list(universalBackground=w$inputForTest[[1]]$PeakID)
          }else{stop("oops")}
          message("seems like all topTable inputs are based on the same background. It will be used for all tests")
        }else{
          w$backgroundCount="multiple"
        }
      }else{
        w$backgroundCount="multiple"
      }
      
      if(w$backgroundCount=="multiple"){ #generating the least number of unique backgrounds
        message("seems like the topTable inputs are based on different backgrounds")
        w$background=list()
        w$testSetToBackgroundMapping=vector()
        for(mySetName in names(w$inputForTest)){
          bgWasFound=F
          for(myBgName in names(w$background)){
            if(identical(w$background[[myBgName]]$PeakID,sort(w$inputForTest[[mySetName]]$PeakID))){
              bgWasFound=T
              w$testSetToBackgroundMapping[mySetName]=myBgName
              break
            }
          }
          if(!bgWasFound){
            newBgName=paste0("bg",length(w$background)+1)
            w$background[[newBgName]]=extractSortedBg(w$inputForTest[[mySetName]])
            w$testSetToBackgroundMapping[mySetName]=newBgName
          }
          rm(bgWasFound,myBgName)
        }
        rm(mySetName)
        if(w$myDataType=="genes"){
          w$background=sapply(w$background,function(x)x$PeakID) #keep only ensemblIDs and not whole table
        }
        message("testSetToBackgroundMapping automatically generated")
      }
    }
    rm(uniqueIdCount)
  }
  
  #should we later combine data across groups (and thus also not plot individual ones)
  w$combineResults=w$combineAcrossMultipleBackgrounds & w$backgroundCount!="single"
  
  
  #convert background(s) to GRanges if necessary, and do some checks
  if(w$myDataType=="peaks"){
    backgroundTypes=unique(sapply(w$background,typeof))
    if(length(backgroundTypes)>1) stop("invalid input for >>background<< or the automatically generated backgrounds doesn't make sense")
    if(any(as.logical(mapply(function(x){anyDuplicated(x$PeakID)},w$background)))) stop("Duplicated PeakID in background")
    if(backgroundTypes=="list"){ #we assume it then is data frames otherwise we assume it is GRanges
      if(any(mapply(function(x){!all(myGRangesCols %in% colnames(x))},w$background))) stop("One or more required column missing from inputted data.frame(s)")
      w$background=mapply(makeGRangesFromDataFrame,w$background,MoreArgs=list(keep.extra.columns=T))
      names(w$background)
    }else{
      if(any(mapply(function(x){!all(names(x)==x$PeakID)},w$background))) stop("Mismatch between genomic ranges names and PeakID collumn")
    }
  }
  
  #make testSetToBackgroundMapping for when there's a single background
  if(w$backgroundCount=="single") {
    w$testSetToBackgroundMapping=setNames(rep("universalBackground",length(w$inputForTest)),names(w$inputForTest))
  } else if(length(w$testSetToBackgroundMapping) == 0) {
    w$testSetToBackgroundMapping = setNames(names(w$background), names(w$background))
  }
  
  #test that we have backgrounds for everything
  if(any(!names(w$inputForTest) %in% names(w$testSetToBackgroundMapping))) stop("one or more test input not mapped to a background")
  
  
  #map background peaks to genes by the aforementioned gene regulatory domains
  if(w$myDataType=="peaks"){
    w$myGeneRegDoms=do.call(geneRegDomains,c(furtherArgs__geneRegDomains,list(alternativeGeneInfo=w$myGenesTable, GENOME_VERSION=w$GENOME_VERSION)))
    message("generated gene regulatory domains genes")
  }
  
  
  #map background(s) to peaks
  if(w$myDataType=="peaks"){
    w$myBackgroundMapping=mapply(geneToPeakMapping, w$background, MoreArgs=list(regDom=w$myGeneRegDoms),SIMPLIFY=F)
    message("mapped peaks in background(s) to genes by the >>gene regulatory domains<<")
  }
  
  
  #filter gene backgrounds to master gene list
  if(w$myDataType=="genes"){
    w$background=sapply(w$background,intersect,w$myGenes,simplify=F)
    message("filtered backgrounds to genes in >>myGenesFallback<<")
  }
  
  
  w$background#get gene meta sets
  if(is.null(w$myGeneMetaSets)){
    message("grabbing geneMetaSets, which may take a while")
    w$myGeneMetaSets=standardFisherGeneSets(GENOME_VERSION=w$GENOME_VERSION,greatGenes=w$myGenesTable)$standardGeneSets
    message("finished grabbing geneMetaSets")
  }else{
    message("used provided geneMetaSets")
  }
  
  
  #map backgrounds to gene sets if test set is peaks. This is done once for each background times each geneMetaSet
  if(w$myDataType=="peaks"){
    message("Doing geneset to peak mapping. The next couple of messages will be jumbled if using parallel processing")
    geneSetToPeakSetWrapper=function(myPeaksGeneToPeakMap){
      cleanNull(mcmapply(
        function(x)do.call(geneSetToPeakSet,c(furtherArgs__geneSetToPeakSet,list(geneToPeakMap=myPeaksGeneToPeakMap, myGeneSets=x)))
        , w$myGeneMetaSets,SIMPLIFY=F))
    }
    w$myBackgroundPeakSets=mapply(geneSetToPeakSetWrapper,w$myBackgroundMapping,SIMPLIFY=F)
    endFlag = sapply(w$myBackgroundPeakSets,function(x){(length(x)==0)})
    if(all(endFlag)) {
      return(list())
    }
    message("mapped gene meta sets to peaks in background(s) through the gene regulatory domains")
  }
  
  #filter gene sets to backgrounds if test set is genes. This is done once for each background times each geneMetaSet
  if(w$myDataType=="genes"){
    minGenesInGeneSet=5 #NOTE: hardcoding
    filterSetToBackground=function(x,bg){
      z=list()
      z$sets=sapply(x$sets,intersect,bg,simplify=F)
      z$sets=z$sets[sapply(z$sets,length)>=minGenesInGeneSet]
      z$metadata=x$metadata[x$metadata$name %in% names(z$sets),]
      if(length(z$sets)==0) z=NULL
      return(z)
    }
    w$filteredGeneMetaSets=mcmapply(function(bg)
      cleanNull(sapply(w$myGeneMetaSets,filterSetToBackground,bg,simplify=F))
      ,w$background,SIMPLIFY=F)
    endFlag = sapply(w$filteredGeneMetaSets,function(x){(length(x)==0)})
    if(all(endFlag)) {
      return(list())
    }
    message("filtered gene meta sets to background sets")
  }
  
  #################################
  #do fisher gseas:
  
  #gsea function
  fisherWrapper=function(testGenes,geneMetaSets,myGenes,outDir){
    fisherHelper=function(metaSet){
      message(metaSet)
      gom.obj=newGOM(testGenes, geneMetaSets[[metaSet]]$sets, length(myGenes))
      #GOM output only has two cols if there's only one group
      if(length(names(testGenes))==1){
        z=cbind(
          names(testGenes),
          sub(paste0("\\.",names(testGenes)),"",rownames(melt(getMatrix(gom.obj, name="pval")))),
          melt(getMatrix(gom.obj, name="pval")),
          melt(getMatrix(gom.obj, name="odds.ratio"))[,1],
          melt(getMatrix(gom.obj, name="intersection"))[,1],
          melt(getMatrix(gom.obj, name="union"))[,1],
          melt(getMatrix(gom.obj, name="Jaccard"))[,1]
        )
      }else{
        z=cbind(
          melt(getMatrix(gom.obj, name="pval")),
          melt(getMatrix(gom.obj, name="odds.ratio"))[,3],
          melt(getMatrix(gom.obj, name="intersection"))[,3],
          melt(getMatrix(gom.obj, name="union"))[,3],
          melt(getMatrix(gom.obj, name="Jaccard"))[,3]
        )
      }
      colnames(z) = c("Set", "Reference", "pval", "odds.ratio", "intersection", "union", "Jaccard")
      addGseaInfoAndPlot(z,geneMetaSets,metaSet,outDir,forceNoPlot=w$combineResults,forceNoOutDirCheck=forceNoOutDirCheck,setMetadata=w$inputForTestMetadata)
    }
    mcmapply(fisherHelper,names(geneMetaSets),SIMPLIFY=F)
  }
  
  #code using gsea function
  if(w$testMethod=="fisher"){
    w$inputForTestComplete=w$inputForTest
    if(w$inputType=="topTable"){
      w$inputForTest=lapply(w$inputForTest,function(x)x$PeakID[x$adj.P.Val<=qValueForFisher]) #grab PeakIDs of significant peaks/genes
    }
    if(w$myDataType=="peaks"){ #discard peaks not mapped to genes
      w$inputForTest=sapply(names(w$inputForTest),function(x)intersect(w$inputForTest[[x]],w$myBackgroundMapping[[w$testSetToBackgroundMapping[x]]]$peaksWithGeneAssignment),simplify=F)
      message("filtered inputForTest to only peaks with one or more assigned genes in the background")
    }else if(w$myDataType=="genes"){ #discard genes not in background
      w$inputForTest=sapply(names(w$inputForTest),function(x)intersect(w$inputForTest[[x]],w$background[[w$testSetToBackgroundMapping[x]]]),simplify=F)
      message("filtered inputForTest to genes in the background(s)")
    }else{stop("oops")}
    nonEmptyInputSets=unname(sapply(w$inputForTest,length)>0)
    if(any(!nonEmptyInputSets)) message("discarded one or more of the inputForTest, because there were no genes/peaks")
    w$inputForTest=w$inputForTest[nonEmptyInputSets] #discard empty test sets
    rm(nonEmptyInputSets)
    w$testSetToBackgroundMapping=w$testSetToBackgroundMapping[names(w$testSetToBackgroundMapping) %in% names(w$inputForTest)] #discard unused testSetToBackgroundMapping
    w$background=w$background[unique(w$testSetToBackgroundMapping)] #discard unused backgrounds
    if(w$myDataType=="peaks"){
      w$myBackgroundPeakSets=w$myBackgroundPeakSets[unique(w$testSetToBackgroundMapping)]#discard unused myBackgroundPeakSets
      w$myBackgroundMapping=w$myBackgroundMapping[unique(w$testSetToBackgroundMapping)] #discard unused background mappings
    }
    
    #this is done seperately for each background:
    w$gseaResults=mapply(
      function(x){
        outSubDir=if(is.null(outDir)){NULL}else{paste0(outDir,"/",x)}
        #the variable name "testGenes" is either the test genes or test peaks depending on what type the input data is.
        if(w$myDataType=="peaks"){
          fisherWrapper(
            testGenes=w$inputForTest[names(w$testSetToBackgroundMapping)[w$testSetToBackgroundMapping==x]],  
            geneMetaSets=w$myBackgroundPeakSets[[x]],
            myGenes=w$myBackgroundMapping[[x]]$peaksWithGeneAssignment,
            outDir=outSubDir
          )
        }else if(w$myDataType=="genes"){
          fisherWrapper(
            testGenes=w$inputForTest[names(w$testSetToBackgroundMapping)[w$testSetToBackgroundMapping==x]],  
            geneMetaSets=w$filteredGeneMetaSets[[x]],
            myGenes=w$background[[x]],
            outDir=outSubDir
          )
        }else{stop("oops")}
      }
      ,names(w$background),SIMPLIFY=F)
  }
  
  #################################
  #do cameraPR gseas:
  
  #gsea function
  cameraPRWrapper=function(topTables,geneMetaSets,outDir){ #for every background
    cameraPRHelper=function(metaSet){ #for every geneMetaSet
      cameraPRSubHelper=function(myInputName){ #for every topTable
        myInput=topTables[[myInputName]]
        if(!w$useRankingForCameraPR){
          if(!("t" %in% colnames(myInput))) stop("required collumn >>t<< in topTable missing but required for cameraPR when useRankingForCameraPR==FALSE. Maybe the resuls are from anova? Consider rerunning with useRankingForCameraPR==TRUE")
          z=cameraPR(setNames(myInput$t,myInput$PeakID), geneMetaSets[[metaSet]]$sets)
        }else if(w$useRankingForCameraPR){
          if(sum(c("F","t") %in% colnames(myInput))!=1) stop("for useRankingForCameraPR==TRUE, topTables must contain either >>F<< or >>t<< collumn")
          if("t" %in% colnames(myInput)){
            if(w$useAbsTinRankingForCameraPR){
              z=cameraPR(setNames(abs(myInput$t),myInput$PeakID), geneMetaSets[[metaSet]]$sets,use.ranks=T)
            }else{
              z=cameraPR(setNames(myInput$t,myInput$PeakID), geneMetaSets[[metaSet]]$sets,use.ranks=T)
            }
          }else{
            z=cameraPR(setNames(myInput$F,myInput$PeakID), geneMetaSets[[metaSet]]$sets,use.ranks=T)
          }
        }else{stop("unexpected useRankingForCameraPR parameter")}
        z=data.frame(
          Set=myInputName,
          Reference=rownames(z),
          pval=z$PValue,
          odds.ratio=NA,
          intersection=NA,
          union=NA,
          Jaccard=NA,
          NGenes=z$NGenes,
          Direction=z$Direction,
          stringsAsFactors=F
        )
        z$odds.ratio[z$Direction=="Up"]=2
        z$odds.ratio[z$Direction=="Down"]=0.5
        z
      }
      message(metaSet)
      myGseas=do.call(rbind,sapply(names(topTables),cameraPRSubHelper,simplify=F))
      addGseaInfoAndPlot(myGseas,geneMetaSets,metaSet,outDir,forceNoPlot=w$combineResults,forceNoOutDirCheck=forceNoOutDirCheck,setMetadata=w$inputForTestMetadata)
    }
    mcmapply(cameraPRHelper,names(geneMetaSets),SIMPLIFY=F)
  }
  
  #code using gsea function
  if(w$testMethod=="cameraPR"){
    #filter topTable rows to only relevant peaks/genes
    w$inputForTestComplete=w$inputForTest
    if(w$myDataType=="peaks"){ #discard peaks not mapped to genes
      w$inputForTest=sapply(names(w$inputForTest),function(x){myTopTable=w$inputForTest[[x]];myTopTable[ myTopTable$PeakID %in% w$myBackgroundMapping[[w$testSetToBackgroundMapping[x]]]$peaksWithGeneAssignment ,]},simplify=F)
      message("filtered inputForTest to only peaks with one or more assigned genes in the background")
    }else if(w$myDataType=="genes"){ #discard genes not in background
      w$inputForTest=sapply(names(w$inputForTest),function(x){myTopTable=w$inputForTest[[x]];myTopTable[ myTopTable$PeakID %in% w$background[[w$testSetToBackgroundMapping[x]]] ,]},simplify=F)
      message("filtered inputForTest to genes in the background(s)")
    }else{stop("oops")}
    
    #this is done seperately for each background:
    w$gseaResults=mapply(
      function(x){
        outSubDir=if(is.null(outDir)){NULL}else{paste0(outDir,"/",x)}
        if(w$myDataType=="peaks"){
          cameraPRWrapper(
            topTables=w$inputForTest[names(w$testSetToBackgroundMapping)[w$testSetToBackgroundMapping==x]],  
            geneMetaSets=w$myBackgroundPeakSets[[x]],
            outDir=outSubDir
          )
        }else if(w$myDataType=="genes"){
          cameraPRWrapper(
            topTables=w$inputForTest[names(w$testSetToBackgroundMapping)[w$testSetToBackgroundMapping==x]],  
            geneMetaSets=w$filteredGeneMetaSets[[x]],
            outDir=outSubDir
          )
        }else{stop("oops")}
      }
      ,names(w$background),SIMPLIFY=F)
  }
  
  #if we have data from multipleBackgrounds and want to combine:
  if(w$combineResults){
    message("Data from multiple backgrounds present. Aggregating these. Clustering in plots are here based on the original gene sets")
    w$individualGseaResults=w$gseaResults
    z=sapply(names(w$myGeneMetaSets),function(x)do.call(rbind.fill,lapply(names(w$gseaResults),function(y)data.frame(w$gseaResults[[y]][[x]],bg=y,stringsAsFactors=F))),simplify=F)
    outSubDir=if(is.null(outDir)){NULL}else{paste0(outDir,"/aggregatedBackground")}
    w$gseaResults=sapply(names(z),function(x)addGseaInfoAndPlot(z[[x]],w$myGeneMetaSets,x,outDir=outSubDir,addGseaInfo=F,forceNoOutDirCheck=forceNoOutDirCheck,setMetadata=w$inputForTestMetadata),simplify=F)
    rm(z,outSubDir)
  }
  
  if(shrinkOutput){
    w$myGeneMetaSets=NULL
    w$inputForTest=NULL
    w$myBackgroundPeakSets=NULL
  }
  return(w)
}


##################################################################
# A traditional overlap GSEA with Fisher's exact test
##################################################################

fisherGsea=function(
    testGenes, #either just a vector of genes or a list of vectors of genes to be tested independently. Use Ensembl IDs
    geneMetaSets=NA, #sets of gene sets we should test. If none provided, we use a standard repetoire. You should filter them to the provided genes you want to test against
    myGenes=NA, #List of all genes to be considered in analysis. We filter analysis down to these genes. If NA we use the GREAT genes
    gseaExcludeChr=c("chrM","chrY"), #only applied if myGenes isnt provided
    myGenesFallback=NA,
    outDir=NULL,
    forceNoOutDirCheck=F,
    cores=NULL, #for parallel if NULL then don't change core settings processing (only partially implemented)
    GENOME_VERSION="hg19"
){ #FIXME:
  if(!is.null(cores)){
    message(paste("Note: setting mc.cores to",cores))
    options(mc.cores=cores)
  }
  
  if(is.na(myGenesFallback)) { myGenesFallback = ifelse(GENOME_VERSION == "hg19", myGenesFallback_hg19, ifelse(GENOME_VERSION == "hg38", myGenesFallback_hg38, myGenesFallback_mm10)) }
  
  if(!is.list(testGenes)) #if just a vector, convert it to a list for easier processesing
    testGenes=list(myTestGenes=testGenes)
  
  
  
  #if no get of complete genes to test against
  if(is.na(myGenes)[1]){
    myGenes=read.delim(myGenesFallback,stringsAsFactors=F,col.names=c("chr","tss","strand","ensembl"))
    myGenes=myGenes[!myGenes$chr %in% gseaExcludeChr,] #since we are running GREAT with a subset of genes (their default) it makes sense go get rid of the other ones up front by using this as a filter
    myGenes=myGenes$ensembl
  }
  
  if(is.na(geneMetaSets)[1]){
    geneMetaSets=c(
      brainSets=grabBrainGeneSets(genesToInclude=myGenes, GENOME_VERSION=GENOME_VERSION),
      grabMsigdbGeneSets2022_1_PLUS_SYNGO(customName="msigdbSets",genesToInclude=myGenes, GENOME_VERSION=GENOME_VERSION),
      #msigdbSets=list(combineGeneSetObjects(grabMsigdbGeneSets6_0(genesToInclude=myGenes, GENOME_VERSION=GENOME_VERSION),includeGroupNameInNewName=F, GENOME_VERSION=GENOME_VERSION)), #we just put cannonical pathways and GO:BP together
      ts7=list(grabTargetScan7(genesToInclude=myGenes))
    )
  }
  
  #Filter provided gene
  testGenes=sapply(names(testGenes),function(x)testGenes[[x]][testGenes[[x]] %in% myGenes],simplify=F)
  
  fisherSubFunction=function(metaSet){
    gom.obj=newGOM(testGenes, geneMetaSets[[metaSet]]$sets, length(myGenes))
    #GOM output only has two cols if there's only one group
    if(length(names(testGenes))==1){
      z = cbind(
        names(testGenes),
        sub(paste0("\\.",names(testGenes)),"",rownames(melt(getMatrix(gom.obj, name="pval")))),
        melt(getMatrix(gom.obj, name="pval")),
        melt(getMatrix(gom.obj, name="odds.ratio"))[,1],
        melt(getMatrix(gom.obj, name="intersection"))[,1],
        melt(getMatrix(gom.obj, name="union"))[,1],
        melt(getMatrix(gom.obj, name="Jaccard"))[,1]
      )
    }else{
      z = cbind(
        melt(getMatrix(gom.obj, name="pval")),
        melt(getMatrix(gom.obj, name="odds.ratio"))[,3],
        melt(getMatrix(gom.obj, name="intersection"))[,3],
        melt(getMatrix(gom.obj, name="union"))[,3],
        melt(getMatrix(gom.obj, name="Jaccard"))[,3]
      )
    }
    
    colnames(z) = c("Set", "Reference", "pval", "odds.ratio", "intersection", "union", "Jaccard")
    z$Bonf_AdjP = p.adjust(z$pval, method = "bonferroni")
    z$BH_AdjP = p.adjust(z$pval, method = "BH")
    z$Z = qnorm(1-(z$pval)/2)
    z$LogP = -log10(z$pval)
    z$FDR_AdjP=fdrtool(z$pval, statistic="pvalue", plot=F,cutoff.method="fndr", verbose=F)$qval
    z=z[order(z$pval),]
    #z$SetSize=sapply(z$Set, function(x) length(testGenes[[x]]))
    #z$ReferenceSize=sapply(z$Reference, function(x) length(geneMetaSets[[metaSet]]$sets[[x]]))
    
    z=cbind(z, geneMetaSets[[metaSet]]$metadata[match(z$Reference,geneMetaSets[[metaSet]]$metadata$name),])
    z$name=NULL
    
    z[order(z$pval),]
  }
  myGseas=mcmapply(fisherSubFunction,names(geneMetaSets),SIMPLIFY=F)
  
  if(!is.null(outDir)){
    #save analysis
    gseaPlotter(myGseas,geneMetaSets,outDir,forceNoOutDirCheck)
    save(list=ls(all=T), file=paste0(outDir,"/gsea.Rdata"), envir=environment())
    
  }
  
  return(myGseas)
}


##################################################################
# plot gsea results
##################################################################
#function to plot gsea results such as those from fisherGSEA and universalGsea

gseaPlotter=function(
    myGseas,
    geneMetaSets, #required for clustering
    outDir,
    doBiclust=F, #sometimes fails, especially if there are NAs so now it's disabled by default
    forceNoOutDirCheck=F,
    doTopSomethingPlots=T,
    customPlotArgs=NULL,
    plotScale=1 #only applies to ggplot
){
  #Check dir
  if(!forceNoOutDirCheck & file.exists(outDir)) stop("the provided output dir already exists. please provide a non-existant dir. Script aborted.")
  dir.create(outDir,recursive=T,showWarnings=F)
  
  #write output datables
  tsvDir=paste0(outDir,"/result_textFiles")
  dir.create(tsvDir,showWarnings=F)
  tsvIndiDir=paste0(outDir,"/result_textFiles_individual")
  dir.create(tsvIndiDir,showWarnings=F)
  niceColumns=c("Set", "pval", "odds.ratio", "intersection", "union", "Jaccard", "Bonf_AdjP", "BH_AdjP", "name_full","group")
  lapply(names(myGseas),function(metaSet){
    mtsv(myGseas[[metaSet]][,niceColumns],tsvDir,gz=T,fileBaseName=make.names(metaSet))
    mcmapply(function(mySet) mtsv(myGseas[[metaSet]][myGseas[[metaSet]]$Set==mySet,niceColumns],tsvIndiDir,gz=T,fileBaseName=make.names(paste0(metaSet,"__",mySet))),unique(myGseas[[metaSet]]$Set))
  })
  
  #create a copy of the data with some changes in variable names for easy reuse of greatR plot code. I admit it is dirty coding
  z=list()
  z$aggGsea=sapply(myGseas,function(x){
    x$regDomEnrichment=x$odds.ratio
    if(is.null(x$SetFullName)){
      x$peaksFullName=x$Set #peak set in previous implementation a list of genes to be analyzed here
    }else{
      x$peaksFullName=x$SetFullName #peak set in previous implementation a list of genes to be analyzed here
    }
    x$name=x$Reference #file name friendly version of gene set name
    x$binomTest=x$pval #binomTest was the name of the pval col in greatR approach
    x$plotLabel=""
    x$plotLabel[x$pval<0.05]="·" #NOTE: different labling than on tradional greatR figs
    x$plotLabel[x$BH_AdjP<0.05]="#"
    x$plotLabel[is.na(x$BH_AdjP)]="NA"
    x
  },simplify=F)
  
  #heatmaps: plot all pathways
  heatmapDir=paste0(outDir,"/heatmaps")
  dir.create(heatmapDir,showWarnings=F,recursive=T)
  mcmapply(function(i){
    if(length(unique(z$aggGsea[[i]]$name))<101){ #we cannot plot like 4k pathways
      mySetName=z$aggGsea[[i]]$geneMetaSets[1]
      subDir=paste0(heatmapDir,"/",mySetName)
      rm(mySetName)
      dir.create(subDir,showWarnings=F,recursive=T)
      annoName=z$aggGsea[[i]]$customPlotName[1]
      if(is.null(annoName)) annoName=i
      aggGseaPlotter(df=z$aggGsea[[i]],annoName=annoName,outDir=subDir,customPlotArgs=customPlotArgs,plotScale=plotScale)
      if(!is.null(geneMetaSets))
        aggGseaPlotter(df=z$aggGsea[[i]],annoName=annoName,outDir=subDir,doCluster=T,geneSets=geneMetaSets[[i]]$sets,customPlotArgs=customPlotArgs,plotScale=plotScale)
    }
  },names(z$aggGsea))
  
  #heatmaps: Top something best pathways from each peakSet
  if(doTopSomethingPlots){
    md=expand.grid(metaSet=names(z$aggGsea),topCount=c(3,5,10,15,20,25),stringsAsFactors=F)
    md$enoughSets=mapply(function(x,y){length(unique(z$aggGsea[[x]]$name))>=y},md$metaSet,md$topCount) #enough sets to plot
    md$myName=paste0(md$metaSet,"_top",md$topCount)
    md=md[md$enoughSets,]
    aggTopPlotter=function(myName,metaSet,topCount){
      w=z$aggGsea[[metaSet]]
      w=w[order(w$binomTest),]
      myPathways=unique(unlist(lapply(unique(w$peaksFullName),function(x){w$name_full[w$peaksFullName==x][1:topCount]})))
      w=w[w$name_full %in% myPathways,]
      mySetName=w$geneMetaSets[1]
      subDir=paste0(heatmapDir,"/",mySetName)
      rm(mySetName)
      dir.create(subDir,showWarnings=F,recursive=T)
      aggGseaPlotter(df=w, annoName=myName, outDir=subDir,customPlotArgs=customPlotArgs,plotScale=plotScale)
      if(!is.null(geneMetaSets))
        aggGseaPlotter(df=w, annoName=myName, outDir=subDir, doCluster=T,  geneSets=geneMetaSets[[metaSet]]$sets,customPlotArgs=customPlotArgs,plotScale=plotScale)
      if(doBiclust & length(unique(w$peaksFullName))>1){ #we can only bicluster with 2+ samples
        plotGseaBiclust(df=w, annoName=myName, outDir=subDir)
      }
      return(w)
    }
    z$topAggGsea=with(md,mcmapply(aggTopPlotter,myName,metaSet,topCount,SIMPLIFY=F))
  }
  
  #barplots: top5
  barplotDir=paste0(outDir,"/barplots")
  lapply(names(z$aggGsea),function(metaSet)singleGseaPlotter(df=z$aggGsea[[metaSet]],annoName=metaSet,barplotDir))
  return(invisible(z))
}


#################################
aggGseaPlotter=function(df,annoName,outDir,doCluster=F,geneSets=NA,customClustMethod=NA,showMissingFields=T,customPlotArgs=NULL,plotScale=1){
  message(paste0("plotting ", annoName))
  
  #init
  myPalette = colorRampPalette(brewer.pal(9, "Greens"), space="Lab")
  myPalette2way=colorRampPalette(rev(c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#F7F7F7","#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")), space="Lab")
  
  #NOTE:ugly hardcoded fix of a typo:
  if(is.ordered(df$name_full)){
    df$name_full=ordered(
      sub("^Endothelian$","Endothelial",df$name_full),
      levels=sub("^Endothelian$","Endothelial",levels(df$name_full))
    )
  }else{
    df$name_full=sub("^Endothelian$","Endothelial",df$name_full)
  }
  
  #add log2 of enrichment
  df$log2regDomEnrichment=log2(df$regDomEnrichment)
  
  #optionally cluster
  if(doCluster){
    if(!is.list(geneSets)) stop("if you wan't to cluster, you must provide the original gene sets as a list")
    df$name_full=as.character(df$name_full)
    myClustering=geneSetClustering(metadata=df,sets=geneSets,plotName=annoName,customClustMethod=customClustMethod)
    if(is.na(customClustMethod)){
      df$name_full=ordered(df$name_full,levels=df$name_full[match(myClustering$orderedGeneSets_Ward.D2,df$name)])
    }else{
      df$name_full=ordered(df$name_full,levels=df$name_full[match(myClustering$orderedGeneSets_custom,df$name)])
    }
  }
  
  #plots
  
  #proper symmetric scaling for log2regDomEnrichment plots
  if(any(is.infinite(df$log2regDomEnrichment)) | all(is.na(df$log2regDomEnrichment))){
    myLims=c(-1,1)
  }else{
    myLims=max(abs(df$log2regDomEnrichment[!is.infinite(df$log2regDomEnrichment)]),na.rm=T)
    if(myLims==0) myLims=1
    myLims=c(-1,1)*myLims
  }
  
  if(showMissingFields){
    allCombs=expand.grid(name_full=unique(df$name_full),peaksFullName=unique(df$peaksFullName)) #all possible name_full and peakFullName combos
    z=rbind(df[,c("name_full","peaksFullName")],allCombs) #combine allCombs with those with those actually found 
    z=z[seq(nrow(z))>nrow(df) & !duplicated(z),] #take only those from that joined table from allCombs, which aren't found in df
    if(nrow(z)){
      z$plotLabel="NA"
      df=rbind.fill(df,z) #add the missing ones to the end. all the other cols will just be NAs
    }
    rm(z,allCombs)
  }
  
  for(flipped in c(F,T)){
    if(flipped){
      df$xCol=df$peaksFullName
      df$yCol=df$name_full
    }else{
      df$xCol=df$name_full
      df$yCol=df$peaksFullName
    }
    
    zz=ggplot(df,aes(xCol,yCol)) +
      ggtitle(annoName) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      ggtitle(annoName) +
      theme_classic() +
      theme(axis.text=element_text(colour="black")) +
      coord_fixed() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(axis.title = element_blank()) +
      theme(plot.title=element_text(face="bold")) +
      theme(legend.title = element_text(size = 12, face = "bold")) + 
      customPlotArgs
    
    zz_p = zz + geom_tile(aes(fill=LogP)) + scale_fill_gradientn(colours = myPalette(100),name="LongSpacer\nMinus\nlog P")
    zz_psq = zz + geom_tile(aes(fill=LogP)) + scale_fill_gradientn(colours = myPalette(100),trans="sqrt",name="LongSpacer\nMinus\nlog P")
    zz_lfc = zz + geom_tile(aes(fill=log2regDomEnrichment)) + scale_fill_gradientn(colours = myPalette2way(100),limits=myLims,name="LongSpacer\nLog2 fold\nchange")
    zz_fc = zz + geom_tile(aes(fill=regDomEnrichment)) + scale_fill_gradientn(colours = myPalette(100),name="LongSpacer\nFold\nchange")
    
    if(doCluster){
      if(is.na(customClustMethod)){
        cl=ggdendrogram(myClustering$clustering_Ward.D2,labels=F,rotate=flipped)+theme_dendro()
      }else{
        cl=ggdendrogram(myClustering$clustering_custom ,labels=F,rotate=flipped)+theme_dendro()
      }
    }else{
      cl=zz_p+geom_blank() + theme_void() + theme(legend.position="none",plot.title=element_blank()) #newer version of grid.arrange/ggplots must have something like zz_p as an argument. This plot is just a filler and isn't used.
    }
    
    filename=paste0(outDir,"/plot_gsea_",gsub(" ","_",annoName),"_clustered_",doCluster,ifelse(flipped,"_flipped.pdf",".pdf"))
    geneSetWidth=5+0.2*length(unique(df$name_full)) 
    peakWidth=3*(5+0.2*length(unique(df$peaksFullName)))
    myHeight=ifelse(flipped,geneSetWidth,peakWidth)
    myWidth=ifelse(flipped,peakWidth,geneSetWidth) + 3*1.7 # to make room for legend
    
    myRow=ifelse(flipped,1,3)
    myCol=ifelse(flipped,3,1)
    
    pdf(filename,height=myHeight*plotScale,width=myWidth*plotScale); #first the two most interesting plots. Then everything systematically
    grid.arrange(zz_psq, zz_lfc+geom_text(aes(label=plotLabel),alpha=0.7),cl,nrow=myRow,ncol=myCol)
    grid.arrange(zz_psq, zz_psq+geom_text(aes(label=plotLabel),alpha=0.7),cl,nrow=myRow,ncol=myCol)
    grid.arrange(zz_p, zz_p+geom_text(aes(label=plotLabel),alpha=0.7),cl,nrow=myRow,ncol=myCol)
    grid.arrange(zz_lfc, zz_lfc+geom_text(aes(label=plotLabel),alpha=0.7),cl,nrow=myRow,ncol=myCol)
    grid.arrange(zz_fc, zz_fc+geom_text(aes(label=plotLabel),alpha=0.7),cl,nrow=myRow,ncol=myCol)
    dev.off()
  }
}


plotGseaBiclust=function(df,outDir,annoName){
  library(reshape2)
  library(gplots)
  df$LogP[is.infinite(df$LogP)]=NA
  myLogPval=acast(df,name_full~peaksFullName,value.var="LogP")
  myPlotLabel=acast(df,name_full~peaksFullName,value.var="plotLabel")
  filename=paste0(outDir,"/plot_gsea_heatmap_biclust_",annoName,".pdf")
  #heatmap.2(...,distfun = function(x) dist(x,method = 'euclidean'),...)
  
  sink("/dev/null") #get rid of extensive printing to output
  pdf(filename,height=11+0.15*nrow(myLogPval),width=7+0.5*ncol(myLogPval))
  print(heatmap.2(myLogPval,margins=c(30,30),notecol="black",trace="none",keysize=1.0,cellnote=myPlotLabel))
  print(heatmap.2(myLogPval,margins=c(30,30),notecol="black",trace="none",keysize=1.0))
  print(heatmap.2(sqrt(myLogPval),margins=c(30,30),notecol="black",trace="none",keysize=1.0,cellnote=myPlotLabel))
  print(heatmap.2(sqrt(myLogPval),margins=c(30,30),notecol="black",trace="none",keysize=1.0))
  dev.off()
  sink()
}


singleGseaPlotter=function(df,annoName,outDir,topResults=5){
  library(ggplot2)
  dir.create(outDir,recursive=T,showWarnings=F)
  pdf(paste0(outDir,"/",make.names(annoName),".pdf"),height=1.5+0.2*topResults,width=9)
  for(mySet in unique(df$peaksFullName)){
    subDf=df[df$peaksFullName==mySet,]
    subDf=subDf[order(-subDf$LogP),]
    subDf=subDf[1:min(nrow(subDf),topResults),]
    subDf$name_full=ordered(subDf$name_full,levels=rev(unique(subDf$name_full)))
    print(
      ggplot(subDf,aes(name_full,LogP)) +
        geom_bar(stat="identity") +
        theme_classic() +
        theme(axis.text=element_text(colour="black")) +
        ggtitle(mySet) +
        xlab("Gene Set") +
        ylab("-logP") +
        scale_y_continuous(expand = c(0, 0)) +
        coord_flip()
    )
  }
  dev.off()
}
##################################################################

