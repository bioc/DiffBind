pv.contrast <- function(pv,group1,group2=!group1,name1,name2,
                        minMembers=3,categories,bMulti=FALSE,block, bNot=FALSE,
                        design=TRUE, contrast, bGetNames=FALSE, reorderMeta) {
   
   numStart <- length(pv$contrasts)
   
   pv <- pv.orderMetaFacs(pv, reorderMeta)
   
   if(bGetNames) {
      if(!is.logical(design)) {
         pv$design <- design
         pv$DESeq2 <- pv$edgeR <- NULL
      }
      return(pv.DESeq2ResultsNames(pv))
   }
   
   if(!missing(reorderMeta) && missing(group1) && 
      missing(contrast) && !is.null(pv$contrasts) &&
      is(design,"logical") ) {
      pv <- pv.CheckReleveledContrasts(pv)
      return(pv)
   }
   
   if(missing(block)) {
      doBlock <- FALSE
   } else {
      doBlock <- TRUE
   }
   
   generateDesign <- FALSE
   if(design!=FALSE) {
      if(design[1] == TRUE) {
         if (!missing(group1)) {
            stop("Can not specify groups when design is TRUE, specify contrast.",call.=FALSE)
         } else {
            if (bMulti) {
               warning("bComplex ignored when design is TRUE.",call.=FALSE)
               bMulti <- FALSE
            }
            if(bNot) {
               warning("bNot ignored when design is TRUE",call.=FALSE)
               bNot <- FALSE
            }
            generateDesign <- TRUE
         }
      }
      
      if(!missing(block)) {
         warning("block ignored when design is used.",call.=FALSE)
         doBlock <- FALSE
      }
      
      if(!missing(group1)) {
         warning("groups ignored when design is used.",call.=FALSE)
      }
      
      if(!generateDesign || !missing(contrast)) {
         res <- pv.contrastDesign(pv=pv, design=design, contrast=contrast, 
                                  name1=name1, name2=name2, 
                                  bGetNames=bGetNames)
         if(!missing(reorderMeta)) {
            res <- pv.CheckReleveledContrasts(res)
         }
         return(res)
      }
   }
   
   if(missing(name1)) {
      name1="Group1"
   }
   if(missing(name2)) {
      name2="Group2"
   }
   
   if(missing(group1)) { # Automatically Generate Contrasts}
      
      if( (sum(pv.mask(pv,PV_CALLER,'source'))==0) &
          (sum(pv.mask(pv,PV_CALLER,'counts'))==0) ) {
         warning('Model must include count data for contrasts.',call.=FALSE)
         return(pv)
      }
      
      if(pv.isDesign(pv$design)) {
         varnames <- all.vars(formula(pv$design))
         categories <- pv.FactorVals[match(varnames,pv.FactorNames)]
      }
      
      if(!is.null(pv$contrasts)) {
         message("Clearing existing contrasts.")
         pv$contrasts <- NULL # clear existing contrasts
      }
      
      olddesign <- pv$design
      oldDESeq2 <- pv$DESeq2
      pv$design <-  pv$DESeq2 <- pv$edgeR <- NULL
      
      if(missing(categories)) {
         attributes <- c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT)
      } else {
         attributes <- categories
      }
      if(!doBlock) {
         res <- pv.getContrasts(pv,minMembers=minMembers,attributes=attributes,
                                bNot=bNot)
         if(bMulti) {
            res <- pv.listaddto(res,pv.contrastPairs(pv,minMembers=minMembers,
                                                     attributes=attributes,
                                                     conlist=res,bNot=bNot))  
         }
      } else {
         res <- pv.getContrasts(pv,minMembers=minMembers,attributes=attributes,
                                block=block,bNot=bNot)
         if(bMulti) {
            res <- pv.listaddto(res,pv.contrastPairs(pv,minMembers=minMembers,
                                                     attributes=attributes,
                                                     block=block,bNot=bNot))   
         }
      }
      if(!is.null(res)) {
         res <- pv.contrastDups(res)
      } else {
         warning("No contrasts added. There must be at least two sample groups with at least three replicates.",
                 call.=FALSE)	
      }
      if(pv.isDesign(olddesign) || generateDesign) {
         res <- pv.checkContrastOrder(res, pv$meta)
      }
      if(doBlock){
         problem <- NULL
         issues <- 0
         for(i in 1:length(res)){
            if(!pv.checkBlock(res[[i]],bWarning=FALSE)){
               #warning("Blocking factor has unmatched sample(s).")
               issues <- issues+1
               problem <- res[[i]]
               res[[i]]$blocklist <- NULL   		
            }
         }
         if(issues > 0) {
            if(issues < length(res)) {
               warning('Blocking factor not used for some contrasts:',call.=FALSE)	
            } else {
               warning('Blocking factor invalid for all contrasts:',call.=FALSE)
            }
            x <- pv.checkBlock(problem)		
         }      	
      }
      
      if(generateDesign && !is.null(res)) {
         pv$contrasts <- res
         pv <- pv.generateDesignFromContrasts(pv, attributes=attributes,
                                              bGetNames=FALSE)
         if(is.null(olddesign)){
            pv$edgeR  <- NULL
            oldDESeq2 <- NULL
         } else {
            if(olddesign != pv$design) {
               pv$design <- pv.compareDesigns(olddesign, pv$design)
               if(olddesign != pv$design) {
                  pv$edgeR  <- NULL
                  oldDESeq2 <- NULL
               }
            }
         }
         pv$DESeq2 <- oldDESeq2
         pv$DESeq2$names <- pv.DESeq2ResultsNames(pv)
         if(!missing(reorderMeta)) {
            pv <- pv.CheckReleveledContrasts(pv)
         }
         if(bGetNames) {
            return(pv$DESeq2$names)
         } else {
            return(pv)
         }
      }
      
   } else {
      res <- pv.addContrast(pv,group1,group2,name1,name2)
      if(doBlock) {
         res$contrasts[[length(res$contrasts)]]$blocklist <- pv.BlockList(pv,block)
         if(!pv.checkBlock(res$contrasts[[length(res$contrasts)]])) {
            #warning("Blocking factor has unmatched sample(s).")
            res$contrasts[[length(res$contrasts)]]$blocklist <- NULL   	
         }
      }
      if(!is.null(res$contrasts)) {
         res$contrasts <- pv.contrastDups(res$contrasts)	
      }
      if(length(res$contrasts) == numStart) {
         warning('Unable to add redundant contrast.',call.=FALSE)  	
      }
      return(res)
   }
   
   pv$contrasts <- pv.listaddto(pv$contrasts,res)
   
   return(pv)
   
}

pv.getContrasts <- function(pv,minMembers=3,
                            attributes=c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT),
                            block,bNot=FALSE){
   
   srcidx <-  pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
   mdata <- pv$class[,srcidx]
   
   if(!missing(block)) {
      #if(block != PV_REPLICATE) {
      #   warning('Unsupported blocking attribute')
      #}
      block <- pv.BlockList(pv,block)
   } else block <- NULL
   
   mcats <- attributes
   jobs <- NULL
   for(mcat in mcats) {
      vals <- unique(mdata[mcat,])
      if (length(vals)>1) {
         members <- NULL
         for(i in 1:length(vals)) {
            members <- c(members,sum(mdata[mcat,]== vals[i]))            
         }
         for(i in 1:length(vals)) {
            if(members[i] >= minMembers) {
               if(i < length(vals)){
                  for(j in (i+1):length(vals)) {
                     if(members[j] >= minMembers) {
                        job <- NULL
                        job$group1 <- pv.mask(pv,mcat,vals[i]) & srcidx
                        job$group2 <- pv.mask(pv,mcat,vals[j]) & srcidx
                        job$name1  <- vals[i]
                        job$name2  <- vals[j]
                        job$blocklist <- block
                        jobs <- pv.listadd(jobs,job)  
                     }
                  }      
               }
               if(bNot && ( (sum(members)-members[i]) >= minMembers) ) {
                  job <- NULL
                  job$group1 <- pv.mask(pv,mcat,vals[i]) & srcidx
                  job$group2 <- pv.mask(pv,mcat,vals[i],merge='nand') & srcidx
                  job$name1  <- vals[i]
                  job$name2  <- sprintf("!%s",vals[i])
                  job$blocklist <- block
                  jobs <- pv.listadd(jobs,job)                 
               }   
            }
         }     
      }
   }
   return(jobs)   
}

pv.contrastPairs <- function(pv,minMembers=3,
                             attributes=c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT),
                             block=NULL,conlist=NULL, bNot=TRUE) {
   
   if(length(attributes)==1) {
      return(pv$contrasts)
   }
   
   clist <- conlist
   srcmask <- pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
   numatts <- length(attributes)
   res <- NULL
   for(a1 in 1:(numatts-1)) {
      att1 <- attributes[a1]
      val1 <- unique(pv$class[att1,])
      if(length(val1)>1) {
         for(i in 1:length(val1)) {
            for (a2 in (a1+1):numatts) {
               att2 <- attributes[a2]
               val2 <- unique(pv$class[att2,])
               if(length(val2)>1) {
                  for(j in 1:length(val2)) {
                     res <- pv.listadd(res,list(a1=att1,v1=val1[i],a2=att2,v2=val2[j]))
                  }
               }
            }
         }
      }
   }
   
   if(!missing(block)) {
      block <- pv.BlockList(pv,block)
   } else block <- NULL
   
   if(!is.null(res)) {
      for(i in 1:length(res)) {
         m1 <- pv.mask(pv,res[[i]]$a1,res[[i]]$v1,mask=srcmask,merge='and')
         m1 <- pv.mask(pv,res[[i]]$a2,res[[i]]$v2,mask=m1,merge='and')
         if(sum(m1)>=minMembers) {
            if(i < length(res)) {
               for(j in (i+1):length(res)) {
                  m2 <- pv.mask(pv,res[[j]]$a1,res[[j]]$v1,mask=srcmask,merge='and')
                  m2 <- pv.mask(pv,res[[j]]$a2,res[[j]]$v2,mask=m2,merge='and') 
                  if( (sum(m2)>=minMembers) && (sum(m1 & m2) == 0) ) {
                     crec <- NULL
                     crec$group1 <- m1
                     crec$group2 <- m2
                     crec$name1  <- sprintf("%s:%s",res[[i]]$v1,res[[i]]$v2)
                     crec$name2  <- sprintf("%s:%s",res[[j]]$v1,res[[j]]$v2)
                     crec$blocklist <- block
                     clist <- pv.listadd(clist,crec)
                  }
               }                      
            }
            if(bNot && (sum((!m1)&srcmask) >= minMembers)) {
               crec <- NULL
               crec$group1 <- m1
               crec$group2 <- (!m1) & srcmask
               crec$name1  <- sprintf("%s:%s",res[[i]]$v1,res[[i]]$v2)
               crec$name2  <- sprintf("!%s:%s",res[[i]]$v1,res[[i]]$v2)
               crec$blocklist <- block
               clist <- pv.listadd(clist,crec)	
            }     
         }
      }	
   }
   #if(!is.null(block)){
   #   for(i in 1:length(clist)){
   #      if(!pv.checkBlock(clist[[i]])){
   #         warning("Unable to add blocking factor: unmatched samples.")
   #         clist[[i]]$blocklist <- NULL   		
   #      }
   #   }      	
   #} 
   return(clist)
}


pv.addContrast <- function(pv,group1,group2=!group1,name1="group1",name2="group2") {
   
   if(is.null(group1) | is.null(group2)) {
      stop('Null group, can not add contrast',call.=FALSE)	
   }
   
   if(!is.logical(group1)) {
      if(max(group1) > length(pv$peaks)) {
         stop('Invalid sample number in first group.',call.=FALSE)	
      }
      temp <- rep(FALSE,length(pv$peaks))
      temp[group1] <- TRUE
      group1 <- temp
   }
   
   if(!is.logical(group2)) {
      if(max(group2) > length(pv$peaks)) {
         stop('Invalid sample number in second group.',call.=FALSE)	
      }  	
      temp <- rep(FALSE,length(pv$peaks))
      temp[group2] <- TRUE
      group2 <- temp
   }
   
   if(sum(group1)==0) {
      return(pv)
   }
   if(sum(group2)==0) {
      return(pv)
   }
   
   if( length(group1) != length(pv$peaks) || length(group2) != length(pv$peaks) ) {
      stop('Length of vector specifying groups greater than number of samples.',call.=FALSE)
   }
   
   crec <- NULL
   crec$name1  <- name1
   crec$name2  <- name2
   crec$group1 <- group1
   crec$group2 <- group2
   pv$contrasts <- pv.listadd(pv$contrasts,crec)
   return(pv)	
}

pv.contrastDups <- function(clist) {
   numc <- length(clist)
   if(numc <= 1) {
      return(clist)
   }
   res <- rep(TRUE,numc)
   for(i  in 1:(numc-1)) {
      crec <- clist[[i]]
      for(j in (i+1):numc) {
         rec2 <- clist[[j]]
         if(!is.null(crec$contrast) && !is.null(rec2$contrast)) {
            if(identical(crec$contrast,rec2$contrast)) {
               res[j] <- FALSE
            }
         } else if (is.null(crec$contrast) && is.null(rec2$contrast)) {
            if ( (identical(crec$group1,rec2$group1) && identical(crec$group2,rec2$group2)) ||
                 (identical(crec$group1,rec2$group2) && identical(crec$group2,rec2$group1)) ){
               res[j] <- FALSE   
            }
         }
      }      	
   }
   newc <- NULL
   for(i in 1:numc) {
      if (res[i]) {
         newc <- pv.listadd(newc,clist[[i]])
      }	
   }
   return(newc)	
}


pv.BlockList <- function(pv,attribute=PV_REPLICATE) {
   
   if(is.numeric(attribute)) {
      if(length(attribute)>1) {
         vec <- rep(FALSE,length(pv$peaks))
         vec[attribute]=TRUE
         attribute=vec	
      }
   }
   if(is(attribute,'numeric')) {   
      vals    <- sort(unique(pv$class[attribute,]))
      attname <- rownames(pv$class)[attribute]
      if(attname == 'Peak caller') {
         attname <- 'Caller'
      }
      res <- NULL
      for(val in vals) {
         newmask <- pv.mask(pv,attribute,val)
         res <- pv.listadd(res,list(attribute=attname,label=val,samples=newmask))   
      }
   } else { #logical vector(s)
      if(is(attribute,'logical')) {
         if(length(attribute)!=length(pv$peaks)) {
            stop('Length of attribute vector must equal total number of peaksets.',call.=FALSE)	
         }
         attribute <- list(true=attribute,false=!attribute)	
      }
      if(!is(attribute,'list')) {
         stop('attribute must be a DBA_ attribute, a logical vector, or a list of logical vectors.',call.=FALSE)	
      }
      attname <- 'Block'
      if(is.null(names(attribute))) {
         names(attribute) <- 1:length(attribute)
      }
      res <- NULL
      hasatt <- rep(FALSE,length(pv$peaks))
      for (i in 1:length(attribute)) {
         att <- attribute[[i]]
         if(is.numeric(att)) {
            att <- rep(FALSE,length(pv$peaks))
            att[attribute[[i]]]=TRUE
         }
         hasatt <- hasatt | att
         res <- pv.listadd(res,list(attribute=attname,label=names(attribute)[i],samples=att))	
      }
      if(sum(!hasatt)>0) {
         res <- pv.listadd(res,list(attribute=attname,label="other",samples=!hasatt))	
      }
   }
   
   return(res)	
}

pv.checkBlock <- function(contrast,bCheckBalanced=FALSE,bCheckMultiple=TRUE,bCheckCross=TRUE,
                          bCheckUnique=TRUE,bCheckAllInOne=TRUE,bWarning=TRUE) {
   
   if(bCheckBalanced){
      if(sum(contrast$group1)!=sum(contrast$group2)) {
         if(bWarning) warning("Blocking factor has unmatched sample(s).",call.=FALSE)
         return(FALSE)	
      }
      
      for(att in contrast$blocklist) {
         if(sum(contrast$group1 & att$samples) != sum(contrast$group2 & att$samples)) {
            if(bWarning) warning("Blocking factor has unmatched sample(s).",call.=FALSE)
            return(FALSE)
         }
      }
   }
   
   if(bCheckMultiple) {
      if(length(contrast$blocklist)<2) {
         if(bWarning) warning('Blocking factor has only one value',call.=FALSE)
         return(FALSE)	
      }	
   }
   
   if(bCheckCross || bCheckAllInOne) {
      cross  <- FALSE
      allone <- FALSE
      for(att in contrast$blocklist) {
         if(sum(contrast$group1 & att$samples) & sum(contrast$group2 & att$samples)) {
            cross <- TRUE
         }
         if(sum(att$samples & (contrast$group1 | contrast$group2)) ==
            sum(contrast$group1 | contrast$group2))
            allone=TRUE
      }
      if(!cross && bCheckCross) {
         if(bWarning) warning('No blocking values are present in both groups',call.=FALSE)	
         return(FALSE)
      }
      if(allone && bCheckAllInOne) {
         if(bWarning) warning('A blocking value applies to all samples in both groups',call.=FALSE)	
         return(FALSE)
      }
   }
   
   if(bCheckUnique) {
      unique <- rep(0,length(contrast$group1)) 	
      for(att in contrast$blocklist) {
         unique <- unique + (contrast$group1 & att$samples)
         unique <- unique + (contrast$group2 & att$samples)	
      }
      if(sum(unique>1)>0) {
         if(bWarning) warning('Some sample(s) have more than one value for blocking factor',call.=FALSE)	
         return(FALSE)
      }
   }
   
   return(TRUE)
}

pv.listContrasts <- function(pv,th=0.05,bUsePval=FALSE) {
   
   if(is.null(pv$contrasts)) {
      return(NULL)
   } else {
      clist <- pv$contrasts
   }
   
   types  <- pv.contrastTypes(clist)
   groups <- pv.contrastGroups(clist)
   coeffs <- pv.contrastCoeffs(clist, pv)
   blocks <- pv.contrastBlocks(clist)
   sdbs   <- pv.contrastSDBs(clist,th=th)
   
   res <- NULL
   if(!is.null(types)) {
      res <- matrix(types,length(types),1)
      colnames(res) <- "Factor"
      res <- cbind(res,groups)
   } else {
      res <- groups
   }
   if(!is.null(coeffs)) {
      res <- cbind(res,coeffs)
   }
   if(!is.null(blocks)) {
      res <- cbind(res,blocks)
   }
   if(!is.null(sdbs)) {
      res <- cbind(res,sdbs)
   }
   
   rownames(res) <- 1:nrow(res)
   return(data.frame(res))
   
}

pv.contrastTypes <- function(contrasts) {
   res <- NULL
   for(cres in contrasts) {
      type <- cres$contrastType
      if(is.null(type)) {
         type <- ""
      } else if (type =="simple") {
         type <- cres$contrast[1]
      } else if (type =="byname") {
         type <- "Name"
      } else if (type =="byresults1" || type=="byresults2") {
         type <- "Results"
      } else if (type =="bycolumn") {
         type <- "Coefficient"
      }
      res <- c(res,type)
   }
   if(sum(res==rep("",length(res))) != length(res)) {
      return(res)
   } else {
      return(NULL)
   }
}

pv.contrastGroups <- function(contrasts) {
   
   res <- NULL
   for(crec in contrasts) {
      simple <- FALSE
      if(!is.null(crec$contrastType)) {
         if(crec$contrastType == 'simple') {
            simple <- TRUE
         }
      } else {
         simple <- TRUE
      }
      
      if(simple) {
         newrec <- c(crec$name1,sum(crec$group1),crec$name2,sum(crec$group2))
      } else { # complex contrast
         if(crec$contrastType=='byname') {
            newrec <- crec$contrast
         } else if(crec$contrastType=="byresults1") {
            newrec <- c(crec$contrast[[1]])
         } else if(crec$contrastType=="byresults2") {
            newrec <- c(crec$contrast[[1]],"",crec$contrast[[2]])
         } else if(crec$contrastType=="bycolumn") {
            newrec <- paste(crec$contrast,collapse=",")
         } else {
            stop('Invalid contrast type',call.=FALSE)
         }
         while(length(newrec)<4){
            newrec <- c(newrec,"")
         }
      }
      
      if(!is.null(res)) {
         if(length(newrec)>ncol(res)) {
            res <- cbind(res,matrix("-",nrow(res),length(newrec)-ncol(res)))	
         } else if (length(newrec)<ncol(res)) {
            newrec <- c(newrec,rep("-",ncol(res)-length(newrec)))	
         }
      }
      res <- rbind(res,newrec)
   }
   
   numCons <- length(contrasts)
   rownames(res) <- 1:numCons
   colnames(res) <- c("Group","Samples","Group2","Samples2")
   
   res <- pv.keepMeta(res, "")
   
   return(res)
}

pv.contrastCoeffs <- function(contrasts, pv) {
   coeffs <- cons <- NULL
   for(i in 1:length(contrasts)) {
      if(!is.null(contrasts[[i]]$contrastType)) {
         if(contrasts[[i]]$contrastType == "bycolumn") {
            coeffs <- c(coeffs,i)
         } else {
            cons <- c(cons,i)
         }
      }   
   }
   if(is.null(coeffs)) {
      return(NULL)
   }
   
   numContrasts <- length(contrasts)
   
   coeffNames <- dba.contrast(pv,bGetCoefficients=TRUE)
   
   res <- matrix("_",numContrasts,length(coeffNames))
   colnames(res) <- coeffNames
   rownames(res) <- 1:numContrasts
   
   for(i in coeffs) {
      res[i,] <- contrasts[[i]]$contrast
   }  
   
   if(!is.null(cons)) {
      edgeRsave <- pv$edgeR
      pv$edgeR$design <- pv.edgeRdesign.matrix(pv,pv$design)
      for(i in cons) {
         con <- contrasts[[i]]
         coeffs <- pv.getCoeffs(con,pv)
         if(!is.null(coeffs)){
            res[i,] <- coeffs
         }
      }
      pv$edgeR <- edgeRsave
   }
   
   return(res)
}

pv.contrastBlocks <- function(contrasts) {
   
   res <- NULL
   maxbvals <- 0
   for(crec in contrasts) {
      bvals <- 0
      newrec <- NULL
      if(!is.null(crec$blocklist)) {
         for(brec in crec$blocklist) {
            bvals <- bvals+1
            newrec <- c(newrec,brec$label,sum(brec$samples&(crec$group1|crec$group2)))
         }
      }
      if(is.null(newrec)) {
         newrec <- "_"
      }
      maxbvals <- max(bvals,maxbvals)
      if(!is.null(res)) {
         if(length(newrec)>ncol(res)) {
            res <- cbind(res,matrix("-",nrow(res),length(newrec)-ncol(res)))
         } else if (length(newrec)<ncol(res)) {
            newrec <- c(newrec,rep("-",ncol(res)-length(newrec)))
         }
      }
      res <- rbind(res,newrec)
   }
   
   if(maxbvals==0) {
      return(NULL)
   } else {
      cvec <- NULL
      for(i in 1:maxbvals) {
         cvec <- c(cvec,paste("Block",i,sep=""),paste("Blk",i,"Samps",sep=""))
      }
      colnames(res) <- cvec
      return(res)
   }
}

EDGER_COL_PVAL <- 4
EDGER_COL_FDR  <- 5
pv.contrastSDBs <- function(contrasts,th=th) {
   res  <- NULL
   cvec <- NULL
   eres <- NULL
   for(crec in contrasts) {
      if(!is.null(names(crec$edgeR))){
         if(!is(crec$edgeR,"DGEList")) {
            eres <- c(eres,sum(crec$edgeR$de$padj <= th))
         } else {
            if(is.null(crec$edgeR$LRT)) {
               eres <- c(eres,sum(edgeR::topTags(crec$edgeR$db,
                                                 nrow(crec$edgeR$db$counts))$table[,EDGER_COL_FDR] <= th,na.rm=TRUE))
            } else {
               eres <- c(eres,sum(edgeR::topTags(crec$edgeR$LRT,
                                                 nrow(crec$edgeR$db$counts))$table[,EDGER_COL_FDR+1] <= th,na.rm=TRUE))
            }
         }
      } else {
         eres <- c(eres,"-")
      }
   }
   res <- cbind(res,eres)
   cvec <- c(cvec,'DB edgeR')
   
   eres <- NULL
   for(crec in contrasts) {
      if(!is.null(names(crec$edgeR$block))){
         eres <- c(eres,sum(edgeR::topTags(crec$edgeR$block$LRT,
                                           nrow(crec$edgeR$block$counts))$table[,EDGER_COL_FDR+1] <= th,na.rm=TRUE))
      } else {
         eres <- c(eres,"-")
      }
   }
   res <- cbind(res,eres)
   cvec <- c(cvec,'DB edgeR-block')
   
   eres <- NULL
   for(crec in contrasts) {
      if(!is.null(names(crec$DESeq2)) && 
         (!is(crec$DESeq2,"try-error")) ){
         eres <- c(eres,sum(crec$DESeq2$de$padj <= th,na.rm=TRUE))
      } else {
         eres <- c(eres,"-")
      }
   }
   res <- cbind(res,eres)
   cvec <- c(cvec,'DB DESeq2')
   
   eres<- NULL
   for(crec in contrasts) {
      if(!is.null(names(crec$DESeq2$block)) && 
         (!is(crec$DESeq2,"try-error")) ){
         eres <- c(eres,sum(crec$DESeq2$block$de$padj <= th,na.rm=TRUE))
      } else {
         eres <- c(eres,"-")
      }
   }
   res <- cbind(res,eres)
   cvec <- c(cvec,'DB DESeq2-block')
   colnames(res) <- cvec
   
   res <- pv.keepMeta(res,"-")
   
   return(res)
}

pv.keepMeta <- function(res, noval="")
{
   numContrasts <- nrow(res)
   cvec <- colnames(res)
   
   keep  <- NULL
   keepc <- NULL
   for(i in 1:ncol(res)) {
      if(sum(res[,i]==rep(noval,numContrasts)) != numContrasts) {
         keep  <- cbind(keep,res[,i])
         keepc <- c(keepc,cvec[i])
      }
   }
   
   if(is.null(keep)) {
      return(NULL)
   }
   
   res <- data.frame(keep)
   rownames(res) <- 1:numContrasts
   colnames(res) <- keepc
   
   return(res)
}

pv.design <- function(DBA,categories=c(DBA_CONDITION,DBA_TREATMENT,DBA_TISSUE,DBA_FACTOR,DBA_REPLICATE)) {
   
   facs <- NULL
   for(cond in categories) {
      if(length(unique(DBA$class[cond,]))>1) {
         facs <- c(facs,cond)
      } else {
         catname <- "UNKNOWN"
         if(cond == DBA_CONDITION) {
            catname <- "DBA_CONDITION"
         }	
         if(cond == DBA_TREATMENT) {
            catname <- "DBA_TREATMENT"
         }	
         if(cond == DBA_TISSUE) {
            catname <- "DBA_TISSUE"
         }	
         if(cond == DBA_FACTOR) {
            catname <- "DBA_FACTOR"
         }	
         if(cond == DBA_REPLICATE) {
            catname <- "DBA_REPLICATE"
         }
         warning(sprintf("Category %s lacks multiple unique values, ignored in design",catname),call.=FALSE)	
      }	
   }
   
   red <- DBA$class[facs,]
   if(is.null(nrow(red))) {
      red <- matrix(red,1,length(red))
      rownames(red) <- rownames(DBA$class)[facs]
      colnames(red) <- colnames(DBA$class)
   }
   dmatrix <- data.frame(t(red))
   terms <- colnames(dmatrix)
   form <- "~"
   for(term in terms) {
      form <- sprintf("%s+%s",form,term)	
   }
   
   design <- model.matrix(as.formula(form),data=dmatrix)
   
   return(design)   	
}

pv.FactorNames <- c("Tissue","Factor","Condition","Treatment","Caller","Replicate")
pv.FactorVals  <- c(PV_TISSUE, PV_FACTOR, PV_CONDITION,PV_TREATMENT,PV_CALLER,PV_REPLICATE)

pv.contrastDesign <- function(pv, design, contrast, name1, name2, bGetNames=FALSE){
   
   # Check if contrast already exists
   if(!missing(contrast)) {
      for(con in pv$contrasts){
         if(identical(con$contrast,contrast)) {
            message("Contrast already present.")
            return(pv)
         }
      }
   }
   
   # Check if design already exists
   newdesign <- FALSE
   if(!missing(design)) {
      if(is.logical(design)) {
         if(is.character(contrast) && length(contrast)==3) {
            design <- pv.generateSimpleDesign(pv,addfactor=contrast[1])
         } else {
            design <- pv$design
         }
      } 
      if(pv.isDesign(pv$design)) {
         if(design != pv$design) {
            if(!is.null(pv$DESeq2)) {
               newdesign <- TRUE
               if(!is.null(pv$DESeq2$DEdata)) {
                  message("Replacing design and removing analysis.")
               } else {
                  message("Replacing design.")
               }
               pv$DESeq2 <- pv$edgeR <- NULL
            } else {
               message("Replacing design.")
            }
         }
      }
   } else {
      if(is.null(pv$design)) {
         stop("No design specified.",call.=FALSE)
      } else {
         design <- pv$design
      }
   }
   
   varnames <-  pv.checkDesign(design)
   pv$design <- design
   
   if(is.null(pv$DESeq2$names)) {
      pv$DESeq2$names <- pv.DESeq2ResultsNames(pv)
   }
   
   # Check if contrasts exist and are compatible with design
   if(!is.null(pv$contrasts)) {
      for(i in 1:length(pv$contrasts)) {
         if(!is.null(pv$contrasts[[i]]$contrast)) {
            pv$contrasts[[i]]$contrastType <- 
               pv.checkContrast(pv,varnames,pv$contrasts[[i]]$contrast)
            if(newdesign) {
               pv$contrasts[[i]]$DESeq2 <- NULL
               pv$contrasts[[i]]$edgeR  <- NULL
            }
         }
      }
   }
   
   if(bGetNames) {
      if(is.null(pv$design)) {
         stop("Can't get names: no design.",call.=FALSE)
      }
      return(pv$DESeq2$names)
   }
   
   # Check new contrast and add   
   if(!missing(contrast)) {
      crec <- NULL
      crec$contrastType <- pv.checkContrast(pv,varnames,contrast)
      if(!is.null(crec$contrastType)) {
         crec$contrast <- contrast
         if(missing(name1)) {
            if(crec$contrastType=='simple') {
               crec$name1 = contrast[2]
            } else if(crec$contrastType=='byname') {
               crec$name1 = contrast[1]
            } else if(crec$contrastType=='byresults1' || crec$contrastType=='byresults2') {
               crec$name1 = contrast[[1]]
            }
         } else crec$name1 <- name1
         if(missing(name2)) {
            if(crec$contrastType=='simple') {
               crec$name2 = contrast[3]
            } else if(crec$contrastType=='byresults2') {
               crec$name2 = contrast[[2]]
            }
         } else crec$name2 <- name2
         if(crec$contrastType=='simple') {
            crec$group1 <- pv.getSamplesForFactorVal(pv,contrast[1],contrast[2])
            crec$group2 <- pv.getSamplesForFactorVal(pv,contrast[1],contrast[3])    
            if(sum(crec$group1)<=1 || sum(crec$group2)<=1) {
               stop("Invalid contrast: no replicates in one group.")
            }
         }   
         pv$contrasts <- pv.listadd(pv$contrasts,crec) 
      } else {
         warning("Empty contrast",call.=FALSE)
      }
   }
   
   # Return DBA object
   return(pv)
}

pv.checkDesign <- function(design) {
   if(is.null(design)) {
      return(NULL)
   }
   varnames <- all.vars(formula(design))
   matches <- varnames %in% pv.FactorNames
   if(sum(matches) != length(varnames)) {
      message <- paste("Invalid factors in design:",paste(varnames[!matches],collapse=','))
      stop(message,call.=FALSE)
   }
   return(varnames)
}

pv.compareDesigns <- function(oldd, newd){
   oldv <-  all.vars(formula(oldd))
   newv <-  all.vars(formula(newd))
   if(all(newv %in% oldv)) {
      return(oldd)
   } else {
      return(newd)
   }
}

pv.DESeq2ResultsNames <- function(pv) {
   if(is.null(pv$design)) {
      stop('No design specified',call.=FALSE)
   }
   if(is.null(pv$DESeq2$names)) {
      message("Computing results names...")
      res <- suppressMessages(pv.DEinitDESeq2(pv,numReads=10))
      DESeq2::sizeFactors(res) <- rep(1,ncol(res))
      res <- suppressWarnings(suppressMessages(
         DESeq2::estimateDispersionsGeneEst(res,quiet=TRUE,maxit=10)))
      DESeq2::dispersions(res) <- S4Vectors::mcols(res)$dispGeneEst
      res <- suppressWarnings(suppressMessages(
         DESeq2::nbinomWaldTest(res, maxit=10, quiet=TRUE,useOptim=FALSE)))      
      names <- DESeq2::resultsNames(res)
   } else {
      names <- pv$DESeq2$names
   }
   return(names)
}

pv.checkContrast <- function(pv,varnames,contrast) {
   
   # What type of contrast description is this?
   simple <- byname <- bycolumn <- byresults1 <- byresults2 <- FALSE
   if(is(contrast,"character")) {
      if(length(contrast) == 3) {
         simple <- TRUE
      } else if (length(contrast)== 1) {
         byname <- TRUE
      } else {
         stop('Invalid contrast. Character vector must be length 1 or 3.',call.=FALSE)
      }
   } else if(is(contrast,"numeric") || is(contrast,"integer")) {
      bycolumn <- TRUE
   } else if(is(contrast,"list")) { 
      if(length(contrast)==1) {
         byresults1 <- TRUE    
      } else if (length(contrast)==2) {
         byresults2 <- TRUE
      } else {
         stop('Invalid contrast. List must be length 1 or 2.',call.=FALSE)
      }
   } else {
      stop('Invalid contrast.',call.=FALSE)
   }
   
   if(simple) {
      facName <- contrast[1]
      if(!facName %in% varnames) {
         message <- paste("Contrast factor not present in design:",contrast[1])
         stop(message,call.=FALSE)
      }
      for(facVal in contrast[2:3]) {
         if (!sum(facVal %in% pv$class[facName,])) {
            message <- paste("Contrast value not valid for factor:",facVal)
         }
      }
      return('simple')
   } else {
      names <- pv.DESeq2ResultsNames(pv)
      if(byname) {
         if (!contrast %in% names ) {
            stop('Invalid contrast. Single name must match one of design column names.',call.=FALSE)
         }
         return('byname')
      } else if(bycolumn) {
         if(length(contrast) != length(names)) {
            stop('Invalid contrast. Numeric vector must be same length as number of columns in design.',call.=FALSE)
         }
         return('bycolumn')
      } else if(byresults1 || byresults2) {
         if(!contrast[[1]] %in% names) {
            stop('Invalid contrast. Names in list must match one of design column names.',call.=FALSE)
         }  
         if(byresults1) {
            return('byresults1')
         } else if(byresults2) {
            if(!contrast[[2]] %in% names) {
               stop('Invalid contrast. Names in list must match one of design column names.',call.=FALSE)
            }  
            return('byresults2')
         } 
      }
   }
   
   return(NULL)
}

pv.getSamplesForFactorVal <- function(pv,facID,facVal){
   samples <-  pv$class[facID,] %in% facVal
   names(samples) <- pv$class[PV_ID,]
   return(samples)
}

pv.generateDesignFromContrasts <- function(pv,attributes=pv.FactorVals,
                                           bGetNames=FALSE,addfac=NULL) {
   factors <- addfac
   for(con in 1:length(pv$contrasts)) {
      cont <- pv$contrasts[[con]]
      factor1 <- NULL
      factor2 <- NULL
      for(fac in attributes) {
         if(cont$name1 %in% pv$class[fac,]) {
            factor1 <- fac
         }
         if(cont$name2 %in% pv$class[fac,]) {
            factor2 <- fac
         }
      }
      if(is.null(factor1) || is.null(factor2)) {
         stop('Unable to generate design, unknown values',call.=FALSE)
      }
      if(factor1 != factor2) {
         stop('Unable to generate design, unknown values',call.=FALSE)            
      }
      factors <- unique(c(factors,factor1))
      pv$contrasts[[con]]$contrast <- c(pv.FactorNames[match(factor1,pv.FactorVals)],
                                        cont$name1,cont$name2)
      pv$contrasts[[con]]$contrastType <- 'simple'
   }
   
   factors <- pv.FactorNames[match(factors,pv.FactorVals)]
   pv$design <- paste("~",paste(factors,collapse =" + "),sep="")
   if(bGetNames) {
      pv$DESeq2$names <- pv.DESeq2ResultsNames(pv)
   }
   return(pv)
}

pv.generateSimpleDesign <- function(pv,addfactor) {
   if(!is.null(pv$design)) {
      varnames <- pv.checkDesign(pv$design)
   } else {
      varnames <- NULL
   }
   if(addfactor %in% varnames) {
      return(pv$design)
   }
   varnames <- unique(c(varnames,addfactor))
   design <- paste("~",paste(varnames,collapse =" + "),sep="")
   return(design)
}

pv.getContrastString <- function(conrec, bFlip=FALSE) {
   if(!is.null(conrec$contrastType)) {
      constr <- conrec$contrast
      if(conrec$contrastType=='bycolumn') {
         constr <- paste(constr,collapse=":")
      } else if(conrec$contrastType=='byresults1') {
         constr <- constr[[1]]
      } else if(conrec$contrastType=='byresults2') {
         if(bFlip==TRUE) {
            constr <- paste(constr[[2]],"vs.",constr[[3]])
         } else {
            constr <- paste(constr[[1]],"vs.",constr[[2]])
         }
      } else if(length(constr) == 2) {
         if(bFlip==TRUE) {
            constr <- paste(constr[2],"vs.",constr[1 ])   
         } else {
            constr <- paste(constr[1],"vs.",constr[2])      
         }
      } else if(length(constr) == 3)  {
         if(bFlip==TRUE) {
            constr <- paste(constr[3],"vs.",constr[2])   
         } else {
            constr <- paste(constr[2],"vs.",constr[3])    
         }
      }
   } else {
      if(bFlip==TRUE) {
         constr <- paste(conrec$name2, "vs.", conrec$name1)         
      } else {
         constr <- paste(conrec$name1, "vs.", conrec$name2)
      }
   }
   return(constr)
}

pv.orderMetaFacs <- function(pv, reorderMeta) {
   meta <- pv$class[pv.FactorVals,]
   meta <- data.frame(t(meta),stringsAsFactors = TRUE)
   
   if(is.null(pv$meta)) { # First time
      for(facname in pv.FactorNames) {
         if(facname=="Caller") {
            facname <- "Peak.caller"
         }
         mfactor <- match(facname,names(meta))
         vals <- as.character(meta[,mfactor])
         pv$meta <- pv.listadd(pv$meta,unique(vals))
      }
      names(pv$meta) <- pv.FactorNames
   }
   
   if(!missing(reorderMeta)) {
      if(is(reorderMeta,"list")) {
         for(i in 1:length(reorderMeta)) {
            facname <- names(reorderMeta)[i]
            mfactor <- match(facname,names(pv$meta))
            if(is.na(mfactor)) {
               stop("Unknown metadata factor: ",facname,call.=FALSE)
            }
            vals <- unique(as.character(meta[,mfactor]))
            facs <- reorderMeta[[i]]
            matches <- match(facs,vals)
            if(sum(is.na(matches))) {
               stop("Invalid values for factor ",facname,call.=FALSE)
            }
            if(length(vals) > length(facs)) {
               vals <- vals[-matches]
               facs <- c(facs,vals)
            }
            pv$meta[[mfactor]] <- facs
         }
         pv$DESeq2 <- pv$edgeR <- NULL
         #pv <- pv.CheckReleveledContrasts(pv)
      } else {
         stop("reorderMeta must be a list.",call.=FALSE)
      }
   }
   return(pv)
}

pv.getMeta <- function(pv) {
   
   if(is.null(pv$meta)) {
      pv <- pv.orderMetaFacs(pv)
   }
   
   meta <- pv$class[pv.FactorVals,]
   meta <- data.frame(t(meta),stringsAsFactors = TRUE) 
   
   for(i in 1:ncol(meta)) {
      meta[,i] <- factor(as.character(meta[,i]),pv$meta[[i]])
   }
   
   return(meta)
}

pv.CheckReleveledContrasts <- function(pv) {
   
   if(is.null(pv$contrasts) || !pv.isDesign(pv$design)) {
      return(pv)
   }
   
   pv$DESeq2$names <- names <- pv.DESeq2ResultsNames(pv)
   toremove <- NULL
   for(i in 1:length(pv$contrasts)) {
      con <- pv$contrasts[[i]]
      if(!is.null(con$contrastType)) {
         if(con$contrastType=="bycolumn") {
            toremove <- c(toremove,i)   
         } else if(con$contrastType=="byname") {
            if (!con$contrast %in% names ) {
               toremove <- c(toremove,i)
            }
         } else if(con$contrastType=="byresults1" || con$contrastType=="byresults2") {
            if(!con$contrast[[1]] %in% names) {
               toremove <- c(toremove,i)
            } else if(con$contrastType=="byresults2") {
               if(!con$contrast[[2]] %in% names) {
                  toremove <- c(toremove,i)
               }
            }
         }
      }
   }
   
   if(length(toremove)>0) {
      toremove <- sort(unique(toremove))
      for(i in toremove) {
         warning("Removing contrast: ",i,call.=FALSE)
      }
      pv$contrasts <- pv$contrasts[-toremove]
   }
   
   return(pv)
} 

pv.isDesign <- function(design) {
   if(is.null(design)) {
      return(FALSE)
   }
   
   if(is(design,"logical")) {
      return(design)
   }
   
   return(TRUE)
}

pv.checkContrastOrder <- function(contrasts, meta) {
   for(i in length(contrasts)) {
      con <- contrasts[[i]]
      for(fac in meta) {
         n1 <- fac %in% con$name1
         n2 <- fac %in% con$name2
         if(sum(n1) && sum(n2)) {
            if(which(n2) > which(n1)) {
               g1 <- con$group1
               n1 <- con$name1
               con$group1 <- con$group2
               con$name1  <- con$name2
               con$group2 <- g1
               con$name2  <- n1
               contrasts[[i]] <- con
            }
         }
      }
   }
   return(contrasts)
}
