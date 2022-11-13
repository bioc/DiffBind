#################GENERIC PARALLEL INTERFACE #################

DBA_PARALLEL_MULTICORE <- 1
DBA_PARALLEL_RLSF      <- 2
DBA_PARALLEL_LSFR      <- 3
DBA_PARALLEL_BIOC      <- 4

### INITIALIZE ###
dba.parallel <- function(DBA) {
   
   if(DBA$config$parallelPackage==DBA_PARALLEL_MULTICORE) {
      DBA$config <-  dba.multicore.init(DBA$config)
      DBA$config$parallelInit <- TRUE
      return(DBA)
   }
   
   if(DBA$config$parallelPackage==DBA_PARALLEL_BIOC) {
      DBA$config <-  dba.biocparallel.init(DBA$config)
      DBA$config$parallelInit <- TRUE
      return(DBA)
   }
   # 
   # if(DBA$config$parallelPackage==DBA_PARALLEL_RLSF) {
   #    DBA$config <-  dba.Rlsf.init(DBA$config)
   #    DBA$config$parallelInit <- TRUE
   #    return(DBA)
   # }
   # 
   # if(DBA$config$parallelPackage==DBA_PARALLEL_LSFR) {
   #    DBA$config <-  dba.lsfR.init(DBA$config)
   #    DBA$config$parallelInit <- TRUE
   #    return(DBA)
   # }
   # 
   if(DBA$config$parallelPackage > 0) {
      warning('UNKNOWN PARALLEL PACKAGE: ',DBA$config$parallelPackage,call.=FALSE)
   }
   return(DBA)
}

dba.parallel.init <- function(config){
   config <- config$initFun(config)
   return(config)
}

### CONFIG ###
dba.parallel.params <- function(config,funlist) {
   params <- config$paramFun(config,funlist)
   return(params)
}

### PARALLEL LAPPLY ###
dba.parallel.lapply <- function(config,params,arglist,fn,...) {
   jobs <- config$lapplyFun(config,params,arglist,fn,...)
   return(jobs)
}

### ADD JOB ###
dba.parallel.addjob <- function(config,params,fn,...) {
   job <- config$addjobFun(config,params,fn,...)
   return(job)
}

### WAIT FOR JOBS TO COMPLETE AND GATHER RESULTS ###
dba.parallel.wait4jobs <- function(config,joblist) {
   res <- config$wait4jobsFun(config,joblist)
   return(res)
}

# dba.Rlsf.init <- function(config){
#    if (length(find.package(package="DiffBindCRI",quiet=TRUE))>0) {
#       config <- dba.CRI.Rlsf.init(config)    
#    } else {
#       warning('Rlsf interface not supported in this version')
#    }
#    
#    return(config)
# }
# 
# dba.lsfR.init <- function(config){
#    if (length(find.package(package="DiffBindCRI",quiet=TRUE))>0) {
#       config <- dba.CRI.lsfR.init(config)    
#    } else {
#       warning('lsfR interface not supported in this version')
#    }
#    
#    return(config)
# }

################# multicore INTERFACE #################

### INITIALIZE ###

dba.multicore.init <- function(config) {
   
   noparallel <- FALSE
   if (requireNamespace("parallel",quietly=TRUE)) {
      if(!exists("mcparallel")) {
         noparallel <- TRUE
      }     
      if(!exists("mclapply")) {
         noparallel <- TRUE
      } 
   } else {
      noparallel <- TRUE
   }
   if(noparallel){
      warning("Parallel execution unavailable: executing serially.",call.=FALSE)
      config$RunParallel <- FALSE
      config$parallelPackage <- 0
      return(config)
   }
   
   config$multicoreInit <- TRUE    
   
   config$initFun      <- dba.multicore.init
   config$paramFun     <- dba.multicore.params
   config$addjobFun    <- dba.multicore.addjob
   config$lapplyFun    <- dba.multicore.lapply
   config$wait4jobsFun <- dba.multicore.wait4jobs
   
   if(is.null(config$cores)) {
      config$cores <- dba.multicore.setCores()
   }
   
   return(config)
}

### CONFIG ###
dba.multicore.params <- function(config,funlist) {
   return(NULL)
}

### PARALLEL LAPPLY ###
dba.multicore.lapply <- function(config,params,arglist,fn,...){
   savecores <- options("cores")$cores
   savemccores <- options("mc.cores")$mc.cores
   options(cores=config$cores)
   options(mc.cores=config$cores)
   res <- mclapply(arglist,fn,...,mc.preschedule=TRUE,mc.allow.recursive=TRUE)
   options(cores=savecores)
   options(mc.cores=savemccores)
   return(res)
}

### ADD JOB ###
dba.multicore.addjob <- function(config,params,fn,...) {
   job <- parallel::mcparallel(fn(...))
   return(job)
}

### WAIT FOR JOBS TO COMPLETE AND GATHER RESULTS ###
dba.multicore.wait4jobs <- function(config,joblist) {
   res <- parallel::mccollect(joblist)
   return(res)
}

### set cores if none specified ###
dba.multicore.setCores <- function(config=NULL) {
  if(!is.null(config$cores)) {
    return(config$cores)
  }
  
  if(!is.null(options("mc.cores")$mc.cores)) {
    return(options("mc.cores")$mc.cores)
  }
  
  if(!is.null(options("cores")$cores)) {
    return(options("cores")$cores)
  }
  
  return(parallel::detectCores(logical=FALSE))
  
}

################# BiocParallel INTERFACE #################

### INITIALIZE ###

dba.biocparallel.init <- function(config) {
   
   noparallel <- FALSE
   if (!requireNamespace("BiocParallel",quietly=TRUE)) {
      noparallel <- TRUE
   }
   if(noparallel){
      warning("Parallel execution unavailable: executing serially.",call.=FALSE)
      config$RunParallel <- FALSE
      config$parallelPackage <- 0
      return(config)
   }
   
   config$biocparallelInit <- TRUE    
   
   config$initFun      <- dba.biocparallel.init
   config$paramFun     <- dba.biocparallel.params
   config$addjobFun    <- dba.biocparallel.addjob
   config$lapplyFun    <- dba.biocparallel.lapply
   config$wait4jobsFun <- dba.biocparallel.wait4jobs
   
   return(config)
}

### CONFIG ###
dba.biocparallel.params <- function(config,funlist) {
   return(NULL)
}

### PARALLEL LAPPLY ###
dba.biocparallel.lapply <- function(config,params,arglist,fn,...){
   if(!is.null(config$bpparam)) {
      res <- bplapply(arglist,fn,...,BPPARAM=config$bpparam)        
   } else {
      res <- bplapply(arglist,fn,...)        
   }
   return(res)
}

### ADD JOB ###
dba.biocparallel.addjob <- function(config,params,fn,...) {
   stop('addjob unimplemented when using BiocParallel.',call.=FALSE)
}

### WAIT FOR JOBS TO COMPLETE AND GATHER RESULTS ###
dba.biocparallel.wait4jobs <- function(config,joblist) {
   stop('wait4jobs unimplemented when using BiocParallel.',call.=FALSE)
}
