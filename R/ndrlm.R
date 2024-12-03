#-----------------------------------------------------------------------------#
#                                                                             #
#  GENERALIZED NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (GNDA)     #
#                                                                             #
#  Written by: Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona,      #
#              Zahid Khan                                                     #
#              *Department of Quantitative Methods                            #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: November 2024                                                #
#-----------------------------------------------------------------------------#
### GENERALIZED NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (GNDR) ##
#' @export
ndrlm<-function(Y,X,optimize=TRUE,cor_method=1,cor_type=1,min_comm=2,Gamma=1,
              null_model_type=4,mod_mode=6,use_rotation=FALSE,
              rotation="oblimin",pareto=TRUE,out_weights=rep(1,ncol(Y)),
              lower.bounds = c(rep(-100,ncol(X)),0,0,0,0),
              upper.bounds = c(rep(100,ncol(X)),0.6,0.6,0.6,0.3),
              popsize = 20, generations = 30, cprob = 0.7, cdist = 5,
              mprob = 0.2, mdist=10, seed=NULL){
  cl<-match.call()
  if (!requireNamespace("mco", quietly = TRUE)) {
    stop(
      "Package \"mco\" must be installed to use this function.",
      call. = FALSE
    )
  }
  Y<-as.data.frame(Y)
  X<-as.data.frame(X)
  hyperparams<-c(rep(1,ncol(X)),0,0,0,0)
  weight<-rep(1,ncol(X))
  cost<-function(hyperparams){
    weight<-hyperparams[1:ncol(X)]
    params<-hyperparams[-c(1:ncol(X))]
    NDA<-try(ndr(X,min_evalue = params[1],
                 min_communality=params[2],
                 com_communalities = params[3],
                 min_R = params[4],weight=weight,covar=FALSE,
                 cor_method=cor_method,
                 cor_type=cor_type,
                 min_comm=min_comm,
                 Gamma=Gamma,
                 null_model_type=null_model_type,
                 mod_mode=mod_mode,
                 use_rotation=use_rotation,
                 rotation=rotation),silent=TRUE)
    res<-c(rep(0,ncol(Y)))
    if (inherits(NDA,"try-error")){
      if (pareto==TRUE){
        return(res)
      }else{
        return(0)
      }
    }else{
      for (i in 1:ncol(Y))
      {
        data<-cbind(Y[,i],NDA$scores)
        colnames(data)[1]<-colnames(Y)[i]
        colnames(data)[-1]<-paste("NDA",1:NDA$factors,sep="")
        data<-as.data.frame(data)
        fit<-stats::lm(str2lang(paste(colnames(data)[1],"~",
                                      gsub(", ","+",
                                           toString(colnames(data)[-1])))),data)
        res[i]<-stats::summary.lm(fit)$adj.r.squared
      }
      if (pareto==TRUE){
        return(res)
      }else{
        return(stats::weighted.mean(res,out_weights))
      }
    }
  }
  if (!is.null(seed))
  {
    set.seed(seed)
  }

  costmin <- function(hyperparams) -cost(hyperparams)
  if (pareto==TRUE){
    ODIM<-ncol(Y)
  }else{
    ODIM<-1
  }
  if (optimize==TRUE)
  {
    NSGA <- mco::nsga2(fn=costmin,idim=length(hyperparams),odim=ODIM,
                       lower.bounds = lower.bounds,
                       upper.bounds = upper.bounds,
                       popsize = popsize,
                       generations = generations, cprob = cprob, cdist = cdist,
                       mprob = mprob, mdist=mdist,vectorized = FALSE)
    hyperparams<-NSGA$par[1,]
  }
  weight<-hyperparams[1:ncol(X)]
  params<-hyperparams[-c(1:ncol(X))]
  NDA<-try(ndr(X,min_evalue = params[1],min_communality=params[2],
               com_communalities = params[3],min_R = params[4],weight=weight,
               covar=FALSE,cor_method=cor_method,
               cor_type=cor_type,min_comm=min_comm,
               Gamma=Gamma,null_model_type=null_model_type,
               mod_mode=mod_mode,use_rotation=use_rotation,
               rotation=rotation),silent=TRUE)
  fits<-list()
  for (i in 1:ncol(Y))
  {
    data<-cbind(Y[,i],NDA$scores)
    colnames(data)[1]<-colnames(Y)[i]
    colnames(data)[-1]<-paste("NDA",1:NDA$factors,sep="")
    data<-as.data.frame(data)
    fits[[i]]<-stats::lm(str2lang(paste(colnames(data)[1],"~",
                                  gsub(", ","+",
                                       toString(colnames(data)[-1])))),data)
  }
  P<-list()
  P$Call<-cl
  P$fval<-cost(hyperparams)
  P$hyperparams<-hyperparams
  P$pareto<-pareto
  P$X<-X
  P$Y<-Y
  P$NDA<-NDA
  P$fits<-fits
  P$NDA_weight<-weight
  P$NDA_min_evalue<-params[1]
  P$NDA_min_communality<-params[2]
  P$NDA_com_communalities<-params[3]
  P$min_R <- params[4]
  if (optimize==TRUE){
    P$NSGA<-NSGA
  }
  P$fn<-"NDRLM"
  class(P)<-c("ndrlm","list")
  return(P)
}

