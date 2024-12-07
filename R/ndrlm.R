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
ndrlm<-function(Y,X,latents="in",dircon=FALSE,optimize=TRUE,cor_method=1,
                cor_type=1,min_comm=2,Gamma=1,
                null_model_type=4,mod_mode=1,use_rotation=FALSE,
                rotation="oblimin",pareto=FALSE,fit_weights=NULL,
                lower.bounds.x = c(rep(-100,ncol(X))),
                upper.bounds.x = c(rep(100,ncol(X))),
                lower.bounds.latentx = c(0,0,0,0),
                upper.bounds.latentx = c(0.6,0.6,0.6,0.3),
                lower.bounds.y = c(rep(-100,ncol(Y))),
                upper.bounds.y = c(rep(100,ncol(Y))),
                lower.bounds.latenty = c(0,0,0,0),
                upper.bounds.latenty = c(0.6,0.6,0.6,0.3),
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

  weight.X<-rep(1,ncol(X))
  weight.Y<-rep(1,ncol(Y))
  latent.X<-c(0,0,0,0)
  latent.Y<-c(0,0,0,0)

  if (("in" %in% latents)==FALSE){ # Pareto-optimiality can be found,
    pareto=FALSE         #  if there are only latent-independent variables
  }
  if (("none" %in% latents)==TRUE){ # If there are no latent variables,
    optimize=FALSE # there is no way to optimize fittings.
  }

  hyperparams=switch(
    latents,
    "in"=c(weight.X,latent.X),
    "out"=c(weight.Y,latent.Y),
    "both"=c(weight.X,latent.X,weight.Y,latent.Y),
    "none"=NULL
  )
  lower.bounds=switch(
    latents,
    "in"=c(lower.bounds.x,lower.bounds.latentx),
    "out"=c(lower.bounds.y,lower.bounds.latenty),
    "both"=c(lower.bounds.x,lower.bounds.latentx,
             lower.bounds.y,lower.bounds.latenty),
    "none"=NULL
  )
  upper.bounds=switch(
    latents,
    "in"=c(upper.bounds.x,lower.bounds.latentx),
    "out"=c(upper.bounds.y,lower.bounds.latenty),
    "both"=c(upper.bounds.x,upper.bounds.latentx,
             upper.bounds.y,upper.bounds.latenty),
    "none"=NULL
  )
  cost<-function(hyperparams){
    if ("in" %in% latents){
      weight.X<-hyperparams[1:ncol(X)]
      params.X<-hyperparams[-c(1:ncol(X))]
    }else{
      if ("out" %in% latents){
        weight.Y<-hyperparams[1:ncol(Y)]
        params.Y<-hyperparams[-c(1:ncol(Y))]
      }else{
        if ("both" %in% latents){
          weight.X<-hyperparams[1:ncol(X)]
          params.X<-hyperparams[(ncol(X)+1):(ncol(X)+4)]
          weight.Y<-hyperparams[(ncol(X)+5):(ncol(X)+4+ncol(Y))]
          params.Y<-hyperparams[(ncol(X)+5+ncol(Y)):(ncol(X)+4+ncol(Y)+4)]
        }
      }
    }
    if (latents %in% c("in","both")){ # For latent-independent variables
      NDA_in<-try(ndr(X,min_evalue = params.X[1],
                      min_communality=params.X[2],
                      com_communalities = params.X[3],
                      min_R = params.X[4],weight=weight.X,covar=FALSE,
                      cor_method=cor_method,
                      cor_type=cor_type,
                      min_comm=min_comm,
                      Gamma=Gamma,
                      null_model_type=null_model_type,
                      mod_mode=mod_mode,
                      use_rotation=use_rotation,
                      rotation=rotation),silent=TRUE)
    }

    if (latents %in% c("out","both")){ # For latent-dependent variables
      NDA_out<-try(ndr(Y,min_evalue = params.Y[1],
                       min_communality=params.Y[2],
                       com_communalities = params.Y[3],
                       min_R = params.Y[4],weight=weight.Y,covar=FALSE,
                       cor_method=cor_method,
                       cor_type=cor_type,
                       min_comm=min_comm,
                       Gamma=Gamma,
                       null_model_type=null_model_type,
                       mod_mode=mod_mode,
                       use_rotation=use_rotation,
                       rotation=rotation),silent=TRUE)
    }
    if (pareto==FALSE){
      res=0
    }else{
      res<-c(rep(0,ncol(Y)))
    }
    error<-FALSE
    if (latents %in% c("in","both")){
      if (inherits(NDA_in,"try-error")){
        error<-TRUE
        return(res)
      }
    }else{
      if (latents %in% c("out","both")){
        if (inherits(NDA_out,"try-error")){
          error<-TRUE
          return(res)
        }
      }
    }

    if (error==FALSE){
      extra_vars<-FALSE
      if (latents %in% c("in","both")){
        if ((dircon==TRUE)&&(sum(NDA_in$membership==0)>0)){
          extra_vars<-TRUE
          dropped_X<-X[,NDA_in$membership==0]
        }
      }else{
        if (latents %in% c("out","both")){
          if ((dircon==TRUE)&&(sum(NDA_out$membership==0)>0)){
            extra_vars<-TRUE
            dropped_Y<-Y[,NDA_out$membership==0]
          }
        }
      }
      dep<-Y
      if (latents %in% c("out","both")){
        if (extra_vars==TRUE){
          dep<-cbind(NDA_out$scores,dropped_Y)
          dep<-as.data.frame(dep)
          colnames(dep)<-c(paste("NDAout",1:NDA_out$factors,sep=""),
                           colnames(dropped_Y))
        }else{
          dep<-NDA_out$scores
          colnames(dep)<-paste("NDAout",1:NDA_out$factors,sep="")
        }
      }
      indep<-X
      if (latents %in% c("in","both")){
        if (extra_vars==TRUE){
          indep<-cbind(NDA_in$scores,dropped_X)
          indep<-as.data.frame(indep)
          colnames(indep)<-c(paste("NDAin",1:NDA_in$factors,sep=""),
                           colnames(dropped_X))
        }else{
          indep<-NDA_in$scores
          colnames(indep)<-paste("NDAin",1:NDA_in$factors,sep="")
        }
      }

      for (i in 1:ncol(dep))
      {
        data<-cbind(dep[,i],indep)
        colnames(data)[1]<-colnames(dep)[i]
        colnames(data)[-1]<-colnames(indep)
        data<-as.data.frame(data)
        fit<-stats::lm(str2lang(paste(colnames(data)[1],"~",
                                      gsub(", ","+",
                                           toString(colnames(data)[-1])))),data)

        res[i]<-stats::summary.lm(fit)$adj.r.squared
      }
      if (pareto==TRUE){
        return(res)
      }else{
        if (is.null(fit_weights)){
          return(mean(res))
        }else{
          return(stats::weighted.mean(res,fit_weights))
        }
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
  if ((optimize==TRUE)&&(latents %in% c("in","out","both"))){
    NSGA <- mco::nsga2(fn=costmin,idim=length(hyperparams),odim=ODIM,
                       lower.bounds = lower.bounds,
                       upper.bounds = upper.bounds,
                       popsize = popsize,
                       generations = generations, cprob = cprob, cdist = cdist,
                       mprob = mprob, mdist=mdist,vectorized = FALSE)
    hyperparams<-NSGA$par[1,]
  }

  if ("in" %in% latents){
    weight.X<-hyperparams[1:ncol(X)]
    params.X<-hyperparams[-c(1:ncol(X))]
  }else{
    if ("out" %in% latents){
      weight.Y<-hyperparams[1:ncol(Y)]
      params.Y<-hyperparams[-c(1:ncol(Y))]
    }else{
      if ("both" %in% latents){
        weight.X<-hyperparams[1:ncol(X)]
        params.X<-hyperparams[(ncol(X)+1):(ncol(X)+4)]
        weight.Y<-hyperparams[(ncol(X)+5):(ncol(X)+4+ncol(Y))]
        params.Y<-hyperparams[(ncol(X)+5+ncol(Y)):(ncol(X)+4+ncol(Y)+4)]
      }
    }
  }
  if (latents %in% c("in","both")){ # For latent-independent variables
    NDA_in<-try(ndr(X,min_evalue = params.X[1],
                    min_communality=params.X[2],
                    com_communalities = params.X[3],
                    min_R = params.X[4],weight=weight.X,covar=FALSE,
                    cor_method=cor_method,
                    cor_type=cor_type,
                    min_comm=min_comm,
                    Gamma=Gamma,
                    null_model_type=null_model_type,
                    mod_mode=mod_mode,
                    use_rotation=use_rotation,
                    rotation=rotation),silent=TRUE)
  }

  if (latents %in% c("out","both")){ # For latent-dependent variables
    NDA_out<-try(ndr(Y,min_evalue = params.Y[1],
                     min_communality=params.Y[2],
                     com_communalities = params.Y[3],
                     min_R = params.Y[4],weight=weight.Y,covar=FALSE,
                     cor_method=cor_method,
                     cor_type=cor_type,
                     min_comm=min_comm,
                     Gamma=Gamma,
                     null_model_type=null_model_type,
                     mod_mode=mod_mode,
                     use_rotation=use_rotation,
                     rotation=rotation),silent=TRUE)
  }
  fits<-list()
  extra_vars<-FALSE
  if (latents %in% c("in","both")){
    if ((dircon==TRUE)&&(sum(NDA_in$membership==0)>0)){
      extra_vars<-TRUE
      dropped_X<-X[,NDA_in$membership==0]
    }
  }else{
    if (latents %in% c("out","both")){
      if ((dircon==TRUE)&&(sum(NDA_out$membership==0)>0)){
        extra_vars<-TRUE
        dropped_Y<-Y[,NDA_out$membership==0]
      }
    }
  }
  dep<-Y
  if (latents %in% c("out","both")){
    if (extra_vars==TRUE){
      dep<-cbind(NDA_out$scores,dropped_Y)
      dep<-as.data.frame(dep)
      colnames(dep)<-c(paste("NDAout",1:NDA_out$factors,sep=""),
                       colnames(dropped_Y))
    }else{
      dep<-NDA_out$scores
      colnames(dep)<-paste("NDAout",1:NDA_out$factors,sep="")
    }
  }
  indep<-X
  if (latents %in% c("in","both")){
    if (extra_vars==TRUE){
      indep<-cbind(NDA_in$scores,dropped_X)
      indep<-as.data.frame(indep)
      colnames(indep)<-c(paste("NDAin",1:NDA_in$factors,sep=""),
                         colnames(dropped_X))
    }else{
      indep<-NDA_in$scores
      colnames(indep)<-paste("NDAin",1:NDA_in$factors,sep="")
    }
  }

  for (i in 1:ncol(dep))
  {
    data<-cbind(dep[,i],indep)
    colnames(data)[1]<-colnames(dep)[i]
    colnames(data)[-1]<-colnames(indep)
    data<-as.data.frame(data)
    fit<-stats::lm(str2lang(paste(colnames(data)[1],"~",
                                  gsub(", ","+",
                                       toString(colnames(data)[-1])))),data)

    fits[[i]]<-fit
  }



  P<-list()
  P$Call<-cl
  P$fval<-cost(hyperparams)
  P$hyperparams<-hyperparams
  P$pareto<-pareto
  P$X<-X
  P$Y<-Y
  P$latents<-latents
  if (latents %in% c("in","both")){
    P$NDAin<-NDA_in
    P$NDAin_weight<-weight.X
    P$NDAin_min_evalue<-params.X[1]
    P$NDAin_min_communality<-params.X[2]
    P$NDAin_com_communalities<-params.X[3]
    P$NDAin_min_R <- params.X[4]
  }
  if (latents %in% c("out","both")){
    P$NDAout<-NDA_out
    P$NDAout_weight<-weight.Y
    P$NDAout_min_evalue<-params.Y[1]
    P$NDAout_min_communality<-params.Y[2]
    P$NDAout_com_communalities<-params.Y[3]
    P$NDAout_min_R <- params.Y[4]
  }
  P$fits<-fits
  P$optimized<-optimize
  if (optimize==TRUE){
    P$NSGA<-NSGA
  }
  P$extra_vars<-extra_vars
  if (latents %in% c("in","both")){
    if (extra_vars==TRUE){
      P$dircon_X<-colnames(dropped_X)
    }
  }
  if (latents %in% c("out","both")){
    if (extra_vars==TRUE){
      P$dircon_Y<-colnames(dropped_Y)
    }
  }
  P$fn<-"NDRLM"
  class(P)<-c("ndrlm","list")
  return(P)
}

