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
# Last modified: December 2024                                                #
#-----------------------------------------------------------------------------#
## PRINT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ##
#' @export
print.ndrlm <- function(x, digits = getOption("digits"), ...) {
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (methods::is(x,"ndrlm")){
    Call<-x$Call
    fval<-x$fval
    pareto<-x$pareto
    X<-x$X
    Y<-x$Y
    NDA<-x$NDA
    fits<-x$fits
    NDA_weight<-x$NDA_weight
    NDA_min_evalue<-x$NDA_min_evalue
    NDA_min_communality<-x$NDA_min_communality
    NDA_com_communalities<-x$NDA_com_communalities
    min_R<-x$min_R
    optimized<-x$optimized
    if (optimized==TRUE){
      NSGA<-x$NSGA
    }
    extra_vars<-x$extra_vars
    if (extra_vars==TRUE){
      dircon_X<-x$dircon_X
    }
    fn<-x$fn
    cat("\nBrief summary of NDRLM:\n")
    cat("\nFunction call: ")
    print(Call)
    cat("\nNumber of independent variables: ",ncol(X))
    cat("\nNumber of dependent variables: ",ncol(Y))
    cat("\nNumber of latent-independent variables: ",ncol(NDA$factors))
    cat("\nNumber of dropped independent variables: ",sum((NDA$membership==0)))
    cat("\n\nSummary of dimensionality reduction\n")
    print.nda(NDA,digits = digits)

    cat("\n\nSummary of fitting\n")
    if (optimized==TRUE){
      cat("\nOptimized fittings\n")
    }else{
      cat("\nNon-optimized fittings\n")
    }
    cat("\nList of dependent variables: ",toString(colnames(Y)))
    cat("\nList of independent variables: ",toString(colnames(X)))
    cat("\nList of latent-independent variables: ",toString(colnames(NDA$scores)))
    if (extra_vars==TRUE){
      cat("\nList of non-groupped independent variables: ",toString(dircon_X))
    }

    for (i in 1:length(fits)){
      cat("\nFitting for variable ",colnames(fits[[i]]$model)[1])
      print(lm.beta::summary.lm.beta(lm.beta::lm.beta(fits[[i]])))
    }
  }else{
    print(x,...)
  }
}
