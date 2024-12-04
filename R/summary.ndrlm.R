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
# Last modified: February 2024                                                #
#-----------------------------------------------------------------------------#
## SUMMARY FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ##
#' @export
summary.ndrlm <- function(object,  digits =  getOption("digits"), ...) {
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (methods::is(object,"ndrlm")){
    Call<-object$Call
    fval<-object$fval
    pareto<-object$pareto
    X<-object$X
    Y<-object$Y
    NDA<-object$NDA
    fits<-object$fits
    NDA_weight<-object$NDA_weight
    NDA_min_evalue<-object$NDA_min_evalue
    NDA_min_communality<-object$NDA_min_communality
    NDA_com_communalities<-object$NDA_com_communalities
    min_R<-object$min_R
    NSGA<-object$NSGA
    fn<-object$fn
    if (!is.null(object$NSGA)){
      results<-list(Call=Call,
                    fval=fval,
                    pareto=pareto,
                    X = X,
                    Y = Y,
                    NDA = NDA,
                    fits = fits,
                    NDA_weight<-NDA_weight,
                    NDA_min_evalue<-NDA_min_evalue,
                    NDA_min_communality<-NDA_min_communality,
                    NDA_com_communalities<-NDA_com_communalities,
                    min_R<-min_R,
                    NSGA<-NSGA,
                    fn<-fn)
    }else{
      results<-list(Call=Call,
                    fval=fval,
                    pareto=pareto,
                    X = X,
                    Y = Y,
                    NDA = NDA,
                    fits = fits,
                    NDA_weight<-NDA_weight,
                    NDA_min_evalue<-NDA_min_evalue,
                    NDA_min_communality<-NDA_min_communality,
                    NDA_com_communalities<-NDA_com_communalities,
                    min_R<-min_R,
                    fn<-fn)
    }
    return(results)
    print.ndrlm(object)
  }else{
    summary(object,...)
  }
}
