\name{ndrlm}
\alias{ndrlm}
\title{Genearlized Network-based Dimensionality Reduction and Regression (GNDR)
}
\description{
The main function of Generalized Network-based Dimensionality Reduction and Regression (GNDR) for supervised learning.
}
\usage{
ndrlm(Y,X,optimize=TRUE,cor_method=1,cor_type=1,min_comm=2,Gamma=1,null_model_type=4,
mod_mode=6,use_rotation=FALSE,rotation="oblimin",pareto=TRUE,out_weights=rep(1,ncol(Y)),
lower.bounds = c(rep(-100,ncol(X)),0,0,0,0),
upper.bounds = c(rep(100,ncol(X)),0.6,0.6,0.6,0.3),
popsize = 20, generations = 30, cprob = 0.7, cdist = 5, mprob = 0.2, mdist=10, seed=NULL)

}
\arguments{
  \item{Y}{A numeric data frame of output variables}
  \item{X}{A numeric data frame of input variables}
  \item{optimize}{Optimization of fittings (default=TRUE)}
  \item{cor_method}{Correlation method (optional).
'1' Pearson's correlation (default),
'2' Spearman's correlation, '3' Kendall's correlation,
'4' Distance correlation}
  \item{cor_type}{Correlation type (optional).
'1' Bivariate correlation (default),
'2' partial correlation, '3' semi-partial correlation}
  \item{min_comm}{Minimal number of indicators per community (default: 2).}
  \item{Gamma}{Gamma parameter in multiresolution null modell (default: 1).}
  \item{null_model_type}{'1' Differential Newmann-Grivan's null model, '2' The null model is the mean of square correlations between indicators, '3' The null model is the specified minimal square correlation, '4' Newmann-Grivan's modell (default)}
  \item{mod_mode}{Community-based modularity calculation mode: '1' Louvain modularity, '2' Fast-greedy modularity, '3' Leading Eigen modularity, '4' Infomap modularity, '5' Walktrap modularity, '6' Leiden modularity (default)}
  \item{use_rotation}{FALSE no rotation (default), TRUE the rotation is used.}
  \item{rotation}{"none", "varimax", "quartimax", "promax", "oblimin", "simplimax", and "cluster" are possible rotations/transformations of the solution. "oblimin" is the default, if use_rotation is TRUE.}
  \item{pareto}{in the case of multiple objectives TRUE (default value) provides pareto-optimal solution, while FALSE provides weighted mean of objective functions (see out_weights)}
  \item{out_weights}{weights of fitting the output variables (weights of means of objectives)}
  \item{lower.bounds}{Lower bounds of weights and hyper parameters of GNDA (last 4 elements).}
  \item{upper.bounds}{Upper bounds of weights and hyper parameters of GNDA (last 4 elements).}
  \item{popsize}{size of population of NSGA-II for fitting betas (default=20)}
  \item{generations}{number of generations to breed of NSGA-II for fitting betas (default=30)}
  \item{cprob}{crossover probability of NSGA-II for fitting betas (default=0.7)}
  \item{cdist}{crossover distribution index of NSGA-II for fitting betas (default=5)}
  \item{mprob}{mutation probability of NSGA-II for fitting betas (default=0.2)}
  \item{mdist}{mutation distribution index of NSGA-II for fitting betas (default=10)}
  \item{seed}{default seed value (default=NULL, no seed)}
}
\details{
  NDRLM is a variable fitting with feature selection based on the tunes of GNDA method with NSGA-II algorithm for parameter fittings.
}
\value{
\item{Call}{Callback function}
\item{fval}{Objective function for fitting}
\item{hyperparams}{optimized hyperparameters}
\item{pareto}{in the case of multiple objectives TRUE (default value) provides pareto-optimal solution, while FALSE provides weighted mean of objective functions (see out_weights)}
\item{Y}{A numeric data frame of output variables}
\item{X}{A numeric data frame of input variables}
\item{NDA}{GNDA object, which is the result of model reduction and features selection}
\item{fits}{List of linear regrassion models}
\item{NDA_weight}{Weights of input variables (used in \code{\link{ndr}})}
\item{NDA_min_evalue}{Optimized minimal eigenvector centrality value (used in \code{\link{ndr}})}
\item{NDA_min_communality}{Optimized minimal communality value of indicators (used in \code{\link{ndr}})}
\item{NDA_com_communalities}{Optimized
minimal common communalities (used in \code{\link{ndr}})}
\item{NDA_min_R}{Optimized
minimal square correlation between indicators (used in \code{\link{ndr}})}
\item{NSGA}{Outpot structure of NSGA-II optimization (list), if the optimization value is true (see in \code{\link[mco:nsga2]{mco::nsga2}})}
\item{fn}{Function (regression) name: NDLM}
}

\references{
Kosztyan, Z. T., Kurbucz, M. T., & Katona, A. I. (2022). Network-based dimensionality reduction of high-dimensional, low-sample-size datasets. Knowledge-Based Systems, 109180. doi:10.1016/j.knosys.2022.109180
}

\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kosztyan.zsolt@gtk.uni-pannon.hu
}

\seealso{\code{\link{ndr}}, \code{\link{plot}}, \code{\link{summary}}, \code{\link[mco:nsga2]{mco::nsga2}}.}

\examples{

# Using GNDA + LM with curve fitting

X<-freeny.x
Y<-freeny.y
NDRLM<-ndrlm(Y,X)
summary(NDRLM$fits[[1]])

# Using GNDA + LM without curve fitting

NDRLM<-ndrlm(Y,X,optimize=FALSE)
summary(NDRLM$fits[[1]])

}

\keyword{multivariate}
\keyword{nonparametric}
\concept{dimensionality}
\concept{reduction}
