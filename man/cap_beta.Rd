\name{cap_beta}
\alias{cap_beta}
\title{ Inference of model coefficients
}
\description{ This function performs inference on the model coefficient \eqn{\beta}.
}
\usage{
cap_beta(Y, X, gamma = NULL, beta = NULL, method = c("asmp", "LLR"), 
    boot = FALSE, sims = 1000, boot.ci.type = c("bca", "perc"), 
    conf.level = 0.95, verbose = TRUE)
}
\arguments{
  \item{Y}{ a data list of length \eqn{n}. Each list element is a \eqn{T\times p} matrix, the data matrix of \eqn{T} observations from \eqn{p} features.
}
  \item{X}{ a \eqn{n\times q} data matrix, the covariate matrix of \eqn{n} subjects with \eqn{q-1} predictors. The first column is all ones.
}
  \item{gamma}{ a \eqn{p}-dimensional vector, the projecting direction \eqn{\gamma}. Default is \code{NULL}. If \code{gamma = NULL}, an error warning will be returned.
}
  \item{beta}{ a \eqn{q}-dimensional vector, the model coefficient \eqn{\beta}. Default is \code{NULL}. If \code{beta = NULL}, when \code{boot = FALSE}, \eqn{\beta} will be estimated using the provided \eqn{\gamma}.
}
  \item{method}{ a character of inference method. If \code{method = "asmp"}, the inference is made based on the asymptotic variance; if \code{method = "LLR"}, the likelihood ratio test is conducted. When \code{boot = TRUE}, this argument is ignored.
}
  \item{boot}{ a logic variable, whether bootstrap inference is performed.
}
  \item{sims}{ a numeric value, the number of bootstrap iterations will be performed.
}
  \item{boot.ci.type}{ a character of the way of calculating bootstrap confidence interval. If \code{boot.ci.type = "bca"}, the bias corrected confidence interval is returned; if \code{boot.ci.type = "perc"}, the percentile confidence interval is returned.
}
  \item{conf.level}{ a numeric value, the designated significance level. Default is \eqn{0.95}.
}
  \item{verbose}{ a logic variable, whether the bootstrap procedure is printed. Default is \code{TRUE}.
}
}
\details{ Considering \eqn{y_{it}} are \eqn{p}-dimensional independent and identically distributed random samples from a multivariate normal distribution with mean zero and covariance matrix \eqn{\Sigma_{i}}. We assume there exits a \eqn{p}-dimensional vector \eqn{\gamma} such that \eqn{z_{it}:=\gamma'y_{it}} satisfies the multiplicative heteroscedasticity:
  \deqn{\log(\mathrm{Var}(z_{it}))=\log(\gamma'\Sigma_{i}\gamma)=\beta_{0}+x_{i}'\beta_{1},}
where \eqn{x_{i}} contains explanatory variables of subject \eqn{i}, and \eqn{\beta_{0}} and \eqn{\beta_{1}} are model coefficients.

The \eqn{\beta} coefficient is estimated by maximizing the likelihood function. The asymptotic variance is obtained based on maximum likelihood estimator theory.
}
\value{ When \code{method = "asmp"}, the output is a \eqn{q \times 6} data frame containing the estimate of \eqn{\beta} coefficient, the asymptotic standard error, the test statistic, the \eqn{p}-value, and the lower and upper bound of the confidence interval.

    When \code{method = "LLR"}, the output is a \eqn{q \times 3} data frame containing the estimate of \eqn{\beta} coefficient, the test statistic, and the \eqn{p}-value.

    When \code{boot = TRUE},
    \item{Inference}{point estimate of the \eqn{\beta} coefficient, as well as the corresponding standard error, test statistic, \eqn{p}-value, and the lower and upper bound of the confidence interval.}
    \item{beta.boot}{the estimate of the \eqn{\beta} coefficient in each iteration.}
}
\references{Zhao et al. (2018) \emph{Covariate Assisted Principal Regression for Covariance Matrix Outcomes} <doi:10.1101/425033>
}
\author{
Yi Zhao, Johns Hopkins University, <zhaoyi1026@gmail.com> 

Bingkai Wang, Johns Hopkins University, <bwang51@jhmi.edu>

Stewart Mostofsky, Johns Hopkins University, <mostofsky@kennedykrieger.org> 

Brian Caffo, Johns Hopkins University, <bcaffo@gmail.com> 

Xi Luo, Brown University, <xi.rossi.luo@gmail.com> 
}
\examples{

#############################################
data(env.example)
X<-get("X",env.example)
Y<-get("Y",env.example)
Phi<-get("Phi",env.example)

# asymptotic variance
re1<-cap_beta(Y,X,gamma=Phi[,2],method=c("asmp"),boot=FALSE)

# likelihood ratio test
re2<-cap_beta(Y,X,gamma=Phi[,2],method=c("LLR"),boot=FALSE)

# bootstrap confidence interval
\donttest{
re3<-cap_beta(Y,X,gamma=Phi[,2],boot=TRUE,sims=500,verbose=FALSE)
}
#############################################
}
\keyword{ models }
