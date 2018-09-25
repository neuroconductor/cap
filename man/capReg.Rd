\name{capReg}
\alias{capReg}
\title{ Covariate Assisted Principal Regression for Covariance Matrix Outcomes
}
\description{ This function identifies the first \eqn{k} projection directions that satisfies the log-linear model assumption.
}
\usage{
capReg(Y, X, nD = 1, method = c("CAP", "CAP-C"), CAP.OC = FALSE, 
  max.itr = 1000, tol = 1e-04, trace = FALSE, score.return = TRUE, 
  gamma0.mat = NULL, ninitial = NULL)
}
\arguments{
  \item{Y}{ a data list of length \eqn{n}. Each list element is a \eqn{T\times p} matrix, the data matrix of \eqn{T} observations from \eqn{p} features.
}
  \item{X}{ a \eqn{n\times q} data matrix, the covariate matrix of \eqn{n} subjects with \eqn{q-1} predictors. The first column is all ones.
}
  \item{nD}{ an integer, the number of directions to be identified. Default is 1.
}
  \item{method}{ a character of optimization method. \code{method = "CAP"} considers a weighted L2-norm on the \eqn{\gamma} vector and solve for the optimizer by block coordinated descent; \code{method = "CAP-C"} assumes the complete common principal component assumption which identifies the common principal component first and then searches for the optimal PC.
}
  \item{CAP.OC}{ a logic variable. Whether the orthogonal constraint is imposed when identifying higher-order PCs. When \code{method = "CAP-C"}, this is ignored. Default is \code{FALSE}.
}
  \item{max.itr}{ an integer, the maximum number of iterations.
}
  \item{tol}{ a numeric value of convergence tolerance.
}
  \item{trace}{ a logic variable. Whether the solution path is reported. Default is \code{FALSE}.
}
  \item{score.return}{ a logic variable. Whether the log-variance in the transformed space is reported. Default is \code{TRUE}.
}
  \item{gamma0.mat}{ a data matrix, the initial value of \eqn{\gamma}. Default is \code{NULL}, and initial value is randomly chosen.
}
  \item{ninitial}{ an integer, the number of different initial value is tested. When it is greater than 1, multiple initial values will be tested, and the one yields the minimum objective function will be reported. Default is \code{NULL}.
}
}
\details{ Considering \eqn{y_{it}} are \eqn{p}-dimensional independent and identically distributed random samples from a multivariate normal distribution with mean zero and covariance matrix \eqn{\Sigma_{i}}. We assume there exits a \eqn{p}-dimensional vector \eqn{\gamma} such that \eqn{z_{it}:=\gamma'y_{it}} satisfies the multiplicative heteroscedasticity:
  \deqn{\log(\mathrm{Var}(z_{it}))=\log(\gamma'\Sigma_{i}\gamma)=\beta_{0}+x_{i}'\beta_{1}},
where \eqn{x_{i}} contains explanatory variables of subject \eqn{i}, and \eqn{\beta_{0}} and \eqn{\beta_{1}} are model coefficients.

Parameters \eqn{\gamma} and \eqn{\beta=(\beta_{0},\beta_{1}')'} are study of interest, and we propose to estimate them by maximizing the likelihood function,
  \deqn{\ell(\beta,\gamma)=-\frac{1}{2}\sum_{i=1}^{n}T_{i}(x_{i}'\beta)-\frac{1}{2}\sum_{i=1}^{n}\exp(-x_{i}'\beta)\gamma'S_{i}\gamma,}
where \eqn{S_{i}=\sum_{t=1}^{T_{i}}y_{it}y_{it}'}. To estimate \eqn{\gamma}, we impose the following constraint
  \deqn{\gamma' H\gamma=1,}
where \eqn{H} is a positive definite matrix. In this study, we consider the choice that
  \deqn{H=\bar{\Sigma}, \quad \bar{\Sigma}=\frac{1}{n}\sum_{i=1}^{n}\frac{1}{T_{i}}S_{i}.}

For higher order projecting directions, an orthogonal constraint is imposed as well.
}
\value{ When \code{method = "CAP"},
  \item{gamma}{the estimate of \eqn{\gamma} vectors, which is a \eqn{p\times nD} matrix.}
  \item{beta}{the estimate of \eqn{\beta} for each projecting direction, which is a \eqn{q\times nD} matrix, where \eqn{q-1} is the number of explanatory variables.}
  \item{orthogonality}{an ad hoc checking of the orthogonality between \eqn{\gamma} vectors.}
  \item{DfD}{output of both average (geometric mean) and individual level of ``deviation from diagonality''.}
  \item{score}{an output when \code{score.return = TRUE}. A \eqn{n\times nD} matrix of \eqn{\log(\hat{\gamma}'S_{i}\hat{\gamma})} value.}

  When \code{method = "CAP-C"},
    \item{gamma}{the estimate of \eqn{\gamma} vectors, which is a \eqn{p\times nD} matrix.}
  \item{beta}{the estimate of \eqn{\beta} for each projecting direction, which is a \eqn{q\times nD} matrix, where \eqn{q-1} is the number of explanatory variables.}
  \item{orthogonality}{an ad hoc checking of the orthogonality between \eqn{\gamma} vectors.}
  \item{PC.idx}{a vector of length \code{nD}, the order index of identified \eqn{\gamma} vectors among all the common principal components.}
  \item{aPC.idx}{the order index of all the principal components that satisfy the log-linear model and the eigenvalue condition.}
  \item{minmax}{a logic output, whether the identified \eqn{\gamma} vectors are estimated from the minmax approach. If \code{FALSE}, indicating the eigenvalue condition is not satisfied for any principal component.}
  \item{score}{an output when \code{score.return = TRUE}. A \eqn{n\times nD} matrix of \eqn{\log(\hat{\gamma}'S_{i}\hat{\gamma})} value.}
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

# method = "CAP"
# without orthogonal constraint
re1<-capReg(Y,X,nD=2,method=c("CAP"),CAP.OC=FALSE)
# with orthogonal constraint
re2<-capReg(Y,X,nD=2,method=c("CAP"),CAP.OC=TRUE)

# method = "CAP-C"
re3<-capReg(Y,X,nD=2,method=c("CAP-C"))
#############################################
}
\keyword{ models }
