\name{finfo_onetype}
\alias{finfo_onetype}
\title{Expected Fisher Info from all CK comparisons of given type}
\description{
This function computes the expected Fisher Information (ie negative of the Hessian) from a set of CK pairwise comparisons all of the same "type", but with many combinations of individual covariates. The "of the same type" means, for example, that this could be applied to cross-cohort half-sib comparisons between place A and place B where both samples have exact age measurements. A different "type" would be needed for, say, parent-offspring comparisons, or for half-sib comparisons where one animal has an inexact age, or...

The two ingredients are:

\itemize{
\item derivs of sqrts of the CK probs for all covariate combinations, with respect to parameters;
\item expected number of comparisons, for all the same covariate combinations.
}

The former comes from calling \code{\link{get_Hbits}} (qv). The number of comparisons is something you (partly) have to set up yourself, though there might be functions around to help.
\subsection{Number of comparisons}{
If the indices of \code{ncomp} are \code{(j1,k1,...1,j2,k2,...2)} where \code{1} pertains to the first sample and \code{2} to the second, then in general

\preformatted{
    ncomp(j1,k1,...,j2,k2,...2) =
      nsamp(j1,k1,...1) * nsamp(j2,k2,...2) * const(j1,......)
}

where \code{const()} is usually 1 (include all such comparisons) or 0 (omit these even if they could be done). Sometimes it is 0.5, to avoid double-counting when \code{j1==j2 & k1==k2 & ...}, etc.
}
\subsection{Note}{\code{finfo_onetype} only deals with \emph{one} type of CK comparison at a time (see its doco) whereas \code{\link{get_Hbits}} deals with all of them. Thus, if you have more than one type of CK comparison (say XHSPs as well as POPs), then you'll need multiple calls to \code{finfo_onetype} to get a Hessian from each type, which can then be added up. You can also add in Hessians stuff from (certain types of) non-CK data, such as age-composition samples and priors on parameters/latent variables.
}
}
\usage{finfo_onetype( dsp, ncomp)
}
\arguments{
\item{ dsp}{An (off)array with indices \code{(i,j,k,...)} where \code{i} refers to the i-th parameter, and \code{j,k,...} to all the covariates. \code{dsp[i,j,k,...]} is \code{D(sqrt(Pr_this_kintype|j,k,...))/D(params(i))}, evaluated at (some guess at) the true parameter values.}
\item{ ncomp}{An (off)array with indices \code{(j,k,...)}, as for \code{dsp}.}
}
\value{A square matrix with dimension \code{(n_params,n_params)}, where \code{n_params=length( dim( dsp)[1]))}.
}
\seealso{\code{\link{get_Hbits}}
}
\examples{
# See 'get_Hbits'
}
\keyword{misc
}
