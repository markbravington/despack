\name{get_Hbits}
\alias{get_Hbits}
\title{Preparation for CK design Hessian calculation}
\description{
\code{get_Hbits} prepares the ingredients needed to compute the expected Fisher Information ("negative Hessian") from any kinship comparison between two samples, depending on the covariate values of the samples and the type of kinship being considered. To calculate the overall Hessian, the results then need to be combined with proposed sample size information and added up, which is done (partly) by \code{\link{finfo_onetype}} (qv).

You can also use it, at the same time or later, to conveniently calculate other numerical derivatives of your CK model, eg of a log-likelihood (if you have real or simulated data) or of "quantities of interest" (for use in Delta-method later on). If your CK stuff is at all complicated, then it is generally cheaper to do a bunch of numerical derivs all at once. However, that's not compulsory. My most general pattern is:

\itemize{
\item 1. call \code{get_Hbits} once just to get the prob derivs;
\item 2. pick a sample-size scenario, and call \code{\link{finfo_onetype}} (as many times as required, once per different prob array) to get the overall Hessian for that scenario;
\item 3. Some numderiv stuff for derivs of quantities-of-interest (for which I might again use \code{get_Hbits});
\item 4. Delta-method to combine #3 with #2.
}

I like the flexibility of being able to deal with quantities-of-interest after-the-fact; they do not affect the parameter Hessian. Also, they can often be calculated very quickly, without needing to compute all CK prob arrays. However, if you are sure you know what your QofIs are, then you \emph{can} merge steps 1 & 3 into a single initial call to \code{get_Hbits}, which will save some time if your model is big and/or has lots of parameters.

In more detail: your CK code should compute at least one probability array, more if there are two types of kin or if some animals have qualitatively different covariates from others. \code{get_Hbits} automatically computes the derivatives of the square-roots of those probabilities. If your code also returns other (must be numerical) quantities, \code{get_Hbits} calculates the derivs of those \emph{untransformed} things (i.e. no square root).
}
\usage{
get_Hbits(
  PARS_FOR_H,
  all_probs_fun, Pargs= list(),
  numderiv_fun= mvbutils::numvbderiv_parallel,
  Dargs= list( eps=1e-6)
)
}
\arguments{
\item{ PARS_FOR_H}{"true" parameter values}
\item{ all_probs_fun}{function taking parameter vector as first argument, and returning a list of all CK prob arrays (or a list of any numeric things, actually; non-numeric elements are ignored). Prob arrays must have names starting "Pr". If you wrote your CK stuff in TMB(O), then \code{all_probs_fun} might just be \code{CK_logalike_from_TMBO_obj(...)$reportees}. Can optionally have other args...}
\item{ Pargs}{... which you set via this argument, eg \code{list(want=TRUE)} for functions obtained from \code{CK_logalike_from_TMB_obj}.}
\item{ numderiv_fun}{function that computes numerical derivatives of its first arg WRTO its second arg (a vector of parameters). The parameter dimension should be the last one in the result. Default should be OK, but don't blame me if it isn't; check its helpfile. See \bold{Details}.}
\item{ Dargs}{optional list of extra args for \code{numderiv_fun}; for the default \code{numderiv_fun}, this might be \code{list( eps=1e-7)} to change the step, or \code{list( PARALLEL=FALSE)} if you have not set up a parallel cluster (which you should; \code{get_Hbits} can be slow).}
}
\value{
A list with components \code{DSP}, \code{Dnonprob}, \code{PARS_FOR_H}, and \code{Prkin}. The latter is the actual CK probabilities at \code{PARS_FOR_H}; it's a list, because there might be more than one type of CK probability. \code{DSP} is also a list, for the same reason; it's the derivs of the square-roots of the probabilities. \code{Dnonprob} is also a list (possibly empty) holdings the derivs of any non-probability returnees from \code{all_probs_fun}, such as an actual log-likelihood or some popdyn quantities; square-roots are \emph{not} applied to such things.
You can't really use the key return value, \code{DSP}, directly; it only makes sense in future calls to \code{\link{finfo_onetype}}.
}
\details{
In this situation, where the dimension of the output (the nubmer of covariate-combinations in the CK probabilities) is large, numerical differentiation is almost as fast as Automatic Differentiation, and is completely general. The nice-sounding idea of using forward-mode AD apparently can't even be done in TMB! So, relax and just use numeric derivs.

The numerical derivatives don't have to be particularly accurate for this application. Currently, \code{get_Hbits} uses my \code{numvbderiv} which is very simple--- and fast, if you use the parallel version, which is the default. However, it not as accurate or robust as if you used a special-purpose R{} package. Note that it will be called via \code{numderiv_fun( fun_to_diff, param_vals, <Dargs>)} and \code{fun_to_diff} must be a \emph{function} taking a parameter vector as its first argument. Thus, if you wanted to use \code{stats::numericDeriv} (which I do not recommend), you'd have to write a wrapper for it, because it expects an \code{expression} argument not a function. Yawn. Boring!
}
\seealso{\code{\link{finfo_onetype}}, \code{CK_logalike_from_TMB_obj}
}
\examples{
\dontrun{
if( require( TMBO))\{ # offsets like offarray; better debugging
  compile( 'myCKex')
  dyn.load( dynlib( 'myCKex'))
  tmbob <- with( env, # so it knows about all the variables
    MakeADFun(
      data= returnList(
          Amat,
          n_MOP_BYA= array(NA),
          n_comp_MOP_BYA= array( NA)
        ),
      parameters= list( log_Nfad_ystart= starto[1], RoI= starto[2]),
      ranges= TMBO_ranges( years, Yad_range, Aad_range, Bju_range),
      DLL="lglk_POP_ideal_mammal",
      silent=TRUE
    ))
  tmbeq <- CK_logalike_funs( tmbob)
  Hbits <- get_Hbits( trupars, tmbeq$reportees)
  design_n_comp_MOP_BYA <- "something or other"
  H_MOPonly <- finfo_onetype(
      Hbits$DSP$DSP_MOP_BYA,
      design_n_comp_MOP_BYA)
  Vpar_MOPonly <- solve( H_MOPonly) # expected parameter covariances
\}
}
}
\keyword{misc
}
