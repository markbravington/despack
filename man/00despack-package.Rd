\name{despack-package}
\alias{despack}
\alias{despack-package}
\docType{package}
\title{How to use the despack package}
\description{
This package provides helper functions for CKMR design. To use it fully, you need:

\itemize{
\item working code for a CKMR model, either in R{} or TMB(O), that can calculate all the kinship probabilities;
\item guessestimates of true parameter values;
\item explicit "designs" ie proposals for sample sizes, disaggregated to whatever level would be used in the "real data" (eg by year and age and sex)
\item a clear definition of some quantity whose precision you're interested in (eg a total abundance in some specified year; a total biomass; the natural mortality rate; replacement yield in some year; ...). This should be a scalar function that takes all the pop dyn parameters as its argument.
}

Note that you do \emph{not} need any datasets, whether real or simulated!

There are currently 4 functions in the package, but \emph{don't} start by looking at their helpfiles! Look at the WORKFLOW below, then have a go at the script/vignette which is \emph{not} currently part of the package. The main "documentation" is really a script, with accompanying functions and data, that doesn't currently live inside this package (I'm not sure what the best arrangement is: maybe a "demo package" that goes alongside this one, or maybe "just" as a vignette, but it's complicated). Of necessity, CK examples tend to be quite application-specific, whereas the functions in this package are deliberately generic.
}
\section{Workflow}{
To calculate the likely SE of your quantity-of-interest, do this (assuming your model is in TMB(O) rather than R{}; see later for all-R comments):

\itemize{
\item \code{TMB(O)::MakeADFun} to set up data ranges and constants in the model (but no data required).
\item \code{\link{CK_logalike_funs}} to make convenient wrappers for calling your TMB(O) function, and for returning REPORTed quantities, such as CK probability arrays and population dynamics numbers;
\item \code{\link{get_Hbits}} to compute all the expected Fisher Information matrices ("negative Hessians") from doing one single comparison between two samples with specified covariates. Results are stored in compressed form.
\item \code{\link{finfo_onetype}} for each \emph{type} of CK comparison in turn (eg POPs, XMHSPs, POPs where one sample only has approximate age, ...). Here you will need to supply a crucial new argument, derived from the design itself: an array holding the relevant \emph{number} of comparisons disaggregated by covariates (eg an \code{offarray} called \code{n_comp_MOP_BYA}). Result is the expected Fisher Info for all comparisons of that \emph{type}.
\item Add up all the \code{\link{finfo_onetype}} results (there might only be one type, if your application is very simple) to get the overall Fisher Info.
\item Invert that to get the parameter covariance matrix
\item Double-dot that with the numerical derivative of your quantity-of-interest function (ie apply "the delta method"), to get the expected approximate variance of your thing.
}

A simple example is given in the vignette (which I haven't yet included in the package... there's a separate script).

There is also a schematic example (ie it shows plausible R{} code, but doesn't actually work) in the help for \code{\link{get_Hbits}}. Also, there is individual documentation for each of \code{\link{CK_logalike_funs}}, \code{\link{get_Hbits}}, and \code{\link{finfo_onetype}}, plus for the helper function \code{mvbutils::numvbderiv} and \emph{especially} its parallel version \code{mvbutils::numvbderiv_parallel}, but for goodness sake don't \emph{start} by reading their helpfiles. Follow the vignette/script instead, until and unless you get stuck.
\subsection{R only design}{If your CK model is entirely in R--- presumably using \code{offarray} otherwise you will go mad with index confusion--- then you obviously don't need the \code{MakeADFun} step, and the second step can be "hand-written" easily; again, the vignette/script is the place to look, for the function \code{reportees_Ronly}. Following steps are identical. Because there's no need to fit the model to actual data, an all-R version is quite feasible even for pretty complex designs--- though if you ever turn the design into reality, you might need to write a TMB(O) version anyway to be able to actually fit it--- though package \pkg{RTMB} may now obviate C-TMB code altogether. But you'll have several years before you need to do that ;)
}
}
\section{Efficiently investigating different designs}{
Once you have gotten all the steps working "manually" for your application, you will want to consider various possible designs and see how they influence Standard Errors / CVs. Since there are a huge number of individual sample sizes that potentially \emph{could} be controlled, you can't really do that manually. For serious use, you will instead want to write your own "driver" function that groups together several of these steps; \code{varcalcs_delfi_A} (which accompanies the script) is an example. The key part is a "design generator", which:

\itemize{
\item takes a small number of control parameters (total sample sizes by category, trends in sample size by year, how to apply selectivity, ...; it all depends on your application);
\item generates disaggregated sample sizes, ie a complete design specification;
}

Once you have generated the detailed design, the rest of your driver function needs to:

\itemize{
\item runs the design thru \code{\link{finfo_onetype}} (also using results of one previous call to \code{\link{get_Hbits}}), once for each type of CK comparison (eg MOPs, XSHPs);
\item sums those Hessians, and inverts the result to get the parameter covariance matrix;
\item applies the "delta method" to that covariance matrix, using numerical derivatives of your quantity-of-interest WRTO the parameters, to calculate the variance of your Thing.
}

You can also incorporate Hessian components from (some types of) non-CKMR data, such as age-at-length compositional samples, and priors on random effects and parameters--- probably not CPUE though :) I've done that myself in several earlier examples. It's really up to you...
}
\section{Not for now}{"Optimal" design subject to arbitrary constraints, using a "driver" function as per last section. (It actually works, but will take me some time, ie at least months, to port over into this package.)
}
\keyword{misc}
