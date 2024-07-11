\name{CK_logalike_funs}
\alias{CK_logalike_funs}
\title{Make wrappers for TMB(O) functions}
\description{
Takes a TMB(O) object as returned by \code{MakeADFun}, and Returns a list containing three functions: \code{lglk}, \code{Dlglk}, and \code{reportees}. \code{lglk} can be used exactly as in the all-R example, and returns a \bold{positive} log-likelihood! \code{Dlglk} returns the gradient. \code{reportees} returns all the things from TMB REPORT() statements (as a list).

These wrappers are just convenient things that facilitate subsequent use in \code{despack} functions. You could easily enough write functions yourself to call \code{TMBOBJ$fn} and \code{TMBOBJ$report} and so on, and reorganize the output. But you don't have to, because I have...

Before running this, you have to manually call something like \code{boring_data_prep_<blah>()} to set up an environment containing all the data (which will be \code{env} here), and then \code{MakeAdFun} to create the "TMB object" (which will be \code{TMBOBJ} here).

Calling \code{lglk(<pars>,report.=TRUE) will create variables in }environment(lglk)' that can be accessed afterwards, just like what happens in the R{} version.

Calling \code{reportees} will by default just return variables whose name starts with "Pr_". You can make it also return other things, such as pop dyn quantities or expected kin-pairs, via the parameter \code{want="all"}.
}
\usage{CK_logalike_funs(env, TMBOBJ, suffix='')
}
\arguments{
\item{ env}{environment containing all necessary data}
\item{ TMBOBJ}{TMB object (an R{} list) returned by \code{MakeADFun}.}
\item{ suffix}{string to append to the names of each function--- so \code{suffix="_D"} would lead to a list with functions called \code{lglk_D}, \code{Dlglk_D}, and \code{reportees_D}.}
}
\value{A list of three functions.
}
\seealso{\code{\link{get_Hbits}}
}
\keyword{misc
}
