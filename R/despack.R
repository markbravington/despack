# This is package despack 

"CK_logalike_funs" <-
function( env, TMBOBJ, suffix=''){
r"--{
}--"

stopifnot( 
    env %is.an% 'environment',
    TMBOBJ %is.a% 'list',
    all( hasName( TMBOBJ, cq( fn, gr, par, report)))
  )

  lglk <- function( params, report.=FALSE){
    stopifnot( 
        is.numeric( params), 
        length( params)==length( par)
      )  
      val <- fn( params)
      if( report.){ # then store results for later
        list2env( report(), environment( sys.function()))
      }
    return( -val) # NB MINUS-SIGN!
    }
    
  Dlglk <- function( params){ 
    stopifnot( 
        is.numeric( params), 
        length( params)==length( par)
      )  
    return( -gr( params)) # NB MINUS-SIGN!
    }

  reportees <- function( params, want=c( 'probs', 'all')){
    stopifnot( 
        is.numeric( params), 
        length( params)==length( par)
      )  
      fn( params)
      # Make useful quantities available afterwards...
      stuff <- report()
      want <- match.arg( want)
      if( want=='probs'){
        # drop other stuff
        stuff <- stuff[ grepl( '^Pr_', names( stuff))]
      }
    return( stuff)
    }
    
  cenv <- new.env( parent=env)
  list2env( TMBOBJ[ cq( fn, gr, par, report)], cenv)
  environment( lglk) <- cenv
  environment( Dlglk) <- cenv
  environment( reportees) <- cenv
  rets <- returnList( lglk, Dlglk, reportees)
  names( rets) <- names( rets) %&% suffix
return( rets)
}


"finfo_onetype" <-
function( dsp, ncomp){
stopifnot( my.all.equal( 
    unname( dimseq( dsp))[-1], 
    unname( dimseq( ncomp))))
  DIM <- dim( dsp)
  dsp <- c( dsp) # strips all attributes
  npar <- DIM[1]
  dim( dsp) <- c( npar, length( dsp) %/% npar)
  rootio <- sqrt( c( ncomp)) * t( dsp)
  
  # covarset * params
  H <- 4 * crossprod( rootio)
return( H)
}


"get_Hbits" <-
function( 
  PARS_FOR_H, 
  all_probs_fun, Pargs= list(),
  numderiv_fun= mvbutils::numvbderiv_parallel,
  Dargs= list( eps=1e-6)
){
## This is a beautiful bit of code--- it's completely generic!
## Which means it will be incomprehensible :) But it *is* beautiful ;)
## Watch out for eps. With 1e-5 I got NAs for SPC ALB

  fcount <- 0
  problens <- all_probs <- is_prob <- is_num <- NULL # so <<- works
  sqrt_probs_plus <- function( pars){
      fcount <<- fcount + 1
      cat( sprintf( '%6i\r', fcount)); flush.console()

      # Allow other args specific to all_probs_fun
      # instead of just all_probs_fun( xpars, ...)
      all_probs <<- do.call( 'all_probs_fun', 
          c( list( pars), Pargs))

      if( is.null( problens)) {
        # Discard non-nums
        is_num <<- sapply( all_probs, is.numeric)
        all_probs <<- all_probs[ is_num]
        problens <<- lengths( all_probs)
        
        # Decide which are probs
        is_prob <<- startsWith( names( all_probs), 'Pr')
      } else { # drop non-numerics
        all_probs <<- all_probs[ is_num]
      }
      
      ans <- sqrt( unlist( all_probs[ is_prob], use.names=FALSE))
      
      # Anything else (eg lglk, popdyn vals): not sqrted
      ans <- c( ans, unlist( all_probs[ !is_prob], use.names=FALSE))
    return( ans)
    }

  infoid <- t( do.call( 'numderiv_fun', c( 
      list( sqrt_probs_plus, PARS_FOR_H), Dargs)))
 
  if( !all( is.finite( infoid))) {
stop( "NAs in numderiv; try changing something (eg eps)")
  }

  # Ensure all_probs is exact (for return, in case stuff gets set in envir)
  sqrt_probs_plus( PARS_FOR_H)

  # Need to pull out the right cols of the numderiv
  # This assumes that make_H1() is called in the same order as all_probs
  # which it is nowadays
  colind <- 0
  make_DSP1 <- function( Pr_name) {
      # produce offarray with same dims as Pr, except preceded by param-index
      Pr <- all_probs[[ Pr_name]]
      n_types <- length( Pr)
      dsp <- infoid[, colind + 1:n_types]
      colind <<- colind + n_types

      # Pre-extend dimensions. Should now handle non-offarrays
      # ie vectors, matrices, regular arrays. Should test that, really...
      dimatts <- attributes( Pr)
      if( is.null( dimatts)){
        names( dsp) <- names( PARS_FOR_H)
      } else {
        dimatts <- within( dimatts, {
          if( Pr %is.an% 'offarray'){
            offset <- c( 1, offset)
            dim <- c( PARAMS=NPAR, dim)
          } else if( !is.null( dimatts$dim)){
            dim <- c( NPAR, dim)
          } else { # pure vector
            names( Pr) <- PARAMS
          }

          if( !is.null( dimatts$dimnames)){
            dimnames <- c( list( 
                PARAMS=names( PARS_FOR_H)), # maybe NULL, that's fine
                dimnames)
          }
        })
        attributes( dsp) <- dimatts
      }
    return( dsp)
    }

  # Do 'em all...
  NPAR <- length( PARS_FOR_H)
  DSP <- FOR( names( all_probs)[ is_prob], make_DSP1( .))
  names( DSP) <- sub( 'Pr', 'DSP', names( DSP))

  if( !all( is_prob)){
    # then unpack remaining cols of infoid
    Dnonprobs <- FOR( names( all_probs)[ !is_prob], make_DSP1(.))
  } else {
    Dnonprobs <- NULL
  }

returnList( 
    Prkin=all_probs[ is_prob], 
    DSP, 
    Dnonprobs,
    PARS_FOR_H
  )
}

