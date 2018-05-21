\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{emdbook News}
\encoding{UTF-8}

\section{Changes in version 1.3.9}{
  \itemize{
    \item add imports for R CMD check
    \item add labels to curve3d with wireframe plot
    \item change rgl from Depends: to Suggests: (helps MacOS users when
   rgl installation is troublesome)
   }
}

\section{Changes in version 1.3.8}{
   \itemize{
     \item CRAN compatibility tweaks; add link to goby data; update link
   to damselfish data
   }
}

\section{Changes in version 1.3.7}{
   \itemize{
     \item CRAN compatibility tweaks
   }
}
\section{Changes in version 1.3.5}{
  \itemize{
    \item CRAN compatibility tweaks
  }
}

\section{Changes in version 1.3.4}{
  \itemize{
    \item dbetabinom fix (reported by Christian Robert)
    \item more CRAN compatibility tweaks (versions 1.3.2 - 1.3.3)
    \item add 'lims' argument to HPDregionplot
    \item add plyr dependency, use_plyr option in curve3d
  }
}
\section{Changes in version 1.3.2.1}{
     \itemize{
       \item remove EISPACK=TRUE in dmvnorm
     }
}

\section{Changes in version 1.3.2}{
     \itemize{
       \item small fix for NAs in lambertW
       \item tweaks to avoid NOTES in CRAN check, cleanup
     }
}
\section{Changes in version 1.3.1}{
    \itemize{
      \item asymptotic extension to lambertW (remove possibility of
    vectorized branch specification)
      \item update NEWS to NEWS.Rd
  }}

\section{Changes in version 1.3}{
    \itemize{
      \item minor bugfixes for HPDregionplot 
      \item tweaks to pass checks on R 2.13
  }}

  \section{Changes in version 1.2.2.1 (24/05/2010)}{
    \itemize{
      \item add @data to calcslice
  }}
  
  \section{Changes in version 1.2.2 (17/03/2010)}{
    \itemize{
      \item add (trivial) NAMESPACE
      \item add contour3d
      \item add xlim/ylim arguments to curve3d
  }}
  
  \section{Changes in version 1.2.1 (11/02/2010)}{
    \itemize{
      \item remove adapt from install.packages()
      \item add warning about adapt
      \item move CITATION file
  }}

  \section{Changes in version 1.2 (17/08/2009)}{
    \itemize{
      \item change is.na() test in lambertW to is.finite()
      \item added ggplot2 to install.packages(), streamline
      \item REMOVED "traceplot" code: nearly identical to xyplot.mcmc etc. in coda package
  }}

  \section{Changes in version 1.1.1.1 (16/10/2008)}{
    \itemize{ 
      \item add firfec0.bug to BUGS
      \item reverse date order of NEWS
  }}

  \section{Changes in version 1.1.1}{
    \itemize{
      \item tweak to curve3d from Robert Burns
      \item moved ChangeLog to NEWS
      \item added axis.scinot
  }}
			
  \section{Changes in version 1.1.0}{
    \itemize{
      \item argh, another tweak to get.emdbook.packages
      \item tweaked deltavar to be slightly more general
  }}

  \section{Changes in version 1.0.9}{
    \itemize{
      \item changed "Wrasse" to "Damsel[fish]" throughout (oops)
      \item tweaked get.emdbook.packages (do.call rather than sapply)
      \item tweaked LambertW NA behavior (still not optimal)
}}
\section{Changes in version 1.0.8.1}{
  \itemize{
    \item ANOTHER big bug in dmvnorm (sigh): spurious factor of n before
    log determinant (should consider removing dmvnorm and just referring
    users to mvtnorm::dmvnorm)
    \item	added cross-ref to mvtnorm::dmvnorm in dmvnorm.Rd
}}
\section{Changes in version 1.0.8}{
  \itemize{
    \item fixed major bugs in dmvnorm
    \item added bbmle to get.emdbook.packages
    \item added an extended version of trcoef (transform coefficients)
    \item added (preliminary) perturb.params
    \item tweaked some keywords/aliases
    \item added traditional parameterization for dbetabinom
    \item added note about SuppDists to ?dbetabinom
    \item removed coda dependency
    \item put myxo1.bug into BUGS directory
}}
\section{Changes in version 1.0.7}{
  \itemize{
    \item tweaked dbetabinom documentation
    \item fixed (major) bug in get.emdbook.packages
    \item added aliases to data .Rd files
    \item    added reference to gsl package in lambertW Rd file
    \item    added gridsearch2d
    \item    added metropSB (Metropolis-Szymura-Barton)
    }}
	

    \section{Changes in version 1.0.6}{
      \itemize{
	\item internal tweaks
	\item tiny doc fix for Myxo (sum instead of bsum)
	}}

	\section{Changes in version 1.0.5}{
	  \itemize{
	    \item removed qcredint (HPDinterval[coda] is better)
	    \item added calcslice, HPDregionplot, *chibarsq functions
	    \item fixed deltavar
	    \item tweaks to scinot() 
	    \item added lump.mcmc.list()
	    \item changed dataset names, moved lily and goby data sets to emdbookx
	    \item added require(coda) to as.mcmc.bugs
	    \item added {d,r}betabinom
	    \item added SeedPred_mass data (numeric vector of seed masses)
	    \item added numeric derivatives to deltavar
	    \item fixed bug in calcslice (not referring to correct function)
	    \item added Lily data (moved from emdbookx)
	    }}

	    \section{Changes in version 1.0.4}{
	      \itemize{
		\item added r/dzinbinom (zero-inflated negative binomial) functions
		\item lseq
		\item scinot
		\item added myxomatosis data
		\item added "none" option to curve3d
		\item added credible interval functions, as.mcmc.bugs
}}

\section{Changes in version 1.0.3}{
  \itemize{
    \item added more data (fir and fish survivorship)
    \item delta-method
    functions
    \item get-packages function
    }}

    \section{Changes in version 1.0.1}{
      \itemize{
	\item added Vonesh data
	\item lambertW function (from bbmisc)
	}}
