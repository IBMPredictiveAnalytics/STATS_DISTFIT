#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2011, 2014
#US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# Version: 1.1.2
# Author: JKP, IBM SPSS

# history
# 26-Nov-11 Add support for initial parameter estimates
# 23-Nov-12 Add logistic to distribution lists
# 23-Nov-12 Stop silliness of specifying both discrete and continuous distributions
# 08-Dec-14 localization changes, support split files

helptext="The STATS DISTFIT command requires the R Integration Plug-in
and the R MASS package.

STATS DISTFIT  VARIABLE=variable list 
DISTRIBUTION=distribution names
[/OPTIONS [SIMULATEP={NO*|YES}] [MINVALUE=value] [MAXVALUE=value]
[QQPLOT={NO*|YES}] PARM1=value PARM2=value].

VARIABLE is a list of numeric variables for which distributions should be estimated.

DISTRIBUTION is a list of one or more distributions to fit.  The distribution keywordss are
beta, cauchy, chisquared, exponential, f, gamma, geometric, lognormal, logistic, negativebinomial,
normal, poisson, t, and weibull.

For certain distributions starting values for estimating the parameters must be supplied using
the PARM1 and PARM2 keywords.

distribution    PARM1   PARM2
beta              shape1    shape2
chi-squared     df
f                   df1         df2

For some continuous distributions, a Kolmogorov-Smirnow test of fit is carried out, and the results are
displayed in a table for each variable.

For some discrete distributions, a chi-squared test is computed and displayed in a table for each
variable.  By default, asymptotic significance is reported.  Specify SIMULATEP=YES to use
Monte Carlo simulation.

For discrete distributions, counts for cells between the observed minimum and maximum 
values are included as zeros in the chi-squared test.  
Use  MINVALUE and/or MAXVALUE to override these values.

Use QQPLOT=YES to produce a Q-Q plot for each variable against each fitted distribution.

Split files and weight are not honored by this command.

STATS DISTFIT /HELP prints this information and does nothing else.

Example:
STATS DISTFIT VARIABLE=x y z DISTRIBUTION=lognormal negativebinomial.
"
# distribution functions for tests
cdfs=list("normal"="pnorm", "lognormal"="plnorm", "uniform"="punif", "Poisson"="ppois", "exponential"="pexp",
"weibull"="pweibull", "chi-squared"="pchisq", "cauchy"="pcauchy", "gamma"="pgamma", "negative binomial"="pnbinom",
"t"="pt", "beta"="pbeta", "f"="pf", "logistic"="plogis")
# quantile functions
qfs=list("normal"="qnorm", "lognormal"="qlnorm", "uniform"="qunif", "Poisson"="qpois", "exponential"="qexp",
"weibull"="qweibull", "chisq"="qchisq", "cauchy"="qcauchy", "gamma"="qgamma", "negative binomial"="qnbinom",
"t"="qt", "beta"="qbeta", "f"="qf", "logistic"="qlogis")

discretepdfs = list("Poisson"=dpois, "chisq"=dchisq, "negative binomial"=dnbinom, "geometric"=pgeom)

fitit<-function(variables, distribution, usermin=NULL, usermax=NULL, simulatep=FALSE, qplot=FALSE,
    shape1=NULL, shape2=NULL, df=NULL, df1=NULL, df2=NULL) {

    procname=gtxt("Distribution Fit")
    warningsprocname = gtxt("Distribution Fit: Warnings")
    omsid="STATSDISTFIT"
    warns = Warn(procname=warningsprocname,omsid=omsid)
    setuplocalization("STATS_DISTFIT")
    tryCatch(library(MASS), error=function(e) {
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "MASS"),
            dostop=TRUE)
        }
    )

    splitvars = spssdata.GetSplitVariableNames()
    StartProcedure(procname,"STATS FITDIST")

    while (!spssdata.IsLastSplit()) {
        dta<-spssdata.GetSplitDataFromSPSS(c(variables, splitvars),missingValueToNA=TRUE)
        numcases = nrow(dta)
        # if only one split variable, reference returns a simple scalar
        splitlabel = makesplitlabel(splitvars, as.data.frame(dta[1, splitvars]))

        # variables loop
        i = 1
        for (var in variables) {
            ksres = matrix(NA, length(distribution), 2)   # fit result structures for continuous and discrete
            chires = matrix(NA, length(distribution), 6)
            jcont = 0
            jdisc = 0
            
            dtai = dta[!is.na(dta[,i]), i]
            numcasesi = length(dtai)
            qpts = seq(0, 1, .01)    # for quantile plots
            
            # distributions loop
            continuousgiven = FALSE
            discretegiven = FALSE
            for (dis in distribution) {
                if (dis == "negativebinomial") {dis = "negative binomial"}
                if (dis == "chisquared") {dis = "chi-squared"}
                if (dis == "poisson") {dis = "Poisson"}
                cdf = cdfs[dis][[1]]
                if (is.null(cdf)) {
                    warns$warn(gtxtf("An unsupported distribution was specified: %s", dis),
                        dostop=TRUE)
                }
                discretedis = discretepdfs[[dis]]   # NULL if not a discrete distribution
                if (is.null(discretedis)) {
                    continuousgiven = TRUE
                }
                else {
                    discretegiven = TRUE
                }
                if (discretegiven && continuousgiven) {
                    warns$warn("Both continuous and discrete distributions were specified.", 
                        dostop=TRUE)
                }
                start = NULL            
                if (dis == "beta") {
                    if (is.null(shape1) || is.null(shape2)) {
                        warns$warn(gtxtf("shape1 and shape2 parameters must be supplied for %s", dis), 
                            dostop=TRUE)
                    }
                    start=list(shape1 = shape1, shape2 = shape2)
                }
                if (dis == "f") {
                    if (is.null(df1) || is.null(df2)) {
                        warns$warn(gtxtf("df1 and df2 parameters must be supplied for %s distribution", dis), 
                            dostop=TRUE)
                    }
                    start=list(df1 = df1, df2 = df2)
                }
                if (dis == "chi-squared") {
                    if (is.null(df)) {
                        warns$warn(gtxtf("df parameter must be supplied for %s distribution", dis), 
                            dostop=TRUE)
                    }
                    start = list(df = df)
                }
                discretedis = discretepdfs[[dis]]   # NULL if not a discrete distribution
                
                # fit the distribution after excluding missing values
                # start parameters are only used for a few distributions
                res = tryCatch(fitdistr(dtai, dis, start),
                    error=function(e) {
                    warns$warn(gtxtf("The distribution could not be fit.  Variable: %s, Distribution: %s", var, dis),
                        dostop=FALSE)
                    return(NULL)
                    }
                )
                    
                if (!is.null(res)) {  # got a result?
                    # display pivot table for distribution parameters for current variable
                    df = data.frame(rbind(coef(res), res$sd), row.names=c(gtxt("Estimate"), gtxt("Std. Error")))
                    spsspivottable.Display(df, 
                        title=gtxtf("Variable: %s Distribution: %s", var, dis),
                        templateName="DISTPARMS",
                        caption = gtxtf("Number of valid cases: %s, out of %s total cases", numcasesi, numcases),
                        outline=gtxt("Distribution Parameters"),
                        rowdim=gtxt("Parameters"),
                        coldim=gtxt("Statistics")
                    )
                    if (qplot) {
                        parmlist = list(qpts)
                        for (ii in 1:length(res$estimate)) parmlist[ii+1] = res$estimate[[ii]]
                        ow = options("warn")
                        options(warn=-1)
                        qtiles = do.call(qfs[dis][[1]], parmlist)
                        qqplot(qtiles, dtai, 
                            main=gtxtf("Q-Q Plot. Variable: %s  Distribution: %s", var, dis),
                            sub=splitlabel,
                            xlab = gtxt("Theoretical Quantiles"), ylab=gtxt("Sample Quantiles"),
                            cex=1.5, col="blue", pch=16)  # larger, solid, blue points
                        options(ow)
                    }
                    # accumulate fit statistics
                    if (is.null(discretedis)) {  #continuous
                        # Accumulate K-S test statistics for continuous distributions
                        parmlist = list(dtai, cdf)
                        for (ii in 1:length(res$estimate)) parmlist[ii+2] = res$estimate[[ii]]
                        
                        ksr = tryCatch(do.call(ks.test, parmlist),  error=function(e) {return(NULL)})
                        if (is.null(ksr)) {
                            ks = NA
                            ksp = NA
                        } else {
                            ks = ksr$statistic[[1]]
                            ksp = ksr$p.value
                        }
                        jcont = jcont+1
                        ksres[jcont,] = c(ks, ksp)
                    } else {
                        # chisq for discrete dis
                        tbl = table(dtai)   # accumulate distribution
                        tbllevels = as.numeric(names(tbl))  # get the levels as numeric
                        tblvalues = as.vector(tbl)  # get the counts
                        
                        # fill in any gaps in range as extended by user spec
                        tblmin = min(tbllevels, usermin)
                        tblmax = max(tbllevels, usermax)
                        othervalues = setdiff(tblmin:tblmax, tbllevels)
                        tbllevels = append(tbllevels, othervalues)
                        tblvalues = append(tblvalues, rep(0, length(othervalues)))
                        if (dis == "negative binomial") { # fitdistr returns the alternate parameterization for this
                            probs = dnbinom(tbllevels, size=res$estimate[1], mu=res$estimate[2])
                        } else {
                            parmlist = list(tbllevels)
                            for (ii in 1:length(res$estimate)) parmlist[ii+1] = res$estimate[[ii]]
                            probs = do.call(discretedis, parmlist)
                        }
                        ow = options("warn")
                        options(warn=-1) # suppress warnings from chisq.test - handled via table caption
                        tryCatch(
                            cres <- chisq.test(tblvalues, p=probs, rescale.p=TRUE, simulate.p.value=simulatep))
                        options(ow)
                        jdisc = jdisc+1
                        chires[jdisc,] = c(cres$statistic[[1]], cres$parameter[[1]], cres$p.value[[1]], sum(probs), tblmin, tblmax)
                        method = cres$method
                    }
                }
            }
                # display fit pivot table for all continuous distributions requested
                if (jcont > 0) {
                    ksdf = data.frame(ksres, row.names=distribution)
                    spsspivottable.Display(ksdf, title = gtxtf("Fit Statistics for Variable: %s - Continuous Distributions", var),
                        rowdim = gtxt("Distribution"),
                        collabels=c(gtxt("K-S Statistic"), gtxt("P Value")),
                        templateName="DISTSTATSCONT",
                        outline = gtxt("Continuous Distributions Fit Summary")
                    )
                }
                
                # display fit pivot table for all discrete distributions requested
                if (jdisc > 0) {
                    chires = data.frame(chires, row.names=distribution)  # need to filter
                    smallcounts = sum(tblvalues < 5)
                    spsspivottable.Display(chires, 
                        title = gtxtf("Fit Statistics for Variable: %s - Discrete Distributions", var),
                        rowdim = gtxt("Distribution"),
                        caption=gtxtf("Chi-squared calculated over all integer values between minimum and maximum.\nCell counts < 5: %s",
                            smallcounts),
                        collabels=c(gtxt("Chi-Squared Statistic"), gtxt("D.F."), gtxt("P Value"), 
                            gtxt("Sum of Probabilities"), gtxt("Minimum Value"), gtxt("Maximum Value")
                        ),
                        templateName="DISTSTATSDISC",
                        outline = gtxt("Discrete Distributions Fit Summary")
                    )
                }
            i = i+1  # variable index
        }
    }
    spssdata.CloseDataConnection()
    spsspkg.EndProcedure()
    # clean up workspace
    res <- tryCatch(rm(list=ls()),warning=function(e){return(NULL)})
}

# override for api to account for extra parameter in V19 and beyond
StartProcedure<-function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
       spsspkg.StartProcedure(procname,omsid)
       }
    else {
       spsspkg.StartProcedure(omsid)
       }
}

setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

gtxt <- function(...) {
    return(gettext(...,domain="STATS_DISTFIT"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_DISTFIT"))
}
Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment
    
    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.
        
        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 
        
        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any
        
        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spss.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                },
                error = function(e) {
                    FALSE
                }
                )
            } else {
                procok = TRUE
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                                               gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)
                
                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                                                spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

makesplitlabel = function(splitvars, splitvalues) {
    # return a possibly empty label for the split
    # splitvars is the list of split variables and may be NULL
    # splitvalues is a row of data for the split

    if (is.null(splitvars)) {
        return("")
    }
    vals = list()
    for (i in 1:length(splitvars)) {
        if (is.factor(splitvalues[1,i])) {
            vals[i] = levels(splitvalues[1,i])[splitvalues[1,i]]
        } else {
            vals[i] = splitvalues[1,i]
        }
    }
    return(paste(splitvars, vals, sep="=", collapse=" "))
}

Run<-function(args){
    
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
            spsspkg.Template("VARIABLE", subc="",  ktype="existingvarlist", var="variables", islist=TRUE),
            spsspkg.Template("DISTRIBUTION", subc="",  ktype="str", var="distribution", islist=TRUE),
            spsspkg.Template("SIMULATEP", subc="OPTIONS", ktype="bool", var="simulatep"),
            spsspkg.Template("MINVALUE", subc="OPTIONS", ktype="int", var="usermin"),
            spsspkg.Template("SHAPE1", subc="", ktype="float", var="shape1"),
            spsspkg.Template("SHAPE2", subc="", ktype="float", var="shape2"),
            spsspkg.Template("DF", subc="", ktype="float", var="df"),
            spsspkg.Template("DF1", subc="", ktype="float", var="df1"),
            spsspkg.Template("DF2", subc="", ktype="float", var="df2"),
            spsspkg.Template("QQPLOT", subc="OPTIONS", ktype="bool", var="qplot"),
            spsspkg.Template("MAXVALUE", subc="OPTIONS", ktype="int", var="usermax")
            ))
                
    if ("HELP" %in% attr(args,"names"))
        #writeLines(helptext)
        helper(cmdname)
    else
        res <- spsspkg.processcmd(oobj,args,"fitit")
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
