stabsel <- function(formula, data, family = "gaussian",
                    q = ceiling(sqrt(ncol(data))), B = 100,
                    thr = .9, ...) {

    ## Scale data
    for(i in seq(ncol(data))) {
        if(inherits(data[,i], "numeric"))
            data[,i] <- scale(data[,i])
    }

    ## Stability Selection
    ## TODO: parallel option
    stabselection <- NULL
    for(i in seq(B)) {
        cat(sprintf("Stability selection boosting run %d / %d \r", i, B))
        xx <- StabStep(formula = formula, data = data,
                       family  = family, q = q, seed = i, ...)
        stabselection <- c(stabselection, xx$sel)
    }
    cat("\n")
    formula <- xx$formula
    family  <- xx$family

    ## Re-build formula
    tabsel <- sort(table(stabselection), decreasing = FALSE)
    f <- StabFormula(tabsel, formula, family, thr, B)

    ## Comupte per-family-error-rate
    p <- 0
    for(i in seq(length(formula))) {
        p <- p + length(attr(terms(formula[[i]]$formula), "term.labels"))
    }
    PFER <- (q^2) / ((2*thr - 1)*p)

    ## Return
    rval <- list("table"       = tabsel, "raw" = stabselection,
                 "formula.org" = formula,
                 "formula.new" = f,
                 "family"      = family,
                 "parameter"   = list("q" = q, "B" = B, "thr" = thr,
                                      "p" = p, "PFER" = PFER))
    class(rval) <- c("stabsel", "list")
    return(rval)
}

StabStep <- function(formula, data, family = "gaussian", q, seed = NULL, ...) {
    if(!is.null(seed))
        set.seed(seed)
    d <- data[sample(nrow(data), size = round(nrow(data)/2)), ]
    b <- bamlss(formula, data = d, family = family, optimizer = boost,
                sampler = FALSE,
                binning = TRUE, maxq = q, scale.d = FALSE,
                plot = FALSE, verbose = FALSE, ...)

    rval <- NULL
    for (model in b$family$names) {
        boosum <- b$model.stats$optimizer$boost.summary$summary[[model]]
        selterms <- rownames(boosum)[boosum[,1] > 0]

        ## ignore intercept
        selterms <- selterms[!(selterms == "(Intercept)")]

        ## grep nets
        nterms   <- selterms[grepl("^n\\(", selterms)]
        selterms <- selterms[!grepl("^n\\(", selterms)]
        for(nterm in nterms) {
            cn <- colnames(b$parameters)[grepl(nterm, colnames(b$parameters), fixed = TRUE)]
            nr <- nrow(b$parameters)
            nsel <- cn[which(b$parameters[nr, cn] != 0)]
            selterms <- c(selterms, regmatches(nsel, regexpr("n\\([^\\)]*\\)\\..*$", nsel)))
        }

        ## grep names of parametric terms
        labels <- attr(b$terms[[model]], "term.labels")
        pID <- !grepl("^[a-z]*\\(", labels)
        pterms <- labels[pID]
        foo <- function(x) {
            if(class(b$model.frame[[x]]) == "factor") {
                rval <- rep(x, nlevels(b$model.frame[[x]]))
                names(rval) <- paste0(x, levels(b$model.frame[[x]]))
            } else {
                rval <- x
                names(rval) <- x
            }
            rval
        }
        facLabels <- NULL
        for (x in pterms) { facLabels <- c(facLabels, foo(x)) }

        ## merge factor dummies
        assign <- attr(b$x[[model]]$model.matrix, "assign")
        names(assign) <- colnames(b$x[[model]]$model.matrix)
        assign <- assign[names(assign) != "(Intercept)"]
        for(j in unique(assign)) {
            labels <- names(assign)[assign == j]
            if(any(selterms %in% labels)) {
                selterms <- selterms[!(selterms %in% labels)]
                res <- facLabels[names(facLabels) %in% labels][1]
                names(res) <- NULL
                selterms <- c(selterms, res)
            }
        }
        if (length(selterms) > 0)
            rval <- c(rval, paste(selterms, model, sep = "."))
    }

    rval <- list("sel"     = rval,
                 "family"  = family(b),
                 "formula" = formula(b))
    return(rval)
}

StabFormula <- function(tabsel, formula, family, thr, B) {
    ## Select terms
    p <- names(tabsel)[tabsel > (thr*B)]
    ## grep model identifier
    modelID <- substring(p, regexpr("\\.[a-z0-9]+$", p) + 1)
    ## skip model identifier
    p <- gsub("\\.[a-z0-9]+$", "", p)

    family <- bamlss.family(family)
    models <- names(formula)
    f <- list()
    for (i in seq_along(models)) {
        p2 <- p[modelID == models[i]]

        ## replace variables with term.labels
        labels <- attr(terms(formula)[[models[i]]], "term.labels")
        if(is.null(labels) & length(models) == 1)
            labels <- attr(terms(formula), "term.labels")

        ## smooth terms
        sID <- grepl("^[a-z]+[(]", labels)
        pure <- sapply(labels[sID], function(x) { eval(parse(text=x))$label })
        names(pure) <- NULL
        rhs <- labels[sID][pure %in% p2]
        ## parametric terms
        p2para <- p2[!(p2 %in% pure[pure %in% p2])]
        paral <- unlist(lapply(labels[!sID], function(pat) { any(grepl(pat, p2para)) }))
        rhs <- c(rhs, labels[!sID][paral])
        ## stack formula
        rhs <- paste(rhs, collapse = " + ", sep="")
        if (rhs == "")
            rhs <- "1"
        
        ## response
        resp <- all.vars(formula[[models[i]]]$formula)[1]
        lhs <- ifelse(i == 1, resp, "")
        if(family$family == "multinomial")
            lhs <- formula[[models[i]]]$response

        f[[models[i]]] <- eval(parse(text = paste(lhs, "~", rhs)))
        environment(f[[models[i]]]) <- NULL
    }

    return(f)
}

## methods for stabsel object
plot.stabsel <- function(x, show = NULL,
    pal = function(n) gray.colors(n, start = 0.9, end = 0.3), ...) {

    ## keep graphic parameters
    hold <- par(no.readonly = TRUE)
    on.exit(par(hold))

    tabsel <- x$table
    thr    <- x$parameter$thr
    B      <- x$parameter$B
    models <- x$family$names

    p <- names(tabsel)
    modelID <- substring(p, regexpr("\\.[a-z0-9]+$", p) + 1)        
    modelID <- factor(modelID, levels = models)
    names(tabsel) <- gsub("\\.[a-z0-9]+$", "", p)
    
    n <- length(tabsel)
    start <- ifelse(is.null(show), 1 , n - show + 1)

    col <- pal(nlevels(modelID))

    par(mar = c(5, 12, 4, 2) + .1)
    bp <- barplot(tabsel[start:n], horiz = TRUE, las = 1,
                  col = col[modelID[start:n]], ...)
    abline(v = thr*B, col = 1, lty = 3, lwd = 2)
    legend("bottomright", fill = col, legend = levels(modelID), bg = "white")
    title(xlab = "Frequency")

    invisible(bp)
}

print.stabsel <- function(x, ...) {
    cat("Stability selection:\n")
    cat("---\n")
    cat(sprintf("Family: %s\n", x$family$family))
    cat("---\n")
    cat(sprintf("Per family error rate: %.3f\n", x$parameter$PFER))
    cat("---\n")
    cat("Selected formula:\n")
    print(x$formula.new)
}

#summary.stabel <- function(x, ...) {}
#print.summary.stabel <- function(x, ...) {}

