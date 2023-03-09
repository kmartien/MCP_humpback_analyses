# This file has modified versions of the gelato() and gelatoPermFunc() functions
# from strataG. The strataG versions wouldn't run on my computer. Changes from
# the strataG version are commented.

gelato.KKM <- function (g, unknown.strata, nrep = 1000, min.sample.size = 5, 
          num.cores = 1) 
{
  all.strata <- getStrata(g)
  unknown.strata <- unique(as.character(unknown.strata))
  if (!all(unknown.strata %in% all.strata)) {
    stop("Some 'unknown.strata' could not be found in 'g'")
  }
  knowns <- sort(setdiff(all.strata, unknown.strata))
  #KKM change - just because I did't have access to this function and am only usig one core
#  num.cores <- .getNumCores(num.cores)
  strata.freq <- table(all.strata)
  
  print(unknown.strata)
  
  result <- sapply(unknown.strata, function(unknown) {
    unknown.ids <- names(all.strata)[all.strata == unknown]
    unknown.result <- sapply(knowns, function(known) {
      can.run <- (strata.freq[known] - length(unknown.ids)) >= 
        min.sample.size
      if (!can.run) 
        return(NULL)
      known.ids <- names(all.strata)[all.strata == known]
      known.g <- g[known.ids, , , drop = TRUE]
      fst.dist <- do.call(rbind, lapply(1:nrep, gelatoPermFunc, 
                                        known.ids = known.ids, unknown.ids = unknown.ids, 
                                        g = g, known.g = known.g))
#KKM change - I do't know why there is a second do.call(rbind...) call. Didn't work with it, worked without it
#      fst.dist <- do.call(rbind, fst.dist)
      fst.dist <- fst.dist[apply(fst.dist, 1, function(x) all(!is.na(x))), 
      ]
      if (nrow(fst.dist) < 2) {
        NULL
      }
      else {
        null.mean <- mean(fst.dist[, "null"])
        null.sd <- stats::sd(fst.dist[, "null"])
        null.lik <- stats::dnorm(fst.dist[, "obs"], null.mean, 
                                 null.sd)
        log.Lik <- sum(log(null.lik), na.rm = T)
        obs.median <- stats::median(fst.dist[, "obs"], 
                                    na.rm = T)
        obs.mean <- mean(fst.dist[, "obs"], na.rm = T)
        list(fst.dist = fst.dist, log.Lik.smry = c(log.Lik = log.Lik, 
                                                   mean.nreps = log.Lik/length(log.Lik), median = log(stats::dnorm(obs.median, 
                                                                                                                   null.mean, null.sd)), mean = log(stats::dnorm(obs.mean, 
                                                                                                                                                                 null.mean, null.sd))), norm.coefs = c(mean = null.mean, 
                                                                                                                                                                                                       sd = null.sd))
      }
    }, simplify = FALSE)
    log.Lik <- sapply(unknown.result, function(x) {
      if (is.null(x)) 
        NA
      else x$log.Lik.smry["median"]
    })
    lik <- exp(log.Lik - max(log.Lik, na.rm = T))
    assign.prob <- lik/sum(lik, na.rm = T)
    names(assign.prob) <- knowns
    list(assign.prob = assign.prob, likelihoods = unknown.result)
  }, simplify = F)

  print(unknown.strata)
  
  assign.prob <- as.data.frame(t(sapply(result, function(x) x$assign.prob)))
  assign.prob$assignment <- apply(assign.prob, 1, function(x) {
    colnames(assign.prob)[which.max(x)]
  })
  result <- lapply(result, function(x) x$likelihoods)
  list(assign.prob = assign.prob, likelihoods = result)
}

gelatoPermFunc <- function(i, known.ids, unknown.ids, g, known.g) {
  # select samples to self assign
  ran.knowns <- sample(known.ids, length(unknown.ids))
  
  # extract gtypes of base known strata
  known.to.keep <- setdiff(known.ids, ran.knowns)
  
  # gtypes for observed Fst
  obs.g <- g[c(unknown.ids, known.to.keep), , , drop = TRUE]
  
  # gtypes for null Fst
  st <- getStrata(known.g)
  st[ran.knowns] <- "<gelato.unknown>"
  setStrata(known.g) <- st
  
# KKM change - here are the lines I changed (originals commented out)  
  c(obs = unname(overallTest(obs.g, nrep = 0)$result[2,1]), 
    null = unname(overallTest(known.g, nrep = 0)$result[2,1])
  )
  #  c(obs = unname(statFst(obs.g, nrep = 0)$result["estimate"]), 
  #    null = unname(statFst(known.g, nrep = 0)$result["estimate"])
  #  )
}