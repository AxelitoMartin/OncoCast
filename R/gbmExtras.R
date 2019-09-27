#' gbm.perf.noprint
#' @export
gbm.perf.noprint <- function (object, plot.it = F, oobag.curve = FALSE, overlay = TRUE,
                              method)
{
  if (missing(method)) {
    method <- guess_error_method(object)
  }
  best.iter <- best_iter(object, method = method)
  return(best.iter)
}
environment(gbm.perf.noprint) <- asNamespace('gbm')



#########################

#' relative.influence.noprint
#' @export
relative.influence.noprint <- function (object, n.trees, scale. = FALSE, sort. = FALSE)
{
  if (missing(n.trees)) {
    if (object$train.fraction < 1) {
      n.trees <- gbm.perf.noprint(object, method = "test", plot.it = FALSE)
    }
    else if (!is.null(object$cv.error)) {
      n.trees <- gbm.perf.noprint(object, method = "cv", plot.it = FALSE)
    }
    else {
      n.trees <- length(object$train.error)
    }
    cat(paste("n.trees not given. Using", n.trees, "trees.\n"))
    if (object$distribution == "multinomial") {
      n.trees <- n.trees * object$num.classes
    }
  }
  get.rel.inf <- function(obj) {
    lapply(split(obj[[6]], obj[[1]]), sum)
  }
  temp <- unlist(lapply(object$trees[1:n.trees], get.rel.inf))
  rel.inf.compact <- unlist(lapply(split(temp, names(temp)),
                                   sum))
  rel.inf.compact <- rel.inf.compact[names(rel.inf.compact) !=
                                       "-1"]
  rel.inf <- rep(0, length(object$var.names))
  i <- as.numeric(names(rel.inf.compact)) + 1
  rel.inf[i] <- rel.inf.compact
  names(rel.inf) <- object$var.names
  if (scale.) {
    rel.inf <- rel.inf/max(rel.inf)
  }
  if (sort.) {
    rel.inf <- rev(sort(rel.inf))
  }
  return(rel.inf = rel.inf)
}
environment(relative.influence.noprint) <- asNamespace('gbm')
