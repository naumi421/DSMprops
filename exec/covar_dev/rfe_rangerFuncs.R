#### Ranger functions for RFE in caret
## FROM:
## https://github.com/topepo/caret/issues/555#issuecomment-268623585
##

rangerFuncs <-  list(summary = defaultSummary,
                     fit = function(x, y, first, last, ...) {
                       loadNamespace("ranger")
                       dat <- if(is.data.frame(x))
                         x else as.data.frame(x)
                       dat$.outcome <- y
                       ranger::ranger(.outcome ~ ., data = dat,
                                      importance = if(first) "impurity" else "none",
                                      probability = is.factor(y),
                                      write.forest = TRUE,
                                      ...)
                     },
                     pred = function(object, x)  {
                       if(!is.data.frame(x)) x <- as.data.frame(x)
                       out <- predict(object, x)$predictions
                       if(object$treetype == "Probability estimation") {
                         out <- cbind(pred = colnames(out)[apply(out, 1, which.max)],
                                      out)
                       }
                       out
                     },
                     rank = function(object, x, y) {
                       if(length(object$variable.importance) == 0)
                         stop("No importance values available")
                       imps <- ranger:::importance(object)
                       vimp <- data.frame(Overall = as.vector(imps),
                                          var = names(imps))
                       rownames(vimp) <- names(imps)

                       vimp <- vimp[order(vimp$Overall, decreasing = TRUE),, drop = FALSE]
                       vimp
                     },
                     selectSize = pickSizeBest,
                     selectVar = pickVars)
