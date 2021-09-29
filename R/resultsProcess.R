#' Processing peeling results
#'
#'
#' @param peel.results peeling results
#' @param posDT a data frame containing gene annotation information; a list
#'   component created by \code{dataPrep}.
#'
#' @return processed peeling results with a list of genes corresponding to each
#'   peeled region
#'
#' @seealso \code{\link{dataPrep}}
#'
#' @import plyr
#' @export
resultsProcess = function(peel.results,posDT) {
  posDT$start = as.numeric(posDT$start)
  posDT$end = as.numeric(posDT$end)
  inner.fun = function(x) {
    x = as.numeric(x)
    chr.rows = which(posDT$chr == x[1])
    gene.rows1 = chr.rows[((posDT$start[chr.rows]-x[2])*(posDT$end[chr.rows]-x[3]) <= 0)]
    gene.rows2 = chr.rows[((posDT$start[chr.rows]-x[2])*(posDT$start[chr.rows]-x[3]) <= 0)]
    gene.rows3 = chr.rows[((posDT$end[chr.rows]-x[2])*(posDT$end[chr.rows]-x[3]) <= 0)]
    gene.rows = unique(c(gene.rows1,gene.rows2,gene.rows2))
    gene.result = data.frame(t(sort(rownames(posDT)[gene.rows])))
    colnames(gene.result) = paste0("G",1:length(gene.rows))
    return(gene.result)
  }
  gene.list = plyr::alply(peel.results[,c("chr","peelStart","peelStop")],1,inner.fun)
  result = t(cbind(peel.results,do.call(plyr::rbind.fill,gene.list)))
  return(result)
}
