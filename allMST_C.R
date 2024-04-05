## Copyright 2024 Emmanuel Paradis

## if there is no tie in the distances, the MST is *not* returned:
## use pegas::mst() instead
allMST2 <- function(d, n = NULL, log = FALSE)
{
    if (!is.loaded("all_mst_C")) dyn.load("all_mst.so")
    if (is.null(n)) n <- attr(d, "Size")
    if (is.null(n)) stop("'d' is not a \"dist\" object")
    z <- .Call("all_mst_C", d, n, log)
    links <- cbind(z[[1]], z[[2]])
    ll <- z[[3]]
    alts <- z[[5]] == 1
    mst <- links[which(!alts), ]
    alt <- links[which(alts), ]
### <FIXME> the 'alt' vector is not always well set
### This does not work for the woodmouse data (see C code).
### </FIXME>
    nt <- cbind(mst, ll[which(!alts)])
    attr(nt, "alter.links") <- cbind(alt, ll[which(alts)])
    labs <- attr(d, "Labels")
    if (is.null(labs)) labs <- as.character(1:n)
    attr(nt, "labels") <- labs
    attr(nt, "N_MSTs") <- z[[4]] # nMST
    class(nt) <- "haploNet"
    nt
}
