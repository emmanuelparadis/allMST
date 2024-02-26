## from pegas:
getIndex <- function(i, j, n)
{
    if (i < j) n * (i - 1) - i * (i - 1)/2 + j - i else n * (j - 1) - j * (j - 1)/2 + i - j
}

.getClusters <- function(x)
{
    nodes <- unique(as.vector(x))
    n <- length(nodes)
    forest <- setNames(1:n, nodes)
    for (i in 1:nrow(x)) {
        h1 <- as.character(x[i, 1L])
        h2 <- as.character(x[i, 2L])
        k1 <- forest[h1]
        k2 <- forest[h2]
        forest[forest == k2] <- k1
    }
    UF <- unique(forest)
    N <- length(UF)
    L <- vector("list", N)
    for (i in 1:N) L[[i]] <- nodes[forest == UF[i]]
    L
}

allMST <- function (d, log = FALSE)
{
    if (is.matrix(d)) d <- as.dist(d)
    uD <- sort(unique.default(d))
    n <- attr(d, "Size")
    forest <- seq_len(n)
    m <- matrix(NA_real_, n * (n - 1)/2, 3L)
    alt <- integer()
    i <- 0L

###    L <- list()
    nMST <- as.integer(!log)

    for (x in uD) {

###        MAT <- matrix(NA_real_, 0, 2)

        pairs <- sapply(which(d == x), getIandJ, n = n)
        todrop <- forest[pairs[1L, ]] == forest[pairs[2L, ]]
        if (all(todrop)) next
        if (any(todrop)) pairs <- pairs[, !todrop, drop = FALSE]
        Npairs <- ncol(pairs)

        MAT <- cbind(forest[pairs[1L, ]], forest[pairs[2L, ]])
        MAT <- t(apply(MAT, 1, sort))
        uMAT <- unique(MAT)
        tmp <- paste(MAT[, 1L], MAT[, 2L], sep = "-")
        w <- tabulate(match(tmp, unique(tmp)))
        stopifnot(length(w) == nrow(uMAT))
        Clusters <- .getClusters(uMAT)
        for (cl in Clusters) {
            L <- diag(length(cl))
            ## dimnames(L) <- list(cl, cl)
            diag(L) <- 0
            K <- integer()
            for (k in 1:nrow(uMAT)) {
                ii <- match(uMAT[k, 1L], cl)
                jj <- match(uMAT[k, 2L], cl)
                if (!is.na(ii) && !is.na(jj)) {
                    L[jj, ii] <- L[ii, jj] <- -1
                    L[ii, ii] <- L[ii, ii] + 1
                    L[jj, jj] <- L[jj, jj] + 1
                    K <- c(K, k)
                }
            }
            if (log) {
                det.cl <- determinant(L[-1, -1, drop = FALSE])
                nMST.cl <- det.cl$modulus
                ## maybe need to check the sign
                if (!is.finite(nMST.cl)) next
                if (nMST.cl == 0) {
                    nMST <- nMST + sum(log(w[K]))
                } else { # there are loops
                    ## <FIXME> NEED TO ELABORATE THIS A BIT MORE: what if
                    ## several loops (see notes)
                    nMST <- nMST + nMST.cl
                    ## </FIXME>
                }
            } else {
                nMST.cl <- det(L[-1, -1, drop = FALSE])
                if (nMST.cl == 0) next
                if (nMST.cl == 1) {
                    nMST <- nMST * prod(w[K])
                } else { # there are loops
                    ## <FIXME> NEED TO ELABORATE THIS A BIT MORE: what if
                    ## several loops (see notes)
                    nMST <- nMST * nMST.cl
                    ## </FIXME>
                }
            }
        }

        ## back to the normal RMST code!
        for (j in seq_len(Npairs)) {
            h1 <- pairs[1L, j]
            h2 <- pairs[2L, j]
            k1 <- forest[h1]
            k2 <- forest[h2]

###            MAT <- rbind(MAT, c(k1, k2))

            if (k1 != k2) {
                forest[forest == k2] <- k1
                i <- i + 1L
                m[i, ] <- c(h1, h2, x)
            }
            else {
                i <- i + 1L
                m[i, ] <- c(h1, h2, x)
                alt <- c(alt, i)
            }
        }

###        L <- c(L, list(MAT))
###        names(L)[length(L)] <- as.character(x)

        if (length(unique(forest)) == 1L) break
    }

###    L <<- L

    m <- m[seq_len(i), , drop = FALSE]
    colnames(m) <- c("", "", "step")
    if (length(alt)) {
        ALT <- m[alt, , drop = FALSE]
        m <- m[-alt, , drop = FALSE]
        attr(m, "alter.links") <- ALT
    }
    attr(m, "labels") <- attr(d, "Labels")
    attr(m, "nMST") <- nMST
    class(m) <- "haploNet"
    m
}
