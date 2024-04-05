## Copyright 2024 Emmanuel Paradis

## from pegas:
getIandJ <- function(ij, n)
{
    b <- n - 0.5
    i <- ceiling(b - sqrt(b * b - 2 * ij))
    j <- n * (1 - i) + (i + 1) * i/2 + ij
    as.integer(c(i, j))
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

allMST <- function (d, log = FALSE, quiet = FALSE)
{
    if (is.matrix(d)) d <- as.dist(d)
    uD <- sort(unique.default(d))
    n <- attr(d, "Size")
    forest <- seq_len(n)
    m <- matrix(NA_real_, n * (n - 1)/2, 3L)
    alt <- integer()
    i <- 0L

    nMST <- as.integer(!log)

    for (x in uD) {
        if (!quiet) cat(sprintf("\rdist = %f (max = %f)", x, uD[length(uD)]))
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
                ## maybe need to check the sign?
                if (!is.finite(nMST.cl)) next
                if (nMST.cl == 0) {
                    nMST <- nMST + sum(log(w[K]))
                } else { # there are loops
                    ## <FIXME> NEED TO ELABORATE THIS MORE?? (what if
                    ## several loops; see notes)
                    if (all(w[K] == 1)) {
                        nMST <- nMST + nMST.cl
                    } else {
                        ## nMST <- nMST + log(exp(nMST.cl) + sum(w[K] - 1))
                        ln_b <- log(sum(w[K] - 1)) # = ln(b) in the ln(a + b) when only ln(a) and ln(b) are known
                        ## so no need to compute exp(nMST.cl)
                        nMST <- nMST + nMST.cl + log1p(exp(ln_b - nMST.cl))
                    }
                    ## </FIXME>
                }
            } else {
                nMST.cl <- round(det(L[-1, -1, drop = FALSE]))
                if (nMST.cl == 0) next
                if (nMST.cl == 1) {
#                    cat(" == nMST.cl = ", nMST.cl, " ", dim(L)[1] -1, " ", length(K), "\n", sep = "")
                    nMST <- nMST * prod(w[K])
                } else { # there are loops
#                    cat("boucle...", nMST.cl, "\n")
                    ## <FIXME> NEED TO ELABORATE THIS MORE?? (what if
                    ## several loops; see notes)
                    nMST <- nMST * (nMST.cl + sum(w[K] - 1))
                    ## </FIXME>
#                    cat(" nMST =", nMST, "\n")
                }
#                if (!quiet) cat(" nMST =", nMST, "\n")
            }
        }

        ## back to the normal RMST code!
        for (j in seq_len(Npairs)) {
            h1 <- pairs[1L, j]
            h2 <- pairs[2L, j]
            k1 <- forest[h1]
            k2 <- forest[h2]
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

        diags <- length(unique(forest))
        if (!quiet) cat(sprintf(" Nb of edges = %d\n", i)) # cat(sprintf("Nb of clusters = %d\n", diags))
        if (diags == 1L) break
        ## if (i == n - 1) break
    }

    if (!quiet) cat("\n")

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

bigallMST <- function (X, log = FALSE, quiet = FALSE)
{
    n <- nrow(X)
    n_buff <- 10 * n
    forest <- seq_len(n)
    m <- matrix(NA_real_, n_buff, 3L)
    alt <- integer()
    im <- 0L

    nMST <- as.integer(!log)

    if (!is.loaded("Euclidean_dist_some_pairs_Call"))
        dyn.load("~/data/Projets/Sentinel2/HMM/Euclidean_dist_some_pairs.so")

    v <- seq_len(n)
    lmind <- numeric(n - 1)
    Lb <- vector("list", n - 1)
    lastmind <- -Inf

    repeat {
        for (i in 1:(n - 1)) {
            ##if (!quiet && !(i %% 1)) cat("\r  ", i, " / ", n, sep = "")
            j <- v[-seq_len(i)]
            ij <- cbind(i, j)
            d <- .Call("Euclidean_dist_some_pairs_Call", X, ij - 1L)
            dsel <- which(d > lastmind)
            if (!length(dsel)) {
                lmind[i] <- Inf
            } else {
                mind <- min(d[dsel])
                lmind[i] <- mind
                Lb[[i]] <- j[which(d == mind)]
            }
        }
        ##if (!quiet) cat("\n")
        lastmind <- min(lmind) # the global minimum

        ii <- which(lmind == lastmind)
        jj <- Lb[ii]
        lx <- lengths(jj)
        jj <- unlist(jj)
        if (any(lx > 1)) ii <- rep(ii, lx)
        pairs <- rbind(ii, jj)

        ## similar to above code:
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
                ## maybe need to check the sign?
                if (!is.finite(nMST.cl)) next
                if (nMST.cl == 0) {
                    nMST <- nMST + sum(log(w[K]))
                } else { # there are loops
                    if (all(w[K] == 1)) {
                        nMST <- nMST + nMST.cl
                    } else {
                        ## nMST <- nMST + log(exp(nMST.cl) + sum(w[K] - 1))
                        ln_b <- log(sum(w[K] - 1)) # = ln(b) in the ln(a + b) when only ln(a) and ln(b) are known
                        ## so no need to compute exp(nMST.cl)
                        nMST <- nMST + nMST.cl + log1p(exp(ln_b - nMST.cl))
                    }
                }
            } else {
                nMST.cl <- round(det(L[-1, -1, drop = FALSE]))
                if (nMST.cl == 0) next
                if (nMST.cl == 1) {
                    nMST <- nMST * prod(w[K])
                } else { # there are loops
                    nMST <- nMST * (nMST.cl + sum(w[K] - 1))
                }
            }
        }

        ## back to the normal RMST code!
        for (j in seq_len(Npairs)) {
            h1 <- pairs[1L, j]
            h2 <- pairs[2L, j]
            k1 <- forest[h1]
            k2 <- forest[h2]
            if (k1 != k2) {
                forest[forest == k2] <- k1
                im <- im + 1L
                m[im, ] <- c(h1, h2, lastmind)
            } else {
                im <- im + 1L
                m[im, ] <- c(h1, h2, lastmind)
                alt <- c(alt, im)
            }
        }

        diags <- length(unique(forest))
        if (!quiet) cat(sprintf(" Nb of edges = %d\tNb of clusters = %d\n", im, diags))
        if (diags == 1L) break
    }

    if (!quiet) cat("\n")

    m <- m[seq_len(im), , drop = FALSE]
    colnames(m) <- c("", "", "step")
    if (length(alt)) {
        ALT <- m[alt, , drop = FALSE]
        m <- m[-alt, , drop = FALSE]
        attr(m, "alter.links") <- ALT
    }
    attr(m, "labels") <- as.character(v)
    attr(m, "nMST") <- nMST
    class(m) <- "haploNet"
    m
}
