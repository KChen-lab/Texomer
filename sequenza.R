read.seqz<-function (file, nrows = -1, fast = FALSE, gz = TRUE, header = TRUE, 
    colClasses = c("character", "integer", "character", "integer", 
        "integer", "numeric", "numeric", "numeric", "character", 
        "numeric", "numeric", "character", "character", "character"), 
    chr.name = NULL, n.lines = NULL, ...) 
{
    if (!is.null(n.lines) & is.null(chr.name)) 
        fast <- FALSE
    if (fast && nrows == -1) {
        if (!is.null(chr.name)) {
            wc <- system(paste(paste("grep -c \"^", chr.name, 
                "\t\"", sep = ""), file, sep = " "), intern = TRUE)
        }
        else {
             wc <- system(paste("wc", file), intern = TRUE)
        }
        if (is.null(chr.name)) {
            wc <- sub("^ +", "", wc)
            wc <- strsplit(wc, " ")[[1]][1]
        }
        nrows <- max(as.integer(wc), 1)
        message("Reading ", nrows, " lines...")
    }
    if (!is.null(chr.name)) {
        grep.part <- paste("grep '^", chr.name, "\t'", sep = "")
        seqz.data <- read.delim(pipe(grep.part), nrows = nrows, 
            colClasses = colClasses, header = FALSE, ...)
        if (header == TRUE) {
            head <- colnames(read.table(file, header = TRUE, 
                nrows = 1))
            colnames(seqz.data) <- head
        }
        seqz.data
    }
    else {
        if (!is.null(n.lines)) {
            if (!is.numeric(n.lines) | length(n.lines) != 2) 
                stop("n.lines must be a vector of 2 integers")
            n.lines <- round(sort(n.lines), 0)
            if (header == TRUE) {
                n.lines <- n.lines + 1
            }
            if (gz) {
                seqz.data <- read.delim(pipe(paste("gzip -d -c", 
                  file, "| sed -n '", paste(n.lines[1], n.lines[2], 
                    sep = ","), "p'")), colClasses = colClasses, 
                  nrows = 1 + n.lines[2] - n.lines[1], header = FALSE, 
                  ...)
            }
            else {
                seqz.data <- read.delim(pipe(paste("sed -n '", 
                  paste(n.lines[1], n.lines[2], sep = ","), "p'", 
                  file)), colClasses = colClasses, nrows = 1 + 
                  n.lines[2] - n.lines[1], header = FALSE, ...)
            }
            if (header == TRUE) {
                head <- colnames(read.table(file, header = TRUE, 
                  nrows = 1))
                colnames(seqz.data) <- head
            }
            seqz.data
        }
        else {
            read.delim(file, nrows = nrows, colClasses = colClasses, 
                header = header, ...)
        }
    }
}

gc.norm<-function (x, gc) 
{
    dr.by.gc <- split(x = x, f = gc)
    raw <- t(sapply(dr.by.gc, quantile, probs = c(0.25, 0.5, 
        0.75), na.rm = TRUE))
    dr.by.gc.mean <- sapply(dr.by.gc, FUN = mean, na.rm = TRUE)
    dr.by.gc.median <- sapply(dr.by.gc, FUN = median, na.rm = TRUE)
    adj <- sweep(raw, 1, dr.by.gc.median, "/")
    list(raw = raw, adj = adj, gc.values = as.numeric(names(dr.by.gc)), 
        raw.mean = dr.by.gc.mean, raw.median = dr.by.gc.median)
}

gc.sample.stats<-function (seqz.data) 
{
    gc.stats <- gc.norm(x = seqz.data$depth.ratio, gc = seqz.data$GC.percent)
    chr.ord <- unique(seqz.data$chromosome)
    chr.dim <- lapply(X = split(seqz.data$chromosome, seqz.data$chromosome), 
            FUN = length)
    chr.dim <- data.frame(chr = chr.ord, n.lines = do.call(rbind, 
                    chr.dim[chr.ord]))
    chr.dim$start <- cumsum(c(1, chr.dim$n.lines[-length(chr.dim$n.lines)]))
    chr.dim$end <- chr.dim$start + chr.dim$n.lines - 1
    gc.stats$file.metrics <- chr.dim
    gc.stats
}
windowValues<-function (x, positions, chromosomes, window = 1e+06, overlap = 0, 
    weight = rep.int(x = 1, times = length(x)), start.coord = 1) 
{
    weight <- sqrt(weight)
    overlap <- as.integer(overlap)
    window.offset <- as.integer(window - round(window * (overlap/(overlap + 
        1))))
    chr.ordered <- unique(chromosomes)
    data.splitByChr <- split(data.frame(pos = positions, x = x, 
        weight = weight), f = factor(chromosomes, levels = chr.ordered))
    lapply(data.splitByChr, function(data.oneChr) {
        range.pos <- range(data.oneChr$pos, na.rm = TRUE)
        if (!is.null(start.coord)) {
            range.pos[1] <- as.integer(start.coord)
        }
        beam.coords <- seq.int(range.pos[1], range.pos[2], by = window.offset)
        if (max(beam.coords) != range.pos[2]) {
            beam.coords <- c(beam.coords, range.pos[2])
        }
        nWindows <- length(beam.coords) - overlap - 1
        pos.cut <- cut(data.oneChr$pos, breaks = beam.coords)
        x.split <- split(data.oneChr$x, f = pos.cut)
        weight.split <- split(data.oneChr$weight, f = pos.cut)
        window.starts <- beam.coords[1:nWindows]
        window.ends <- beam.coords[(1:nWindows) + 1 + overlap]
        idx.list <- lapply(1:nWindows, function(ii) ii + (0:overlap))
        x.window <- lapply(idx.list, function(idx) unlist(x.split[idx], 
            use.names = FALSE))
        weight.window <- lapply(idx.list, function(idx) unlist(weight.split[idx], 
            use.names = FALSE))
        window.means <- mapply(weighted.mean, x = x.window, w = weight.window)
        window.quantiles <- sapply(x.window, quantile, probs = c(0.25, 
            0.75), na.rm = TRUE, names = FALSE)
        window.counts <- sapply(x.window, length)
        data.frame(start = window.starts, end = window.ends, 
            mean = window.means, q0 = window.quantiles[1, ], 
            q1 = window.quantiles[2, ], N = window.counts)
    })
}

find.breaks<-function (seqz.baf, gamma = 80, kmin = 10, baf.thres = c(0, 0.5), 
    verbose = FALSE, seg.algo = "aspcf", ...) 
{
    chromosome <- gsub(x = seqz.baf$chromosome, pattern = "chr", 
        replacement = "")
    logR = data.frame(chrom = chromosome, pos = seqz.baf$position, 
        s1 = log2(seqz.baf$adjusted.ratio))
    logR.wins <- copynumber::winsorize(logR, verbose = verbose)
    if (seg.algo == "aspcf") {
        BAF = data.frame(chrom = chromosome, pos = seqz.baf$position, 
            s1 = seqz.baf$Bf)
        allele.seg <- copynumber::aspcf(logR = logR.wins, BAF = BAF, 
            baf.thres = baf.thres, verbose = verbose, gamma = gamma, 
            kmin = kmin, ...)
    }
    else if (seg.algo == "pcf") {
        allele.seg <- copynumber::pcf(data = logR.wins, verbose = verbose, 
            gamma = gamma, kmin = kmin, ...)
    }
    else {
        stop("Supported segmentation algorithms are only 'aspcf' or 'pcf' from the copynumber package.")
    }
    if (length(grep("chr", seqz.baf$chromosome)) > 0) {
        allele.seg$chrom <- paste("chr", allele.seg$chrom, sep = "")
    }
    breaks <- allele.seg[, c("chrom", "start.pos", "end.pos", 
        "arm")]
    not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1], 
        0))
    breaks$end.pos[not.uniq] <- breaks$end.pos[not.uniq] - 1
    breaks
}
segment.breaks<-function (seqz.tab, breaks, min.reads.baf = 1, weighted.mean = TRUE) 
{
    if (weighted.mean) {
        w.r <- sqrt(seqz.tab$depth.normal)
        rw <- seqz.tab$adjusted.ratio * w.r
        w.b <- sqrt(seqz.tab$good.reads)
        bw <- seqz.tab$Bf * w.b
        seqz.tab <- cbind(seqz.tab[, c("chromosome", "position", 
            "zygosity.normal", "good.reads")], rw = rw, w.r = w.r, 
            bw = bw, w.b = w.b)
    }
    chr.order <- unique(seqz.tab$chromosome)
    seqz.tab <- split(seqz.tab, f = seqz.tab$chromosome)
    segments <- list()
    for (i in 1:length(seqz.tab)) {
        seqz.b.i <- seqz.tab[[i]][seqz.tab[[i]]$zygosity.normal == 
            "het", ]
        seqz.b.i <- seqz.b.i[seqz.b.i$good.reads >= min.reads.baf, 
            ]
        breaks.i <- breaks[breaks$chrom == names(seqz.tab)[i], 
            ]
        nb <- nrow(breaks.i)
        breaks.vect <- do.call(cbind, split.data.frame(breaks.i[, 
            c("start.pos", "end.pos")], f = 1:nb))
        unique.breaks <- function(b, offset = 1) {
            while (any(diff(b) == 0)) {
                b[which(diff(b) == 0) + 1] <- b[diff(b) == 0] + 
                  offset
            }
            b
        }
        breaks.vect <- unique.breaks(b = as.numeric(breaks.vect), 
            offset = 1)
        breaks.vect=unique(breaks.vect)
        fact.r.i <- cut(seqz.tab[[i]]$position, breaks.vect)
        fact.b.i <- cut(seqz.b.i$position, breaks.vect)
        seg.i.s.r <- sapply(X = split(seqz.tab[[i]]$chromosome, 
            f = fact.r.i), FUN = length)
        seg.i.s.b <- sapply(X = split(seqz.b.i$chromosome, f = fact.b.i), 
            FUN = length)
        if (weighted.mean) {
            seg.i.rw <- sapply(X = split(seqz.tab[[i]]$rw, f = fact.r.i), 
                FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.w.r <- sapply(X = split(seqz.tab[[i]]$w.r, 
                f = fact.r.i), FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.r.sd <- sapply(X = split(seqz.tab[[i]]$rw/seqz.tab[[i]]$w.r, 
                f = fact.r.i), FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.b.sd <- sapply(X = split(seqz.b.i$bw/seqz.b.i$w.b, 
                f = fact.b.i), FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.bw <- sapply(X = split(seqz.b.i$bw, f = fact.b.i), 
                FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.w.b <- sapply(X = split(seqz.b.i$w.b, f = fact.b.i), 
                FUN = function(a) sum(a, na.rm = TRUE))
            segments.i <- data.frame(chromosome = names(seqz.tab)[i], 
                start.pos = as.numeric(breaks.vect[-length(breaks.vect)]), 
                end.pos = as.numeric(breaks.vect[-1]), Bf = seg.i.bw/seg.i.w.b, 
                N.BAF = seg.i.s.b, sd.BAF = seg.i.b.sd, depth.ratio = seg.i.rw/seg.i.w.r, 
                N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd, stringsAsFactors = FALSE)
        }
        else {
            seg.i.r <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio, 
                f = fact.r.i), FUN = function(a) mean(a, na.rm = TRUE))
            seg.i.b <- sapply(X = split(seqz.b.i$Bf, f = fact.b.i), 
                FUN = function(a) mean(a, na.rm = TRUE))
            seg.i.r.sd <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio, 
                f = fact.r.i), FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.b.sd <- sapply(X = split(seqz.b.i$Bf, f = fact.b.i), 
                FUN = function(a) sd(a, na.rm = TRUE))
            segments.i <- data.frame(chromosome = names(seqz.tab)[i], 
                start.pos = as.numeric(breaks.vect[-length(breaks.vect)]), 
                end.pos = as.numeric(breaks.vect[-1]), Bf = seg.i.b, 
                N.BAF = seg.i.s.b, sd.BAF = seg.i.b.sd, depth.ratio = seg.i.r, 
                N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd, stringsAsFactors = FALSE)
        }
        segments[[i]] <- segments.i[seq(from = 1, to = nrow(segments.i), 
            by = 2), ]
    }
    segments <- do.call(rbind, segments[as.factor(chr.order)])
    row.names(segments) <- 1:nrow(segments)
    len.seg <- (segments$end.pos - segments$start.pos)/1e+06
    segments[(segments$N.ratio/len.seg) >= 2, ]
}
mutation.table<-function (seqz.tab, mufreq.treshold = 0.15, min.reads = 40, min.reads.normal = 10, 
    max.mut.types = 3, min.type.freq = 0.9, min.fw.freq = 0, 
    segments = NULL) 
{
    chroms <- unique(seqz.tab$chromosome)
    hom.filt <- seqz.tab$zygosity.normal == "hom"
    seqz.tab <- seqz.tab[hom.filt, ]
    reads.filt <- seqz.tab$good.reads >= min.reads & seqz.tab$depth.normal >= 
        min.reads.normal
    seqz.tab <- seqz.tab[reads.filt, ]
    mufreq.filt <- seqz.tab$Af <= (1 - mufreq.treshold)
    seqz.tab <- seqz.tab[mufreq.filt, ]
    if (!is.null(segments)) {
        for (i in 1:nrow(segments)) {
            pos.filt <- seqz.tab$chromosome == segments$chromosome[i] & 
                seqz.tab$position >= segments$start.pos[i] & 
                seqz.tab$position <= segments$end.pos[i]
            seqz.tab$adjusted.ratio[pos.filt] <- segments$depth.ratio[i]
        }
    }
    seqz.dummy <- data.frame(chromosome = chroms, position = 1, 
        GC.percent = NA, good.reads = NA, adjusted.ratio = NA, 
        F = 0, mutation = "NA", stringsAsFactors = FALSE)
    if (nrow(seqz.tab) >= 1) {
        mu.fracts <- mut.fractions(AB.tumor = seqz.tab$AB.tumor, 
            Af = seqz.tab$Af, tumor.strand = seqz.tab$tumor.strand)
        mufreq.filt <- mu.fracts$freq >= mufreq.treshold
        type.filt <- mu.fracts$base.count <= max.mut.types
        prop.filt <- mu.fracts$maj.base.freq >= min.type.freq
        if (!is.na(min.fw.freq)) {
            fw.2 = 1 - min.fw.freq
            fw.2 <- sort(c(fw.2, min.fw.freq))
            fw.filt <- mu.fracts$fw.freq > fw.2[1] & mu.fracts$fw.freq < 
                fw.2[2]
            mufreq.filt <- mufreq.filt & type.filt & prop.filt & 
                fw.filt
        }
        else {
            mufreq.filt <- mufreq.filt & type.filt & prop.filt
        }
        mut.type <- paste(seqz.tab$AB.normal, mu.fracts$base, 
            sep = ">")
        seqz.tab <- seqz.tab[, c("chromosome", "position", "GC.percent", 
            "good.reads", "adjusted.ratio")]
        seqz.tab <- cbind(seqz.tab, F = mu.fracts$freq, mutation = mut.type)
        rbind(seqz.tab[mufreq.filt, ], seqz.dummy)
    }
    else {
        seqz.dummy
    }
}

sequenza.extract<-function (seqz.data,gc.stats, chroso,window = 1e+06, overlap = 1, gamma = 80, 
        kmin = 10, gamma.pcf = 140, kmin.pcf = 40, mufreq.treshold = 0.1, 
        min.reads = 40, min.reads.normal = 10, min.reads.baf = 1, 
        max.mut.types = 1, min.type.freq = 0.9, min.fw.freq = 0, 
        verbose = TRUE, chromosome.list = NULL, breaks = NULL, breaks.method = "het", 
        assembly = "hg19", weighted.mean = TRUE) 
{
    chr.vect <- chroso
    windows.baf <- list()
    windows.ratio <- list()
    mutation.list <- list()
    segments.list <- list()
    coverage.list <- list()
    if (is.null(dim(breaks))) {
        breaks.all <- NULL
    }else {
        breaks.all <- breaks
    }
    if (is.null(chromosome.list)) {
        chromosome.list <- chr.vect
    }else {
        chromosome.list <- chromosome.list[chromosome.list %in% 
                        chr.vect]
    }
    for (chr in chromosome.list) {
        if (verbose) {
            message("Processing ", chr, ": ", appendLF = FALSE)
        }
        file.lines <- gc.stats$file.metrics[which(chr.vect == 
                                chr), ]
        seqz.hom <- seqz.data$zygosity.normal[seqz.data$chr==chr] == "hom"
        seqz.het <- seqz.data[seqz.data$chr==chr, ][!seqz.hom,]
        het.filt <- seqz.het$good.reads >= min.reads.baf
        seqz.r.win <- windowValues(x = seqz.data$adjusted.ratio[seqz.data$chr==chr], 
                positions = seqz.data$position[seqz.data$chr==chr], chromosomes = seqz.data$chromosome[seqz.data$chr==chr], 
                window = window, overlap = overlap, weight = seqz.data$depth.normal[seqz.data$chr==chr])
        if (nrow(seqz.het) > 0) {
            breaks.method.i <- breaks.method
            seqz.b.win <- windowValues(x = seqz.het$Bf[het.filt], 
                    positions = seqz.het$position[het.filt], chromosomes = seqz.het$chromosome[het.filt], 
                    window = window, overlap = overlap, weight = seqz.het$good.reads[het.filt])
            if (is.null(breaks.all)) {
                if (breaks.method.i == "full") {
                    breaks <- find.breaks(seqz.data[seqz.data$chr==chr,], gamma = gamma.pcf, 
                            assembly = assembly, kmin = kmin.pcf, seg.algo = "pcf")
                    breaks.het <- try(find.breaks(seqz.het, gamma = gamma, 
                                    assembly = assembly, kmin = kmin, baf.thres = c(0, 
                                            0.5)), silent = FALSE)
                    if (!is.null(breaks.het)) {
                        merge.breaks <- function(breaks, breaks.het) {
                            merged.breaks <- unique(sort(c(breaks$start.pos, 
                                                    breaks$end.pos, breaks.het$start.pos, 
                                                    breaks.het$end.pos)))
                            merged.breaks <- merged.breaks[diff(merged.breaks) > 
                                            1]
                            merged.start <- merged.breaks
                            merged.start[-1] <- merged.start[-1] + 
                                    1
                            breaks <- data.frame(chrom = unique(breaks$chrom), 
                                    start.pos = merged.start[-(length(merged.start))], 
                                    end.pos = merged.breaks[-1])
                        }
                        chr.p <- merge.breaks(breaks[breaks$arm == 
                                                "p", ], breaks.het[breaks.het$arm == "p", 
                                ])
                        chr.q <- merge.breaks(breaks[breaks$arm == 
                                                "q", ], breaks.het[breaks.het$arm == "q", 
                                ])
                        breaks <- rbind(chr.p, chr.q)
                    }
                }else if (breaks.method.i == "het") {
                    breaks <- try(find.breaks(seqz.het, gamma = gamma, 
                                    assembly = assembly, kmin = kmin, baf.thres = c(0, 
                                            0.5)), silent = FALSE)
                }else if (breaks.method.i == "fast") {
                    BAF <- data.frame(chrom = chr, pos = c(seqz.b.win[[1]]$start, 
                                    tail(seqz.b.win[[1]]$end, n = 1)), s1 = c(seqz.b.win[[1]]$mean, 
                                    tail(seqz.b.win[[1]]$mean, n = 1)))
                    BAF$s1[is.na(BAF$s1)] <- 0
                    logR <- data.frame(chrom = chr, pos = c(seqz.r.win[[1]]$start, 
                                    tail(seqz.r.win[[1]]$end, n = 1)), s1 = c(log2(seqz.r.win[[1]]$mean), 
                                    log2(tail(seqz.r.win[[1]]$mean, n = 1))))
                    not.cover <- is.na(logR$s1)
                    BAF <- BAF[!not.cover, ]
                    logR <- logR[!not.cover, ]
                    logR.wins <- copynumber::winsorize(logR, verbose = FALSE)
                    allele.seg <- copynumber::aspcf(logR = logR.wins, 
                            BAF = BAF, baf.thres = c(0, 0.5), verbose = FALSE, 
                            gamma = gamma, kmin = kmin)
                    if (length(grep("chr", chr)) > 0) {
                        allele.seg$chrom <- paste("chr", allele.seg$chrom, 
                                sep = "")
                    }
                    breaks <- allele.seg[, c("chrom", "start.pos", 
                                    "end.pos")]
                    not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1], 
                                    0))
                    breaks$end.pos[not.uniq] <- breaks$end.pos[not.uniq] - 
                            1
                }else {
                    stop("The implemented segmentation methods are 'full', 'het' and 'fast'.")
                }
            }
            else {
                breaks <- breaks.all[breaks.all$chrom == chr, 
                ]
            }
            if (!is.null(breaks) & nrow(breaks) > 0) {
                seg.s1 <- segment.breaks(seqz.tab = seqz.data[seqz.data$chr==chr,], 
                        breaks = breaks, min.reads.baf = min.reads.baf, 
                        weighted.mean = weighted.mean)
            }
            else {
                seg.s1 <- segment.breaks(seqz.data[seqz.data$chr==chr,], breaks = data.frame(chrom = chr, 
                                start.pos = min(seqz.data$position[seqz.data$chr==chr], na.rm = TRUE), 
                                end.pos = max(seqz.data$position[seqz.data$chr==chr], na.rm = TRUE)), 
                        weighted.mean = weighted.mean)
            }
        }
        else {
            seqz.b.win <- list()
            seqz.b.win[[1]] <- data.frame(start = min(seqz.data$position[seqz.data$chr==chr], 
                            na.rm = TRUE), end = max(seqz.data$position[seqz.data$chr==chr], 
                            na.rm = TRUE), mean = 0.5, q0 = 0.5, q1 = 0.5, 
                    N = 1)
            if (breaks.method == "full") {
                breaks <- find.breaks(seqz.data[seqz.data$chr==chr,], gamma = gamma.pcf, 
                        assembly = assembly, kmin = kmin.pcf, seg.algo = "pcf")
            }
            else {
                breaks = data.frame(chrom = chr, start.pos = min(seqz.data$position[seqz.data$chr==chr], 
                                na.rm = TRUE), end.pos = max(seqz.data$position[seqz.data$chr==chr], 
                                na.rm = TRUE))
            }
            seg.s1 <- segment.breaks(seqz.data[seqz.data$chr==chr,], breaks = breaks, 
                    weighted.mean = weighted.mean)
        }
        mut.tab <- mutation.table(seqz.data[seqz.data$chr==chr,], mufreq.treshold = mufreq.treshold, 
                min.reads = min.reads, min.reads.normal = min.reads.normal, 
                max.mut.types = max.mut.types, min.type.freq = min.type.freq, 
                min.fw.freq = min.fw.freq, segments = seg.s1)
        windows.baf[[which(chromosome.list == chr)]] <- seqz.b.win[[1]]
        windows.ratio[[which(chromosome.list == chr)]] <- seqz.r.win[[1]]
        mutation.list[[which(chromosome.list == chr)]] <- mut.tab
        segments.list[[which(chromosome.list == chr)]] <- seg.s1
        coverage.list[[which(chromosome.list == chr)]] <- data.frame(sum = sum(as.numeric(seqz.data$depth.tumor[seqz.data$chr==chr]), 
                        na.rm = TRUE), N = length(seqz.data$depth.tumor[seqz.data$chr==chr]))
        if (verbose) {
            message(nrow(mut.tab), " variant calls; ", nrow(seqz.het), 
                    " heterozygous positions; ", sum(seqz.hom), " homozygous positions.")
        }
    }
    names(windows.baf) <- chromosome.list
    names(windows.ratio) <- chromosome.list
    names(mutation.list) <- chromosome.list
    names(segments.list) <- chromosome.list
    coverage.list <- do.call(rbind, coverage.list)
    coverage <- sum(coverage.list$sum)/sum(coverage.list$N)
    return(list(BAF = windows.baf, ratio = windows.ratio, mutations = mutation.list, 
                    segments = segments.list, chromosomes = chromosome.list, 
                    gc = gc.stats, avg.depth = round(coverage, 0)))
}
mclapply<-function (X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, 
    mc.silent = FALSE, mc.cores = 1L, mc.cleanup = TRUE, mc.allow.recursive = TRUE) 
{
    cores <- as.integer(mc.cores)
    if (cores < 1L) 
        stop("'mc.cores' must be >= 1")
    if (cores > 1L) 
        stop("'mc.cores' > 1 is not supported on Windows")
    lapply(X, FUN, ...)
}

mclapplyPb <- function(X, FUN, mc.cores = getOption("mc.cores", 2L), ...) {
   env <- environment()
   pb_Total <- length(X) / mc.cores
   counter <- 0
   pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
   wrapper <- function(...){
      curVal <- get("counter", envir = env)
      assign("counter", curVal + 1 ,envir = env)
      setTxtProgressBar(get("pb", envir = env), curVal + 1)
      FUN(...)
   }
   res <- mclapply(X, wrapper, mc.cores = mc.cores, ...)
   close(pb)
   res
}
model.points <- function(cellularity, ploidy,
                         types, avg.depth.ratio) {
   mufreqs     <-  theoretical.mufreq(cellularity = cellularity , CNn = types[, 1], CNt = types[, 2], Mt = types[, 3])
   depth.ratio <-  theoretical.depth.ratio(cellularity = cellularity, ploidy = ploidy,
                                       CNt = types[, 2], CNn = types[, 1],
                                       avg.depth.ratio = avg.depth.ratio)
   cbind(mufreqs, depth.ratio)
}
theoretical.mufreq <- function(Mt, CNt, CNn = 2, cellularity) {
   normal.alleles <- (CNt - Mt) * cellularity + CNn * (1 - cellularity)
   all.alleles    <- (CNt * cellularity) + CNn * (1 - cellularity)
   1 - (normal.alleles / all.alleles)
}
theoretical.depth.ratio <- function(CNt, CNn = 2, cellularity, ploidy, normal.ploidy = 2, avg.depth.ratio = 1) {
   cellu.copy.term   <- (1 - cellularity) + (CNt / CNn * cellularity)
   ploidy.cellu.term <- (ploidy / normal.ploidy * cellularity) + 1 - cellularity
   avg.depth.ratio * cellu.copy.term / ploidy.cellu.term
}
expected.baf <- function(sd, ...) {
   baf      <- theoretical.baf(...)
   baf.t2 <- function(BAF, sd, by = 0.001){
      bafs   <- seq(0,1,0.001)
      b.b    <- dt2(x = bafs, mean = BAF, sd = sd, df = 5)
      b.a    <- dt2(x = bafs, mean = 1-BAF, sd = sd, df = 5)
      half.b <- bafs[bafs <= 0.5]
      b <- (b.b+b.a)[bafs <= 0.5]
      weighted.mean(half.b,b)
   }
   BAF <- mapply(FUN = baf.t2, baf$BAF,
                 MoreArgs = list(sd = sd))
   wgh <- dt2(x = baf$BAF, mean = 0.5, sd = sd, df = 5)
   wgh <- wgh/max(wgh)
   mean.bf <- function(x) {
      weighted.mean(x=c(x["BAF"], x["eBAF"]), w = c((1-x["wgh"]), x["wgh"]))
   }
   baf$BAF <- apply(cbind(BAF = baf$BAF, eBAF = BAF, wgh = wgh), 1, FUN = mean.bf)
   baf
}
theoretical.baf <- function(CNt, CNn = 2, cellularity) {
   alleles       <- seq(from = 1, to = CNt, by = 1)
   max.b <- function(CNt) {
      max.b.alleles <- CNt / 2
      if (CNt %% 2 != 0 ) {
         max.b.alleles <- trunc(max.b.alleles)
      }
      max.b.alleles
   }
   fract.normal.alleles <- (1 - cellularity)
   res                  <- list()
   for (i in 1:length(alleles)) {
      max.b.alleles <- max.b(alleles[i])
      max.a.alleles <- alleles[i] - max.b.alleles
      decrements.b  <- seq(from = max.b.alleles, to = 0, by = -1)
      res[[i]]     <- list()
      for (n in 1:length(decrements.b)) {
         A.i <- (max.a.alleles + decrements.b[n])
         B.i <- (max.b.alleles - decrements.b[n])
         BAF <- (fract.normal.alleles + (cellularity * B.i)) / ((alleles[i] * cellularity) + (CNn * fract.normal.alleles))
         res[[i]][[n]] <- cbind(A = A.i, B = B.i, BAF = BAF, CNt = alleles[i])
      }
   }
   for (i in 1:length(res)) {
      res[[i]] <- do.call(rbind,res[[i]])
   }
   as.data.frame(do.call(rbind,res))
}
dt2 <- function(x, df, ncp, log = FALSE, mean, sd) {
   x2 <- (x - mean) / sd
   dt(x2, df = df, ncp = ncp, log = log)
}
baf.bayes <- function(Bf, depth.ratio, cellularity, ploidy, avg.depth.ratio,
                      sd.Bf = 0.1, sd.ratio = 0.5, weight.Bf = 1, weight.ratio = 1, CNt.min = 0,
                      CNt.max = 7, CNn = 2, priors.table = data.frame(CN = CNt.min:CNt.max,
                      value = 1), ratio.priority = FALSE) {

   mufreq.tab <- data.frame(Bf = Bf, ratio = log(depth.ratio),
                            sd.Bf = sd.Bf, sd.ratio = sd.ratio,
                            weight.Bf = weight.Bf,
                            weight.ratio = weight.ratio)
   mufreq.depth.ratio <- model.points(cellularity = cellularity, ploidy = ploidy,
                                      types = cbind(CNn = CNn, CNt = CNt.min:CNt.max, Mt = 0),
                                      avg.depth.ratio = avg.depth.ratio)
   model.d.ratio      <- cbind(CNt = CNt.min:CNt.max, depth.ratio = log(mufreq.depth.ratio[, 2]))
   model.baf          <- expected.baf(sd = mean(sd.Bf[Bf > 0], na.rm = TRUE), CNn = CNn, CNt = CNt.max, cellularity = cellularity)
   #model.baf          <- theoretical.baf(CNn = CNn, CNt = CNt.max, cellularity = cellularity)
   # B-allele freq are never 0.5, always smaller. just a work around on this.. to be better fixed!
   #model.baf$BAF[model.baf$A==model.baf$B] <- quantile(rep(mufreq.tab$Bf, times = mufreq.tab$N.Bf),
   #                                                    na.rm = TRUE, probs = 0.95)
   if(CNt.min == 0) {
     model.baf          <- as.data.frame(rbind(c(0, 0, max(model.baf$BAF), 0), model.baf))
   }
   #                                                na.rm = TRUE, probs = skew.baf)
   model.pts          <- merge(model.baf, model.d.ratio)
   # model.pts          <- cbind(baf.type = apply(model.pts[, 1:3], 1, FUN = function(x) paste(x, collapse = "_")),
   #                             model.pts[, 4:5])
   rows.x             <- 1:nrow(mufreq.tab)

   priors <- rep(1, nrow(model.pts))
   for (i in 1:nrow(priors.table)) {
      priors[model.pts$CNt == priors.table$CN[i]] <- priors.table$value[i]
   }
   priors <- priors / sum(priors)

   bayes.fit <- function (x, mat, model.pts, priors, ratio.priority) {
      test.ratio <- model.pts$depth.ratio
      test.baf   <- model.pts$BAF
      min.offset <- 1e-323
      #score.r    <- depth.ratio.dbinom(size = mat[x,]$sd.ratio, depth.ratio = mat[x,]$ratio, test.ratio)
      score.r    <- dt2(sd = mat[x,]$sd.ratio/sqrt(mat[x,]$weight.ratio), mean = mat[x,]$ratio, x = test.ratio, df = 5, log = TRUE)
      score.r    <- score.r + log(priors)
      if (!is.na(mat[x,]$Bf) & !is.na(mat[x,]$sd.Bf/sqrt(mat[x,]$weight.Bf))) {
         #score.b    <- baf.dpois(baf = mat[x,]$Bf, depth.t = mat[x,]$N.Bf, test.baf, log = TRUE)
         score.b    <- dt2(mean = mat[x,]$Bf, sd = mat[x,]$sd.Bf/sqrt(mat[x,]$weight.Bf), x = test.baf, df = 5, log = TRUE)
         post.model <- score.r + score.b
      } else {
         post.model <- score.r         
      }

      #post.model[post.model == 0] <- min.offset
      if (ratio.priority == FALSE) {
         max.lik <-  which.max(post.model)
         max.post <- c(as.numeric(model.pts[max.lik,1:3]), post.model[max.lik])
      } else {
         res.cn     <- model.pts$CNt[which.max(score.r)]
         idx.pts    <- model.pts$CNt == res.cn
         model.lik  <- cbind(model.pts[idx.pts, 1:3], post.model[idx.pts])
         if (is.null(dim(model.lik))) {
            max.post <- model.lik
         } else {
            max.post   <- model.lik[which.max(model.lik[,4]),]
         }
      }
      if (is.na(mat[x,]$Bf)) {
         max.post[2:3] <- NA
      }
      max.post
   }
   bafs.L           <- mapply(FUN = bayes.fit, rows.x,
                         MoreArgs = list(mat = mufreq.tab,
                                         model.pts = model.pts,
                                         priors = priors,
                                         ratio.priority = ratio.priority),
                                         SIMPLIFY = FALSE)
   bafs.L           <- do.call(rbind, bafs.L)
   colnames(bafs.L) <- c("CNt", "A", "B", "LPP")
   bafs.L
}
baf.model.fit<-function (cellularity = seq(0.3, 1, by = 0.01), ploidy = seq(1, 
    7, by = 0.1), mc.cores = getOption("mc.cores", 2L), ...) 
{
    result <- expand.grid(ploidy = ploidy, cellularity = cellularity, 
        KEEP.OUT.ATTRS = FALSE)
    fit.cp <- function(ii) {
        L.model <- baf.bayes(cellularity = result$cellularity[ii], 
            ploidy = result$ploidy[ii], ...)
        sum(L.model[, 4])
    }
    bayes.res <- mclapplyPb(X = 1:nrow(result), FUN = fit.cp, 
        mc.cores = mc.cores)
    result$LPP <- unlist(bayes.res)
    z <- tapply(result$LPP, list(result$ploidy, result$cellularity), 
        mean)
    x <- as.numeric(rownames(z))
    y <- as.numeric(colnames(z))
    max.lik <- max(result$LPP, na.rm = TRUE)
    LogSumLik <- log(sum(exp(result$LPP - max.lik))) + max.lik
    znorm <- exp(z - LogSumLik)
    list(ploidy = x, cellularity = y, lpp = znorm)
}
mufreq.model.fit<-function (cellularity = seq(0.3, 1, by = 0.01), ploidy = seq(1, 
    7, by = 0.1), mc.cores = getOption("mc.cores", 2L), ...) 
{
    result <- expand.grid(ploidy = ploidy, cellularity = cellularity, 
        KEEP.OUT.ATTRS = FALSE)
    fit.cp <- function(ii) {
        L.model <- mufreq.bayes(cellularity = result$cellularity[ii], 
            ploidy = result$ploidy[ii], ...)
        sum(L.model[, 4])
    }
    bayes.res <- mclapplyPb(X = 1:nrow(result), FUN = fit.cp, 
        mc.cores = mc.cores)
    result$LPP <- unlist(bayes.res)
    z <- tapply(result$LPP, list(result$ploidy, result$cellularity), 
        mean)
    x <- as.numeric(rownames(z))
    y <- as.numeric(colnames(z))
    max.lik <- max(result$LPP, na.rm = TRUE)
    LogSumLik <- log(sum(exp(result$LPP - max.lik))) + max.lik
    znorm <- exp(z - LogSumLik)
    list(ploidy = x, cellularity = y, lpp = znorm)
}

sequenza.fit<-function (sequenza.extract, female = TRUE, N.ratio.filter = 10, 
        N.BAF.filter = 1, segment.filter = 3e+06, mufreq.treshold = 0.1, 
        XY = c(X = "X", Y = "Y"), cellularity = seq(0.1, 1, 0.01), 
        ploidy = seq(1, 7, 0.1), ratio.priority = FALSE, method = "baf", 
        priors.table = data.frame(CN = 2, value = 2), chromosome.list = 1:24, 
        mc.cores = getOption("mc.cores", 2L)) 
{
    if (method == "baf") {
        if (is.null(chromosome.list)) {
            segs.all <- do.call(rbind, sequenza.extract$segments)
        }else {
            segs.all <- do.call(rbind, sequenza.extract$segments[chromosome.list])
        }
        segs.len <- segs.all$end.pos - segs.all$start.pos
        avg.depth.ratio <- 1
        avg.sd.ratio <- sum(segs.all$sd.ratio * segs.all$N.ratio, 
                na.rm = TRUE)/sum(segs.all$N.ratio, na.rm = TRUE)
        avg.sd.Bf <- sum(segs.all$sd.BAF * segs.all$N.BAF, na.rm = TRUE)/sum(segs.all$N.BAF, 
                na.rm = TRUE)
        segs.all$sd.BAF[segs.all$sd.BAF == 0] <- max(segs.all$sd.BAF, 
                na.rm = TRUE)
        segs.all$sd.ratio[segs.all$sd.ratio == 0] <- max(segs.all$sd.ratio, 
                na.rm = TRUE)
        segs.filt <- segs.all$N.ratio > N.ratio.filter & segs.all$N.BAF > 
                N.BAF.filter
        segs.filt <- segs.len >= segment.filter & segs.filt
        if (female) {
            segs.is.xy <- segs.all$chromosome == XY["Y"]
        }else {
            segs.is.xy <- segs.all$chromosome %in% XY
        }
        filt.test <- segs.filt & !segs.is.xy
        seg.test <- segs.all[filt.test, ]
        seg.len.mb <- segs.len[filt.test]/1e+06
        baf.model.fit(Bf = seg.test$Bf, depth.ratio = seg.test$depth.ratio, 
                sd.ratio = seg.test$sd.ratio, weight.ratio = seg.len.mb, 
                sd.Bf = seg.test$sd.BAF, weight.Bf = seg.len.mb, 
                avg.depth.ratio = avg.depth.ratio, cellularity = cellularity, 
                ploidy = ploidy, priors.table = priors.table, mc.cores = mc.cores, 
                ratio.priority = ratio.priority)
    }else if (method == "mufreq") {
        if (is.null(chromosome.list)) {
            mut.all <- do.call(rbind, sequenza.extract$mutations)
            mut.all <- na.exclude(mut.all)
            segs.all <- do.call(rbind, sequenza.extract$segments)
        }
        else {
            mut.all <- do.call(rbind, sequenza.extract$mutations[chromosome.list])
            mut.all <- na.exclude(mut.all)
            segs.all <- do.call(rbind, sequenza.extract$segments[chromosome.list])
        }
        mut.filt <- mut.all$F >= mufreq.treshold
        segs.len <- segs.all$end.pos - segs.all$start.pos
        avg.depth.ratio <- 1
        if (female) {
            mut.is.xy <- mut.all$chromosome == XY["Y"]
        }
        else {
            mut.is.xy <- mut.all$chromosome %in% XY
        }
        filt.test <- mut.filt & !mut.is.xy
        mut.test <- mut.all[filt.test, ]
        w.mufreq <- round(mut.test$good.reads, 0)
        mufreq.model.fit(mufreq = mut.test$F, depth.ratio = mut.test$adjusted.ratio, 
                weight.ratio = 2 * w.mufreq, weight.mufreq = w.mufreq, 
                avg.depth.ratio = avg.depth.ratio, cellularity = cellularity, 
                ploidy = ploidy, priors.table = priors.table, mc.cores = mc.cores)
    }
    else {
        stop("The only available methods are \"baf\" and \"mufreq\"")
    }
}
cp.plot.contours <- function(cp.table, likThresh = c(0.95), alternative = TRUE,
                             col = palette(), legend.pos = 'bottomright', pch = 18, alt.pch = 3, ...) {
   znormsort <- sort(cp.table$lpp, decreasing = TRUE)
   znormcumLik <- cumsum(znormsort)
   n <- sapply(likThresh, function(x) sum(znormcumLik < x) + 1)
   LikThresh <- znormsort[n]
   names(LikThresh) <- paste0(likThresh * 100, '%')
   contour(x = cp.table$ploidy, y = cp.table$cellularity, z = cp.table$lpp,
           levels = znormsort[n], col = col, drawlabels = FALSE,
           xlab = "Ploidy", ylab = "Cellularity", ...)
   max.xy <- which(cp.table$lpp == max(cp.table$lpp), arr.ind = TRUE)
   points(x = cp.table$ploidy[max.xy[, "row"]],
          y = cp.table$cellularity[max.xy[, "col"]], pch = pch)
   if (alternative == TRUE){
      alt.sol <- alternative.cp.solutions(cp.table)
      alt.sol <- alt.sol[-1, ]
      points(x = alt.sol$ploidy,
             y = alt.sol$cellularity, pch = alt.pch)
   }
   if(!is.na(legend.pos)) {
      if (alternative == FALSE) {
         legend(legend.pos, legend = c(paste("C.R.", names(LikThresh), sep = " "), "Point estimate"),
                col = c(col[1:length(LikThresh)], "black"), lty = c(rep(1, length(LikThresh)), NA),
                pch = c(rep(NA, length(LikThresh)), pch), border = NA, bty = "n")
      } else {
         legend(legend.pos, legend = c(paste("C.R.", names(LikThresh), sep = " "),
                                       "Point estimate", "Alternative solutions"),
                col = c(col[1:length(LikThresh)], "black", "black"), lty = c(rep(1, length(LikThresh)), NA, NA),
                pch = c(rep(NA, length(LikThresh)), pch, alt.pch), border = NA, bty = "n")         
      }
   }
   invisible(LikThresh)
}
alternative.cp.solutions <- function(cp.table) {
   ci <- get.ci(cp.table)
   p.alt <- which(diff(sign(diff(ci$values.ploidy$y))) == -2) + 1
   get.alt <- function(idx.p, cp.table) {
      idx.c <- which.max(cp.table$lpp[idx.p,])
      c(cellularity = cp.table$cellularity[idx.c],
        ploidy = cp.table$ploidy[idx.p],
        SLPP = cp.table$lpp[idx.p, idx.c])
   }
   res <- lapply(p.alt, FUN = function (x) get.alt(x, cp.table))
   res <- as.data.frame(do.call(rbind, res))
   if (nrow(res) > 0 ){
      res[order(res$SLPP, decreasing = TRUE), ]
   } else {
      data.frame(cellularity = ci$max.cellularity, 
                 ploidy = ci$max.ploidy,
                 SLPP =  cp.table$lpp[which(cp.table$ploidy == ci$max.ploidy),
                                      which(cp.table$cellularity == ci$max.cellularity)])
   }
}

get.ci <- function(cp.table, level = 0.95) {
  znormsort <- sort(cp.table$lpp, decreasing = TRUE)
  znormcumLik <- cumsum(znormsort)
  n <- sapply(level, function(x) sum(znormcumLik < x) + 1)
  LikThresh <- znormsort[n]
  values.x <- data.frame(x = cp.table$ploidy, y = apply(cp.table$lpp, 1, max))
  values.y <- data.frame(x = apply(cp.table$lpp, 2, max), y = cp.table$cellularity)
  up.x  <- max(values.x$x[values.x$y >= LikThresh])
  low.x <- min(values.x$x[values.x$y >= LikThresh])
  max.x <- values.x$x[which.max(values.x$y)]
  up.y  <- max(values.y$y[values.y$x >= LikThresh])
  low.y <- min(values.y$y[values.y$x >= LikThresh])
  max.y <- values.y$y[which.max(values.y$x)]
  results <- list()
  values.x$y <- values.x$y/sum(values.x$y)
  values.y$x <- values.y$x/sum(values.y$x)
  results$values.ploidy <- values.x
  results$confint.ploidy <- c(low.x, up.x)
  results$max.ploidy <- max.x
  results$values.cellularity <- values.y
  results$confint.cellularity <- c(low.y, up.y)
  results$max.cellularity <- max.y
  results
}
