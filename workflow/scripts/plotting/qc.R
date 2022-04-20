#
# Copyright (C) 2017 Sascha Meiers
# Distributed under the MIT software license, see the accompanying
# file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.



suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(cowplot))


add_overview_plot <- T


args <- commandArgs(trailingOnly = T)
if (length(args) < 2 || length(args) > 4 || !grepl("\\.pdf$", args[length(args)]) || any(!file.exists(args[1:(length(args) - 1)]))) {
    warning("Usage: Rscript R/qc.R input-file [SCE-file] [cell-info-file] output-pdf")
    quit(status = 1)
}
f_in <- args[1]
pdf_out <- args[length(args)]

# Detect info or SCE file.
info <- NULL
sces <- NULL
is_sce_file <- function(x) {
    all(c("sample", "cell", "chrom", "start", "end", "state") %in% colnames(x))
}
is_info_file <- function(x) {
    all(c("sample", "cell", "pass1", "dupl", "mapped", "nb_p", "nb_r", "nb_a") %in% colnames(x))
}
if (length(args) > 2) {
    x <- fread(args[2])
    if (is_sce_file(x)) {
        message("* Using SCE file ", args[2])
        sces <- x
    } else if (is_info_file(x)) {
        message("* Using INFO file ", args[2])
        info <- x
    }
    if (length(args) > 3) {
        x <- fread(args[3])
        if (is_sce_file(x)) {
            message("* Using SCE file ", args[3])
            sces <- x
        } else if (is_info_file(x)) {
            message("* Using INFO file ", args[3])
            info <- x
        }
    }
}
if (!is.null(sces)) {
    sces[, chrom := sub("^chr", "", chrom)]
    sces[, chrom := factor(chrom, levels = as.character(c(1:22, "X", "Y")), ordered = T)]
}



format_Mb <- function(x) {
    paste(comma(x / 1e6), "Mb")
}


# if gzip
zcat_command <- "zcat"
if (substr(f_in, nchar(f_in) - 2, nchar(f_in)) == ".gz") {
    f_in <- paste(zcat_command, f_in)
}

# Read counts & filter chromosomes (this is human-specific)
d <- fread(f_in)

# Check that correct files are given:
invisible(assert_that(
    "chrom" %in% colnames(d),
    "start" %in% colnames(d) && is.integer(d$start),
    "end" %in% colnames(d) && is.integer(d$end),
    "sample" %in% colnames(d),
    "cell" %in% colnames(d),
    "w" %in% colnames(d) && is.numeric(d$w),
    "c" %in% colnames(d) && is.numeric(d$c),
    "class" %in% colnames(d)
))

# Re-name and -order chromosomes - this is human-specific
d <- d[, chrom := sub("^chr", "", chrom)][]
d <- d[grepl("^([1-9]|[12][0-9]|X|Y)$", chrom), ]
d <- d[, chrom := factor(chrom, levels = as.character(c(1:22, "X", "Y")), ordered = T)]


message("* Writing plot ", pdf_out)

if (add_overview_plot == T) {
    message("* Plotting an overview page")

    n_samples <- nrow(unique(d[, .(sample)]))
    n_cells <- nrow(unique(d[, .(sample, cell)]))
    n_bins <- nrow(unique(d[, .(chrom, start, end)]))
    mean_bin <- unique(d[, .(chrom, start, end)])[, mean(end - start)]
    n_excl <- nrow(d[, .N, by = .(chrom, start, end, class)][class == "None" & N == n_cells, ])

    # Bin sizes
    ov_binsizes <- ggplot(unique(d[, .(chrom, start, end)])) +
        geom_histogram(aes(end - start), bins = 50) +
        theme_minimal() +
        scale_y_log10(breaks = c(1, 10, 100, 1000, 10e3, 100e3, 1e6)) +
        scale_x_log10(labels = comma) +
        ggtitle(paste0("Bin sizes (", n_bins, " bins, mean ", round(mean_bin / 1000, 1), " kb)")) +
        xlab("Bin size (bp)")

    # analyse how many bins are "None"
    ov_excbins <- ggplot(d[, .N, by = .(chrom, start, end, class)][class == "None" & N == n_cells, ]) +
        aes(chrom) +
        geom_bar() +
        theme_minimal() +
        ggtitle(paste0("Excluded bins per chromosome (total = ", n_excl, ")"))


    # coverage
    ov_coverage <- ggplot(d[, .(total = sum(w + c)), by = .(sample, cell)]) +
        geom_histogram(aes(total, fill = sample), bins = 50) +
        scale_x_continuous(
            breaks = pretty_breaks(5),
            labels = comma
        ) +
        xlab("Total number of reads per cell") +
        theme_minimal() +
        theme(legend.position = "bottom") +
        scale_fill_brewer(type = "qual", palette = 6)

    # Overview mean / variance
    d_mv <- d[class != "None", .(mean = mean(w + c), var = var(w + c)), by = .(sample, cell)]
    d_p <- d_mv[, .(p = sum(mean * mean) / sum(mean * var)), by = sample]
    ov_meanvar <- ggplot(d_mv) +
        geom_point(aes(mean, var), alpha = 0.4) +
        facet_wrap(~sample, nrow = 1) +
        theme_minimal() +
        geom_abline(data = d_p, aes(slope = 1 / p, intercept = 0), col = "dodgerblue") +
        geom_label(data = d_p, aes(x = 0, y = Inf, label = paste("p =", round(p, 3))), hjust = 0, vjust = 1) +
        ggtitle("Mean variance relationship of reads per bin") +
        xlab("Mean") +
        ylab("Variance")


    # Arranging overview plot
    content <- ggdraw() +
        draw_plot(ov_binsizes, x = 0, y = .66, width = .5, height = .33) +
        draw_plot(ov_excbins, x = .5, y = .66, width = .5, height = .33) +
        draw_plot(ov_coverage, x = 0, y = .33, width = .5, height = .33) +
        draw_plot(ov_meanvar, x = 0, y = 0, width = min(n_samples / 3, 1), height = .33)


    # Add duplicate rates if available
    if (exists("info")) {
        ov_duplicate <- ggplot(info) +
            aes(dupl / (mapped - suppl), fill = sample) +
            geom_histogram(bins = 50) +
            xlab("Duplicate rate") +
            scale_x_continuous(labels = percent) +
            theme_minimal() +
            theme(legend.position = "bottom") +
            scale_fill_brewer(type = "qual", palette = 6)
        content <- content +
            draw_plot(ov_duplicate, x = 0.55, y = .33, width = .45, height = .33)
    }

    title <- ggdraw() + draw_label(paste("Overview across", n_cells, "cells from", n_samples, "samples"), fontface = "bold")
    side <- ggdraw() + draw_label(label = paste0(args[1], "\n", date()), angle = 90, size = 10, vjust = 1)

    final <- plot_grid(title, content, ncol = 1, rel_heights = c(0.07, 1))
    xxx <- plot_grid(side, final, nrow = 1, rel_widths = c(0.05, 1))
}

cairo_pdf(pdf_out, width = 14, height = 10, onefile = T)
if (add_overview_plot == T) {
    print(xxx)
}


# Plot all cells
for (s in unique(d$sample))
{
    for (ce in unique(d[sample == s, ]$cell))
    {
        message(paste("* Plotting sample", s, "cell", ce))

        e <- d[sample == s & cell == ce, ]


        # Calculate some information
        info_binwidth <- median(e$end - e$start)
        info_reads_per_bin <- median(e$w + e$c)
        if (!is.integer(info_reads_per_bin)) info_reads_per_bin <- round(info_reads_per_bin, 2)
        info_chrom_sizes <- e[, .(xend = max(end)), by = chrom]
        info_num_bins <- nrow(e)
        info_total_reads <- sum(e$c + e$w)
        info_y_limit <- 2 * info_reads_per_bin + 1
        if (!is.integer(info_y_limit)) info_y_limit <- round(info_y_limit, 2)
        info_sample_name <- substr(s, 1, 25)
        if (nchar(s) > 25) info_sample_name <- paste0(info_sample_name, "...")
        info_cell_name <- substr(ce, 1, 25)
        if (nchar(ce) > 25) info_cell_name <- paste0(info_cell_name, "...")

        # start main plot:
        plt <- ggplot(e) +
            aes(x = (start + end) / 2)


        # prepare consecutive rectangles for a better plotting experience
        consecutive <- cumsum(c(0, abs(diff(as.numeric(as.factor(e$class))))))
        e$consecutive <- consecutive
        f <- e[, .(start = min(start), end = max(end), class = class[1]), by = .(consecutive, chrom)][]

        plt <- plt +
            geom_rect(data = f, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = class), inherit.aes = F, alpha = 0.2) +
            scale_fill_manual(values = c(WW = "sandybrown", CC = "paleturquoise4", WC = "yellow", None = NA))

        # Show SCEs
        if (!is.null(sces)) {
            sces_local <- sces[sample == s & cell == ce][, .SD[.N > 1], by = chrom]
            if (nrow(sces_local) > 0) {
                sces_local <- sces_local[, .(pos = (end[1:(.N - 1)] + start[2:(.N)]) / 2), by = chrom]
                plt <- plt + geom_point(data = sces_local, aes(x = pos, y = -info_y_limit), size = 3, shape = 18)
            }
        }


        # Watson/Crick bars
        plt <- plt +
            geom_rect(aes(xmin = start, xmax = end, ymin = -w, ymax = 0), fill = "sandybrown") +
            geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = c), fill = "paleturquoise4") +
            # geom_bar(aes(y = -w, width=(end-start)), stat='identity', position = 'identity', fill='sandybrown') +
            # geom_bar(aes(y = c, width=(end-start)), stat='identity', position = 'identity', fill='paleturquoise4') +
            # Trim image to 2*median cov
            coord_flip(expand = F, ylim = c(-info_y_limit, info_y_limit)) +
            facet_grid(. ~ chrom, switch = "x") +
            ylab("Watson | Crick") + xlab(NULL) +
            scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
            scale_y_continuous(breaks = pretty_breaks(3)) +
            theme_classic() +
            theme(
                panel.spacing = unit(0.2, "lines"),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                strip.background = element_rect(fill = NA, colour = NA)
            ) +
            guides(fill = FALSE) +
            # Dotted lines at median bin count
            geom_segment(
                data = info_chrom_sizes, aes(xend = xend, x = 0, y = -info_reads_per_bin, yend = -info_reads_per_bin),
                linetype = "dotted", col = "darkgrey", size = 0.5
            ) +
            geom_segment(
                data = info_chrom_sizes, aes(xend = xend, x = 0, y = +info_reads_per_bin, yend = +info_reads_per_bin),
                linetype = "dotted", col = "darkgrey", size = 0.5
            ) +
            geom_segment(data = info_chrom_sizes, aes(xend = xend, x = 0), y = 0, yend = 0, size = 0.5)

        # Rename classes:
        labels <- e[, .N, by = class][, label := paste0(class, " (n=", N, ")")][]

        e[, class := factor(class, levels = labels$class, labels = labels$label)]

        # Histogram in upper right corner
        e.melt <- melt.data.table(e, c("chrom", "start", "end", "class"), measure.vars = c("w", "c"), variable.name = "strand", value.name = "coverage")
        plt_hist_xlim <- 10 + 3 * info_reads_per_bin
        plt_hist <- ggplot(e.melt) +
            aes(coverage, fill = strand) +
            geom_histogram(binwidth = 1, position = position_dodge(), alpha = 0.9) +
            scale_x_continuous(limits = c(-1, plt_hist_xlim), breaks = pretty_breaks(5), labels = comma) +
            theme(text = element_text(size = 10), axis.text = element_text(size = 8)) +
            scale_fill_manual(values = c(w = "sandybrown", c = "paleturquoise4")) +
            guides(fill = FALSE, col = FALSE) +
            ylab("bin count") +
            xlab("reads per bin") +
            facet_wrap(~class, nrow = 1, scales = "free")
        if (!is.null(info)) {
            Ie <- info[sample == s & cell == ce, ]
            if (nrow(Ie) != 1 || Ie$pass1 != 1 || !all(c("WW", "WC", "CC") %in% unique(Ie$class))) {
                message("  Problem finding additional info for ", s, " - ", ce)
            } else {
                p <- Ie$nb_p
                r <- Ie$nb_r
                a <- Ie$nb_a
                x <- seq(0, plt_hist_xlim)
                scale_factors <- e.melt[, .N, by = .(class, strand, coverage)][, .(scale = max(N)), by = class]
                nb <- data.table(
                    x = rep(x, 6),
                    strand = rep(c(rep("w", length(x)), rep("c", length(x))), 3),
                    class = c(
                        rep(labels[class == "WW", ]$label, 2 * length(x)),
                        rep(labels[class == "WC", ]$label, 2 * length(x)),
                        rep(labels[class == "CC", ]$label, 2 * length(x))
                    ),
                    scale = c(
                        rep(scale_factors[class == labels[class == "WW", ]$label, ]$scale, 2 * length(x)),
                        rep(scale_factors[class == labels[class == "WC", ]$label, ]$scale, 2 * length(x)),
                        rep(scale_factors[class == labels[class == "CC", ]$label, ]$scale, 2 * length(x))
                    ),
                    y = c(
                        dnbinom(x, a * r, p),
                        dnbinom(x, (1 - a) * r, p),
                        dnbinom(x, r / 2, p),
                        dnbinom(x, r / 2, p),
                        dnbinom(x, (1 - a) * r, p),
                        dnbinom(x, a * r, p)
                    )
                )
                plt_hist <- plt_hist + geom_line(data = nb, aes(x, y * scale, col = strand))
            }
        }


        plot_hst_width <- .03 + .13 * length(unique(e$class))
        all <- ggdraw() + draw_plot(plt) +
            draw_plot(plt_hist, x = .45, y = .76, width = plot_hst_width, height = .23) +
            draw_label(paste("Sample:", info_sample_name), x = .29, y = .97, vjust = 1, hjust = 0, size = 14) +
            draw_label(paste("Cell:", info_cell_name), x = .29, y = .94, vjust = 1, hjust = 0, size = 12) +
            draw_label(paste("Median binwidth:", format(round(info_binwidth / 1000, 0), big.mark = ",", scientific = F), "kb"),
                x = .29, y = .91, vjust = 1, hjust = 0, size = 10
            ) +
            draw_label(paste("Number bins:", format(info_num_bins, big.mark = ",")),
                x = .29, y = .89, vjust = 1, hjust = 0, size = 10
            ) +
            draw_label(paste("Total number of reads:", format(info_total_reads, big.mark = ",")),
                x = .29, y = .87, vjust = 1, hjust = 0, size = 10
            ) +
            draw_label(paste("Median reads/bin (dotted):", info_reads_per_bin),
                x = .29, y = .85, vjust = 1, hjust = 0, size = 10
            ) +
            draw_label(paste0("Plot limits: [-", info_y_limit, ",", info_y_limit, "]"),
                x = .29, y = .83, vjust = 1, hjust = 0, size = 10
            )
        # If available, add additional info like duplicate rate and NB params!
        if (!is.null(info)) {
            Ie <- info[sample == s & cell == ce, ]
            if (nrow(Ie) == 1) {
                all <- all +
                    draw_label(paste0("Duplicate rate: ", round(Ie$dupl / Ie$mapped, 2) * 100, "%"),
                        x = .29, y = .80, vjust = 1, hjust = 0, size = 10
                    )
                if (Ie$pass1 == 1) {
                    all <- all +
                        draw_label(paste0("NB parameters (p,r,a): ", round(Ie$nb_p, 2), ",", round(Ie$nb_r, 2), ",", round(Ie$nb_a, 2)),
                            x = .29, y = .78, vjust = 1, hjust = 0, size = 10
                        )
                }
            }
        }

        # If available, write number of SCEs detected
        if (!is.null(sces)) {
            sces_local <- sces[sample == s & cell == ce][, .SD[.N > 1], by = chrom]
            if (nrow(sces_local) > 0) sces_local <- sces_local[, .(pos = (end[1:(.N - 1)] + start[2:(.N)]) / 2), by = chrom]
            all <- all + draw_label(paste("SCEs detected:", nrow(sces_local)),
                x = .29, y = .76, vjust = 1, hjust = 0, size = 10
            )
        }

        print(all)
    }
}
# dev.off()