#!/usr/bin/env Rscript

# Usage: Rscript copywriter2seg.R <in_fname.Rdata> <out_fname.seg>

(args = commandArgs(TRUE))

in_fname = args[1]
out_fname = args[2]
# in_fname = "segment.Rdata"
# out_fname = "sample.seg"

cat("Loading", in_fname, "\n")
load(in_fname)
seg = segment.CNA.object$output
seg$loc.start = floor(seg$loc.start)
seg$loc.end = floor(seg$loc.end)
write.table(seg, out_fname,
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
cat("Wrote", out_fname, "\n")
