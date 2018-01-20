# Benchmark CNVkit-RNA estimates vs. TCGA segments

# TCGA-generated aCGH segments
tcga_acgh_seg = $(sort $(wildcard tcga-acgh-seg/TCGA-*.txt))
acgh_cns := $(patsubst tcga-acgh-seg/%.txt,acgh/%.acgh.cns,$(tcga_acgh_seg))

# RNA-derived read counts, expression ratios
tcga_rna_counts = $(sort $(wildcard tcga-rna-counts/TCGA-*.txt))
plain_cnr := $(patsubst tcga-rna-counts/%.txt,out/%.cnr,$(tcga_rna_counts))

# Smoothings
tsmooth100_cnr := $(patsubst out/%.cnr,tsmooth100/%.tsmooth100.cnr,$(plain_cnr))
wsmooth200_cnr := $(patsubst out/%.cnr,wsmooth200/%.wsmooth200.cnr,$(plain_cnr))
wsmooth400_cnr := $(patsubst out/%.cnr,wsmooth400/%.wsmooth400.cnr,$(plain_cnr))
# asmooth50_cnr := asmooth50/%

# Segmentations
arm_cns := $(patsubst out/%.cnr,arm/%.cns,$(plain_cnr))
cbs_cns := $(patsubst out/%.cnr,cbs/%.cns,$(plain_cnr))
flasso_cns := $(patsubst out/%.cnr,flasso/%.cns,$(plain_cnr))
haar_cns := $(patsubst out/%.cnr,haar/%.cns,$(plain_cnr))
hmm_cns := $(patsubst out/%.cnr,hmm/%.cns,$(plain_cnr))

# Bias corrections
# TODO create .cnr files without GC and/or tx_length corrections
nogc_cnr :=
notx_cnr :=
nogctx_cnr :=

# Plots
all_labels_bin := plain tsmooth100 wsmooth200 wsmooth400
all_labels_seg := arm cbs cbsgt20 flasso haar hmm hmmgt20
resid_labels := tsmooth100 wsmooth200 arm cbs flasso haar hmm
# bias_labels := plain nogc notx nogctx
bias_labels := plain wsmooth200
heat_labels := arm cbs flasso haar hmm

plots_2d := $(patsubst %,plots/rna-vs-acgh-%.pdf,$(resid_labels))
plot_resid := plots/residuals-vs-acgh.pdf
plot_bias := plots/bias-vs-acgh.pdf
heatmaps :=  $(patsubst %,plots/heatmap-%.pdf,$(heat_labels))


all: build/_thunk $(heatmaps) $(plots_2d) $(plot_resid) $(plot_bias)

.PSEUDO: clean
clean:
	rm -vf build/*


# Generate intermediate intermediate files that aren't explicit targets
build/_thunk: build/tcga_combined_log2ratios.tsv \
	$(firstword $(acgh_cns)) build/tcga-acgh-seg-genes.tsv
	touch $@


# Convert TCGA raw data to CNVkit formats ===

correlate := ~/code/cnvkit/scripts/cnv_expression_correlate.py
corr_out := tcga-skcm.cnv-expr-corr.tsv

$(corr_out): $(correlate) resource/TCGA_SKCM_CNA.tsv resource/TCGA_SKCM_RNA_Seq_v2_expression_median.tsv
	python $^ -o $@


gene_info := ~/code/cnvkit/data/ensembl-gene-info.hg38.tsv

# Makes all $(plain_cnr) as out/*.cnr
build/tcga_combined_log2ratios.tsv: $(corr_out) $(gene_info)
	mkdir -p build/
	cnvkit.py import-rna -f counts -c $< -g $(gene_info) \
		-o $@ --output-dir out/ $(tcga_rna_counts)

# Makes all $(acgh_cns) as acgh/*.acgh.cns
$(firstword $(acgh_cns)): tcga2cns.py $(tcga_acgh_seg)
	mkdir -p acgh/
	python $^ -d acgh/ -g $(gene_info)


# CNVkit smoothing and segmentation ===

# $(asmooth50_cnr): adaptive_smooth.py $(plain_cnr)
#         python $^ -w 50 -d asmooth50

$(tsmooth100_cnr): smooth_tirosh.py $(plain_cnr)
	mkdir -p tsmooth100
	python $^ -w 100 -d tsmooth100

$(wsmooth200_cnr): smooth_cnr.py $(plain_cnr)
	mkdir -p wsmooth200
	python $^ -w 200 -d wsmooth200

$(wsmooth400_cnr): smooth_cnr.py $(plain_cnr)
	mkdir -p wsmooth400
	python $^ -w 400 -d wsmooth400

$(arm_cns): arm/%.cns: out/%.cnr
	cnvkit.py segment -m none $< -o $@

$(cbs_cns): cbs/%.cns: out/%.cnr
	cnvkit.py segment -m cbs -p 2 $< -o $@

$(hmm_cns): hmm/%.cns: out/%.cnr
	cnvkit.py segment -m hmm $< -o $@

$(haar_cns): haar/%.cns: out/%.cnr
	cnvkit.py segment -m haar $< -o $@

$(flasso_cns): flasso/%.cns: out/%.cnr
	cnvkit.py segment -m flasso $< -o $@


# Intermediate data tables (used for plotting) ===

build/tcga-acgh-seg-genes.tsv: collate_by_gene.py $(acgh_cns)
	mkdir -p build/
	python $^ -s -o $@

build/tcga-rna-plain-genes.tsv: collate_by_gene.py $(plain_cnr)
	python $^ -o $@

build/tcga-rna-wsmooth200-genes.tsv: collate_by_gene.py $(wsmooth200_cnr)
	python $^ -o $@

build/tcga-rna-wsmooth400-genes.tsv: collate_by_gene.py $(wsmooth400_cnr)
	python $^ -o $@

# tcga-rna-asmooth50-genes.tsv: collate_by_gene.py $(asmooth50_cnr)
#         python $^ -o $@

build/tcga-rna-tsmooth100-genes.tsv: collate_by_gene.py $(tsmooth100_cnr)
	python $^ -o $@

build/tcga-rna-arm-genes.tsv: collate_by_gene.py $(arm_cns)
	python $^ -s -o $@

build/tcga-rna-cbs-genes.tsv: collate_by_gene.py $(cbs_cns)
	python $^ -s -o $@

build/tcga-rna-cbsgt20-genes.tsv: collate_by_gene.py $(cbs_cns)
	python $^ -s -w 20 -o $@

build/tcga-rna-hmm-genes.tsv: collate_by_gene.py $(hmm_cns)
	python $^ -s -o $@

build/tcga-rna-hmmgt20-genes.tsv: collate_by_gene.py $(hmm_cns)
	python $^ -s -w 20 -o $@

build/tcga-rna-flasso-genes.tsv: collate_by_gene.py $(flasso_cns)
	python $^ -s -o $@


# Plots ===
# - 2D scatter of aCGH vs. each RNA estimator (separately), with correlation
# - Violin or box plot of aCGH vs. each RNA estimator (together), with SD/IQR
# - Violin or box plot of each bias correction/not vs. aCGH (together), with SD/IQR

$(plots_2d): plots/rna-vs-acgh-%.pdf: \
	plot2d.py build/tcga-rna-%-genes.tsv build/tcga-acgh-seg-genes.tsv
	python $^ -o $@

$(plot_resid): \
	plot_residuals.py build/tcga-acgh-seg-genes.tsv \
	$(patsubst %,build/tcga-rna-%-genes.tsv,$(resid_labels))
	python $^ -o $@

$(plot_bias): \
	plot_residuals.py build/tcga-acgh-seg-genes.tsv \
	$(patsubst %,build/tcga-rna-%-genes.tsv,$(bias_labels))
	python $^ -o $@

# $(heatmaps): plots/heatmap-%.pdf: $(foreach $(heat_labels,%/...)
#         cnvkit.py heatmap -d --no-shift-xy $^ -o $@
plots/heatmap-arm.pdf: $(arm_cns)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

plots/heatmap-cbs.pdf: $(cbs_cns)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

plots/heatmap-flasso.pdf: $(flasso_cns)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

plots/heatmap-haar.pdf: $(haar_cns)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

plots/heatmap-hmm.pdf: $(hmm_cns)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

