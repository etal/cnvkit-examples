# Benchmark CNVkit-RNA estimates vs. TCGA segments

# TCGA-generated aCGH segments
tcga_acgh_seg = $(sort $(wildcard tcga-acgh-seg/TCGA-*.txt))
acgh_cns := $(patsubst tcga-acgh-seg/%.txt,acgh/%.acgh.cns,$(tcga_acgh_seg))

# RNA-derived read counts, expression ratios
tcga_rna_counts = $(sort $(wildcard tcga-rna-counts/TCGA-*.txt))
plain_cnr := $(patsubst tcga-rna-counts/%.txt,out/%.cnr,$(tcga_rna_counts))
nocorr_cnr := $(patsubst tcga-rna-counts/%.txt,nocorr/%.cnr,$(tcga_rna_counts))

# Smoothings
tsmooth100_cnr := $(patsubst out/%.cnr,tsmooth100/%.tsmooth100.cnr,$(plain_cnr))
tsmooth200_cnr := $(patsubst out/%.cnr,tsmooth200/%.tsmooth200.cnr,$(plain_cnr))
wsmooth100_cnr := $(patsubst out/%.cnr,wsmooth100/%.wsmooth100.cnr,$(plain_cnr))
wsmooth200_cnr := $(patsubst out/%.cnr,wsmooth200/%.wsmooth200.cnr,$(plain_cnr))
wsmooth400_cnr := $(patsubst out/%.cnr,wsmooth400/%.wsmooth400.cnr,$(plain_cnr))
nocorrsmooth100_cnr := $(patsubst nocorr/%.cnr,nocorrsmooth100/%.nocorrsmooth100.cnr,$(nocorr_cnr))
# asmooth50_cnr := asmooth50/%

# Segmentations
arm_cns := $(patsubst out/%.cnr,arm/%.cns,$(plain_cnr))
cbs_cns := $(patsubst out/%.cnr,cbs/%.cns,$(plain_cnr))
flasso_cns := $(patsubst out/%.cnr,flasso/%.cns,$(plain_cnr))
haar_cns := $(patsubst out/%.cnr,haar/%.cns,$(plain_cnr))
hmm_cns := $(patsubst out/%.cnr,hmm/%.cns,$(plain_cnr))

# Bias corrections
nogc_cnr := $(patsubst tcga-rna-counts/%.txt,nogc/%.cnr,$(tcga_rna_counts))
notx_cnr := $(patsubst tcga-rna-counts/%.txt,notx/%.cnr,$(tcga_rna_counts))
nogctx_cnr := $(patsubst tcga-rna-counts/%.txt,nogctx/%.cnr,$(tcga_rna_counts))

# Plots
all_labels_bin := plain tsmooth100 tsmooth200 nocorrsmooth100 wsmooth100 wsmooth200 wsmooth400
all_labels_seg := arm cbs cbsgt20 flasso haar hmm hmmgt20
resid_labels := tsmooth100 nocorrsmooth100 wsmooth100 arm cbs cbsgt20
bias_labels := nogctx nogc notx plain
heat_labels := arm cbs flasso

heatmaps :=  $(patsubst %,plots/heatmap-%.pdf,$(heat_labels))
plots_2d := $(patsubst %,plots/rna-vs-acgh-%.pdf,$(resid_labels))
plot_resid := plots/residuals-vs-acgh.pdf
plot_bias := plots/bias-vs-acgh.pdf
plot_bias_facet := plots/bias-vs-acgh-facet.pdf
sample_bias_gc := $(patsubst tcga-rna-counts/%.txt,samplebias/%.gc_bias.png,$(tcga_rna_counts))
sample_bias_tx := $(patsubst tcga-rna-counts/%.txt,samplebias/%.tx_length_bias.png,$(tcga_rna_counts))

all: build/_thunk $(heatmaps) \
	$(plots_2d) $(plot_resid) $(plot_bias) $(plot_bias_facet) \
	$(firstword $(sample_bias_gc)) $(firstword $(sample_bias_tx))

.PSEUDO: clean
clean:
	rm -vf build/*


# Generate intermediate intermediate files that aren't explicit targets
build/_thunk: build/tcga_combined_log2ratios.tsv \
	build/tcga_combined_log2ratios_nocorr.tsv \
	build/tcga_combined_log2ratios_nogc.tsv \
	build/tcga_combined_log2ratios_notx.tsv \
	build/tcga_combined_log2ratios_nogctx.tsv \
	$(firstword $(acgh_cns)) build/tcga-acgh-seg-genes.tsv
	touch $@


# Convert TCGA raw data to CNVkit formats ===

correlate := ~/code/cnvkit/scripts/cnv_expression_correlate.py
corr_out := tcga-skcm.cnv-expr-corr.tsv

$(corr_out): $(correlate) resource/TCGA_SKCM_CNA.tsv resource/TCGA_SKCM_RNA_Seq_v2_expression_median.tsv
	python $^ -o $@


gene_info := ~/code/cnvkit/data/ensembl-gene-info.hg38.tsv

# Makes all $(acgh_cns) as acgh/*.acgh.cns
$(firstword $(acgh_cns)): tcga2cns.py $(tcga_acgh_seg)
	mkdir -p $(dir $@)
	python $^ -d $(dir $@) -g $(gene_info)

# Makes all $(plain_cnr) as out/*.cnr
build/tcga_combined_log2ratios.tsv: $(corr_out) $(gene_info)
	mkdir -p $(dir $@)
	cnvkit.py import-rna -f counts -c $< -g $(gene_info) \
		-o $@ --output-dir out/ $(tcga_rna_counts)

# Makes nocorr_cnr (log2 values same as plain_cnr, but different weights)
build/tcga_combined_log2ratios_nocorr.tsv: $(gene_info)
	mkdir -p $(dir $@)
	cnvkit.py import-rna -f counts -g $(gene_info) \
		-o $@ --output-dir nocorr/ $(tcga_rna_counts)

# Makes nogc_cnr, notx_cnr, nogctx_cnr
build/tcga_combined_log2ratios_nogc.tsv: $(corr_out) $(gene_info)
	mkdir -p $(dir $@)
	cnvkit.py import-rna -f counts -c $< -g $(gene_info) \
		--no-gc \
		-o $@ --output-dir nogc/ $(tcga_rna_counts)

build/tcga_combined_log2ratios_notx.tsv: $(corr_out) $(gene_info)
	mkdir -p $(dir $@)
	cnvkit.py import-rna -f counts -c $< -g $(gene_info) \
		--no-tx \
		-o $@ --output-dir notx/ $(tcga_rna_counts)

build/tcga_combined_log2ratios_nogctx.tsv: $(corr_out) $(gene_info)
	mkdir -p $(dir $@)
	cnvkit.py import-rna -f counts -c $< -g $(gene_info) \
		--no-gc --no-tx \
		-o $@ --output-dir nogctx/ $(tcga_rna_counts)


# CNVkit smoothing and segmentation ===

# $(asmooth50_cnr): adaptive_smooth.py $(plain_cnr)
#         python $^ -w 50 -d asmooth50

$(firstword $(tsmooth100_cnr)): smooth_tirosh.py $(plain_cnr)
	mkdir -p $(dir $@)
	python $^ -w 100 -d $(dir $@)

$(firstword $(tsmooth200_cnr)): smooth_tirosh.py $(plain_cnr)
	mkdir -p $(dir $@)
	python $^ -w 200 -d $(dir $@)

$(firstword $(nocorrsmooth100_cnr)): smooth_cnr.py $(nocorr_cnr)
	mkdir -p $(dir $@)
	python $^ -w 100 -d $(dir $@)

$(firstword $(wsmooth100_cnr)): smooth_cnr.py $(plain_cnr)
	mkdir -p $(dir $@)
	python $^ -w 100 -d $(dir $@)

$(firstword $(wsmooth200_cnr)): smooth_cnr.py $(plain_cnr)
	mkdir -p $(dir $@)
	python $^ -w 200 -d $(dir $@)

$(firstword $(wsmooth400_cnr)): smooth_cnr.py $(plain_cnr)
	mkdir -p $(dir $@)
	python $^ -w 400 -d $(dir $@)

$(arm_cns) $(cbs_cns) $(flasso_cns) $(haar_cns) $(hmm_cns): %.cns: %.raw.cns
	cnvkit.py call -m none --center median $< -o $@

$(arm_cns:.cns=.raw.cns): arm/%.raw.cns: out/%.cnr
	cnvkit.py segment -m none $< -o $@

$(cbs_cns:.cns=.raw.cns): cbs/%.raw.cns: out/%.cnr
	mkdir -p $(dir $@)
	cnvkit.py segment -m cbs -t .001 -p 2 $< -o $@

$(flasso_cns:.cns=.raw.cns): flasso/%.raw.cns: out/%.cnr
	mkdir -p $(dir $@)
	cnvkit.py segment -m flasso -t .001 $< -o $@

$(haar_cns:.cns=.raw.cns): haar/%.raw.cns: out/%.cnr
	mkdir -p $(dir $@)
	cnvkit.py segment -m haar $< -o $@

$(hmm_cns:.cns=.raw.cns): hmm/%.raw.cns: out/%.cnr
	mkdir -p $(dir $@)
	cnvkit.py segment -m hmm $< -o $@


# Intermediate data tables (used for plotting) ===

build/tcga-acgh-seg-genes.tsv: collate_by_gene.py $(acgh_cns)
	mkdir -p $(dir $@)
	python $^ -g $(gene_info) -o $@ -s ${@:.tsv=.size.tsv}

build/tcga-rna-plain-genes.tsv: collate_by_gene.py $(plain_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-nogc-genes.tsv: collate_by_gene.py $(nogc_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-notx-genes.tsv: collate_by_gene.py $(notx_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-nogctx-genes.tsv: collate_by_gene.py $(nogctx_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-nocorrsmooth100-genes.tsv: collate_by_gene.py $(nocorrsmooth100_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-wsmooth100-genes.tsv: collate_by_gene.py $(wsmooth100_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-wsmooth200-genes.tsv: collate_by_gene.py $(wsmooth200_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-wsmooth400-genes.tsv: collate_by_gene.py $(wsmooth400_cnr)
	python $^ -g $(gene_info) -o $@

# tcga-rna-asmooth50-genes.tsv: collate_by_gene.py $(asmooth50_cnr)
#         python $^ -g $(gene_info) -o $@

build/tcga-rna-tsmooth100-genes.tsv: collate_by_gene.py $(tsmooth100_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-tsmooth200-genes.tsv: collate_by_gene.py $(tsmooth200_cnr)
	python $^ -g $(gene_info) -o $@

build/tcga-rna-arm-genes.tsv: collate_by_gene.py $(arm_cns)
	python $^ -g $(gene_info) -o $@ -s ${@:.tsv=.size.tsv}

build/tcga-rna-flasso-genes.tsv: collate_by_gene.py $(flasso_cns)
	python $^ -g $(gene_info) -o $@ -s ${@:.tsv=.size.tsv}

build/tcga-rna-cbs-genes.tsv: collate_by_gene.py $(cbs_cns)
	python $^ -g $(gene_info) -o $@ -s ${@:.tsv=.size.tsv}

build/tcga-rna-cbsgt20-genes.tsv: collate_by_gene.py $(cbs_cns)
	python $^ -g $(gene_info) -o $@ -s ${@:.tsv=.size.tsv} -w 20

build/tcga-rna-haar-genes.tsv: collate_by_gene.py $(haar_cns)
	python $^ -g $(gene_info) -o $@ -s ${@:.tsv=.size.tsv}

build/tcga-rna-hmm-genes.tsv: collate_by_gene.py $(hmm_cns)
	python $^ -g $(gene_info) -o $@ -s ${@:.tsv=.size.tsv}

build/tcga-rna-hmmgt20-genes.tsv: collate_by_gene.py $(hmm_cns)
	python $^ -g $(gene_info) -o $@ -s ${@:.tsv=.size.tsv} -w 20


# Plots ===
# - 2D scatter of aCGH vs. each RNA estimator (separately), with correlation
# - Violin or box plot of aCGH vs. each RNA estimator (together), with SD/IQR
# - Violin or box plot of each bias correction/not vs. aCGH (together), with SD/IQR

$(plots_2d): plots/rna-vs-acgh-%.pdf: \
	plot2d.py build/tcga-rna-%-genes.tsv build/tcga-acgh-seg-genes.tsv
	python $^ -o $@ -s $(patsubst %.tsv,%.size.tsv,$(lastword $^))

$(plot_resid): \
	plot_residuals.py build/tcga-acgh-seg-genes.tsv \
	$(patsubst %,build/tcga-rna-%-genes.tsv,$(resid_labels))
	python $^ -s build/tcga-acgh-seg-genes.size.tsv -o $@

$(plot_bias): \
	plot_residuals.py build/tcga-acgh-seg-genes.tsv \
	$(patsubst %,build/tcga-rna-%-genes.tsv,$(bias_labels))
	python $^ -o $@

$(plot_bias_facet): \
	plot_residuals.py build/tcga-acgh-seg-genes.tsv \
	$(patsubst %,build/tcga-rna-%-genes.tsv,$(bias_labels))
	python $^ -s build/tcga-acgh-seg-genes.size.tsv -o $@

# $(heatmaps): plots/heatmap-%.pdf: $(foreach $(heat_labels,%/...)
#         cnvkit.py heatmap -d --no-shift-xy $^ -o $@
plots/heatmap-arm.pdf: $(arm_cns)
	mkdir -p $(dir $@)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

plots/heatmap-cbs.pdf: $(cbs_cns)
	mkdir -p $(dir $@)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

plots/heatmap-flasso.pdf: $(flasso_cns)
	mkdir -p $(dir $@)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

plots/heatmap-haar.pdf: $(haar_cns)
	mkdir -p $(dir $@)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

plots/heatmap-hmm.pdf: $(hmm_cns)
	mkdir -p $(dir $@)
	cnvkit.py heatmap -d --no-shift-xy $^ -o $@

$(firstword $(sample_bias_gc)): plot_bias.py $(plain_cnr)
	mkdir -p $(dir $@)
	python $^ -k gc -d $(dir $@)

$(firstword $(sample_bias_tx)): plot_bias.py $(plain_cnr)
	mkdir -p $(dir $@)
	python $^ -k tx_length -d $(dir $@)
