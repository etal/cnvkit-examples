# Test CNVkit with bigger datasets.
#
# Dependency: cnvkit.py must be in $PATH

cnvkit=cnvkit.py
refgenome_ucsc=~/db/ucsc.hg19.fasta

# ------------------------------------------------------------------------------
# Targeted resequencing samples ("TR")

tr_thin_ref_cnns=$(wildcard tr-thin/TR_*_N_thin.*.cnn)
tr_thin_tcnn=$(wildcard tr-thin/TR*_thin.targetcoverage.cnn)
tr_thin_cnrs=$(patsubst tr-thin/%.targetcoverage.cnn,build/%.cnr,$(tr_thin_tcnn))
tr_thin_segs=$(tr_thin_cnrs:.cnr=.cns)
# Flat reference
tr_thin_cnrs_flat=$(patsubst tr-thin/%.targetcoverage.cnn,build/%_flat.cnr,$(tr_thin_tcnn))
tr_thin_segs_flat=$(tr_thin_cnrs_flat:.cnr=.cns)

tr_wide_ref_cnns=$(wildcard tr-wide/TR_*_N_wide.*.cnn)
tr_wide_tcnn=$(wildcard tr-wide/TR*_wide.targetcoverage.cnn)
tr_wide_cnrs=$(patsubst tr-wide/%.targetcoverage.cnn,build/%.cnr,$(tr_wide_tcnn))
tr_wide_segs=$(tr_wide_cnrs:.cnr=.cns)
# Flat reference
tr_wide_cnrs_flat=$(patsubst tr-wide/%.targetcoverage.cnn,build/%_flat.cnr,$(tr_wide_tcnn))
tr_wide_segs_flat=$(tr_wide_cnrs_flat:.cnr=.cns)


# ------------------------------------------------------------------------------
# Exome samples ("EX")

ex_ref_cnns=$(wildcard exome/EX_*_N.*.cnn)
ex_tcnn=$(wildcard exome/EX*.targetcoverage.cnn)
ex_cnrs=$(patsubst exome/%.targetcoverage.cnn,build/%.cnr,$(ex_tcnn))
ex_segs=$(ex_cnrs:.cnr=.cns)
# Flat reference
# ex_cnrs_flat=$(patsubst exome/%.targetcoverage.cnn,build/%_flat.cnr,$(ex_tcnn))
# ex_segs_flat=$(ex_cnrs:.cnr=.cns)


# ------------------------------------------------------------------------------
#  Action!

.PHONY: tr
tr: heatmap-tr-thin.pdf heatmap-tr-wide.pdf tr-thin-metrics.csv tr-wide-metrics.csv \
	heatmap-tr-thin-flat.pdf heatmap-tr-wide-flat.pdf tr-thin-metrics-flat.csv tr-wide-metrics-flat.csv

.PHONY: ex
ex: $(ex_segs) ex-metrics.csv heatmap-exome.pdf


all: tr ex


.PHONY: clean
clean:
	# Targeted
	rm -vf build/TR_* heatmap-tr.pdf
	# Exome
	rm -vf build/EX_* heatmap-exome.pdf

# ------------------------------------------------------------------------------
# Standard workflow

# == Build pooled references from normal samples

reference-tr-thin.cnn: $(tr_thin_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_ucsc) -y -o $@

reference-tr-wide.cnn: $(tr_wide_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_ucsc) -y -o $@

reference-exome.cnn: $(ex_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_ucsc) -y -o $@


# == Build components

$(tr_thin_cnrs): build/%.cnr: tr-thin/%.targetcoverage.cnn tr-thin/%.antitargetcoverage.cnn reference-tr-thin.cnn
	$(cnvkit) fix $^ -o $@

$(tr_wide_cnrs): build/%.cnr: tr-wide/%.targetcoverage.cnn tr-wide/%.antitargetcoverage.cnn reference-tr-wide.cnn
	$(cnvkit) fix $^ -o $@

$(tr_thin_cnrs_flat): build/%_flat.cnr: tr-thin/%.targetcoverage.cnn tr-thin/%.antitargetcoverage.cnn intervals/reference-tr-flat-thin.cnn
	$(cnvkit) fix $^ -o $@

$(tr_wide_cnrs_flat): build/%_flat.cnr: tr-wide/%.targetcoverage.cnn tr-wide/%.antitargetcoverage.cnn intervals/reference-tr-flat-wide.cnn
	$(cnvkit) fix $^ -o $@

$(ex_cnrs): build/%.cnr: exome/%.targetcoverage.cnn exome/%.antitargetcoverage.cnn reference-exome.cnn
	$(cnvkit) fix $^ -o $@

$(ex_cnrs_flat): build/%_flat.cnr: exome/%.targetcoverage.cnn exome/%.antitargetcoverage.cnn intervals/reference-ex-flat-wide.cnn
	$(cnvkit) fix $^ -o $@

$(tr_thin_segs) $(tr_wide_segs) $(ex_segs) $(tr_thin_segs_flat) $(tr_wide_segs_flat) $(ex_segs_flat): %.cns: %.cnr
	$(cnvkit) segment $< -o $@


# == Results

heatmap-tr-thin.pdf: $(tr_thin_segs)
	$(cnvkit) heatmap -d -o $@ $(filter %_T_thin.cns,$^)

heatmap-tr-wide.pdf: $(tr_wide_segs)
	$(cnvkit) heatmap -d -o $@ $(filter %_T_wide.cns,$^)

heatmap-exome.pdf: $(ex_segs)
	$(cnvkit) heatmap -d -o $@ $(filter %_T.cns,$^)


tr-thin-metrics.csv: $(tr_thin_segs)
	$(cnvkit) metrics $(tr_thin_cnrs) -s $^ -o $@

tr-wide-metrics.csv: $(tr_wide_segs)
	$(cnvkit) metrics $(tr_wide_cnrs) -s $^ -o $@

ex-metrics.csv: $(ex_segs)
	$(cnvkit) metrics $(ex_cnrs) -s $^ -o $@


heatmap-tr-thin-flat.pdf: $(tr_thin_segs_flat)
	$(cnvkit) heatmap -d -o $@ $(filter %_T_thin_flat.cns,$^)

heatmap-tr-wide-flat.pdf: $(tr_wide_segs_flat)
	$(cnvkit) heatmap -d -o $@ $(filter %_T_wide_flat.cns,$^)

heatmap-exome-flat.pdf: $(ex_segs_flat)
	$(cnvkit) heatmap -d -o $@ $(filter %_T_flat.cns,$^)


tr-thin-metrics-flat.csv: $(tr_thin_segs_flat)
	$(cnvkit) metrics $(tr_thin_cnrs_flat) -s $^ -o $@

tr-wide-metrics-flat.csv: $(tr_wide_segs_flat)
	$(cnvkit) metrics $(tr_wide_cnrs_flat) -s $^ -o $@

ex-metrics-flat.csv: $(ex_segs_flat)
	$(cnvkit) metrics $(ex_cnrs) -s $^ -o $@
