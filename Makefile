# Test CNVkit with bigger datasets.
#
# Dependency: cnvkit.py must be in $PATH

cnvkit=cnvkit.py
refgenome_ucsc=~/db/ucsc.hg19.fasta

# ------------------------------------------------------------------------------
# Cell line validation ("CL")

cl_ref_cnns=$(wildcard cell/normal/*_N.*.cnn)
cl_tcnn=cell/CL_seq.targetcoverage.cnn
cl_cnrs=build/CL_seq.cnr build/CL_acgh.cnr
cl_cnrs_flat=build/CL_seq_flat.cnr
cl_segs=$(cl_cnrs:.cnr=.cns) $(cl_cnrs_flat:.cnr=.cns)


# ------------------------------------------------------------------------------
# Targeted resequencing samples ("TR")

tr_ref_cnns=$(wildcard targeted/TR_*_N.*.cnn)
tr_tcnn=$(wildcard targeted/TR*.targetcoverage.cnn)
tr_cnrs=$(patsubst targeted/%.targetcoverage.cnn,build/%.cnr,$(tr_tcnn))
tr_segs=$(tr_cnrs:.cnr=.cns)

tr_thin_ref_cnns=$(wildcard tr-thin/TR_*_N_thin.*.cnn)
tr_thin_tcnn=$(wildcard tr-thin/TR*_thin.targetcoverage.cnn)
tr_thin_cnrs=$(patsubst tr-thin/%.targetcoverage.cnn,build/%.cnr,$(tr_thin_tcnn))
tr_thin_segs=$(tr_thin_cnrs:.cnr=.cns)
# Flat reference
# tr_thin_cnrs_flat=$(patsubst tr-thin/%.targetcoverage.cnn,build/%_flat.cnr,$(tr_thin_tcnn))
# tr_thin_segs_flat=$(tr_thin_cnrs_flat:.cnr=.cns)


# ------------------------------------------------------------------------------
# Exome samples ("EX")

ex_ref_cnns=$(wildcard exome/EX_*_N.*.cnn)
ex_tcnn=$(wildcard exome/EX*.targetcoverage.cnn)
ex_cnrs=$(patsubst exome/%.targetcoverage.cnn,build/%.cnr,$(ex_tcnn))
ex_segs=$(ex_cnrs:.cnr=.cns)


# ------------------------------------------------------------------------------
#  Action!

all: cl tr ex


.PHONY: cl
cl: heatmap-cl.pdf

.PHONY: tr
tr: heatmap-tr.pdf

.PHONY: ex
ex: $(ex_segs) heatmap-exome.pdf

.PHONY: metrics
metrics: ex-metrics.csv tr-metrics.csv cl-metrics.csv


.PHONY: clean
clean:
	# Targeted
	rm -vf build/TR* heatmap-tr*.pdf
	# Exome
	rm -vf build/EX* heatmap-exome.pdf
	# Cell
	rm -vf build/CL* heatmap-cl.pdf


# ------------------------------------------------------------------------------
# Standard workflow

# == Build pooled references from normal samples

reference-cell.cnn: $(cl_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_ucsc) -y -o $@

reference-tr-thin.cnn: $(tr_thin_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_ucsc) -y -o $@

reference-tr.cnn: $(tr_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_ucsc) -y -o $@

reference-exome.cnn: $(ex_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_ucsc) -y -o $@


# == Build components

build/CL_seq.cnr: build/%.cnr: cell/%.targetcoverage.cnn cell/%.antitargetcoverage.cnn reference-cell.cnn
	$(cnvkit) fix $^ -o $@

build/CL_acgh.cnr: cell/CL_acgh.cnr
	cp $< $@

$(cl_cnrs_flat): build/%_flat.cnr: cell/%.targetcoverage.cnn cell/%.antitargetcoverage.cnn intervals/reference-cl-flat.cnn
	$(cnvkit) fix $^ -o $@

$(tr_thin_cnrs): build/%.cnr: tr-thin/%.targetcoverage.cnn tr-thin/%.antitargetcoverage.cnn reference-tr-thin.cnn
	$(cnvkit) fix $^ -o $@

$(tr_cnrs): build/%.cnr: targeted/%.targetcoverage.cnn targeted/%.antitargetcoverage.cnn reference-tr.cnn
	$(cnvkit) fix $^ -o $@

# $(tr_thin_cnrs_flat): build/%_flat.cnr: tr-thin/%.targetcoverage.cnn tr-thin/%.antitargetcoverage.cnn intervals/reference-tr-flat-thin.cnn
#         $(cnvkit) fix $^ -o $@

$(ex_cnrs): build/%.cnr: exome/%.targetcoverage.cnn exome/%.antitargetcoverage.cnn reference-exome.cnn
	$(cnvkit) fix $^ -o $@

$(cl_segs) $(tr_thin_segs) $(tr_segs) $(ex_segs): %.cns: %.cnr
	$(cnvkit) segment $< -o $@


# == Results

heatmap-cl.pdf: $(cl_segs)
	$(cnvkit) heatmap -d -o $@ $^

heatmap-tr-thin.pdf: $(tr_thin_segs)
	$(cnvkit) heatmap -d -o $@ $(filter %_T_thin.cns,$^)

heatmap-tr.pdf: $(tr_segs)
	$(cnvkit) heatmap -d -o $@ $(filter %_T.cns,$^)

heatmap-exome.pdf: $(ex_segs)
	$(cnvkit) heatmap -d -o $@ $(filter %_T.cns,$^)


cl-metrics.csv: $(cl_segs)
	$(cnvkit) metrics $(cl_segs:.cns=.cnr) -s $^ -o $@

tr-thin-metrics.csv: $(tr_thin_segs)
	$(cnvkit) metrics $(tr_thin_cnrs) -s $^ -o $@

tr-metrics.csv: $(tr_segs)
	$(cnvkit) metrics $(tr_cnrs) -s $^ -o $@

ex-metrics.csv: $(ex_segs)
	$(cnvkit) metrics $(ex_cnrs) -s $^ -o $@


# heatmap-tr-thin-flat.pdf: $(tr_thin_segs_flat)
#         $(cnvkit) heatmap -d -o $@ $(filter %_T_thin_flat.cns,$^)

# tr-thin-metrics-flat.csv: $(tr_thin_segs_flat)
#         $(cnvkit) metrics $(tr_thin_cnrs_flat) -s $^ -o $@

