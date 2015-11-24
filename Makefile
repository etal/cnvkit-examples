# Test CNVkit with bigger datasets.
#
# Dependency: cnvkit.py must be in $PATH

refgenome_ucsc := ~/db/ucsc.hg19.fasta

# ------------------------------------------------------------------------------
# Cell line validation ("CL")

cl_ref_cnns := $(wildcard cell/normal/*_N.*.cnn)
cl_tcnn := cell/CL_seq.targetcoverage.cnn
cl_cnrs := build/CL_seq.cnr build/CL_acgh.cnr
cl_cnrs_flat := build/CL_seq_flat.cnr
cl_segs := $(cl_cnrs:.cnr=.cns) $(cl_cnrs_flat:.cnr=.cns)


# ------------------------------------------------------------------------------
# Targeted resequencing samples ("TR")

tr_ref_cnns := $(wildcard targeted/TR_*_N.*.cnn)
tr_tcnn := $(wildcard targeted/TR*.targetcoverage.cnn)
tr_cnrs := $(patsubst targeted/%.targetcoverage.cnn,build/%.cnr,$(tr_tcnn))
tr_segs := $(tr_cnrs:.cnr=.cns)

tr_thin_ref_cnns := $(wildcard tr-thin/TR_*_N_thin.*.cnn)
tr_thin_tcnn := $(wildcard tr-thin/TR*_thin.targetcoverage.cnn)
tr_thin_cnrs := $(patsubst tr-thin/%.targetcoverage.cnn,build/%.cnr,$(tr_thin_tcnn))
tr_thin_segs := $(tr_thin_cnrs:.cnr=.cns)


# ------------------------------------------------------------------------------
# Exome samples ("EX")

ex_ref_cnns := $(wildcard exome/EX_*_N.*.cnn)
ex_tcnn := $(wildcard exome/EX*.targetcoverage.cnn)
ex_cnrs := $(patsubst exome/%.targetcoverage.cnn,build/%.cnr,$(ex_tcnn))
ex_segs := $(ex_cnrs:.cnr=.cns)


# ------------------------------------------------------------------------------
#  Action!

all: cl tr ex


.PHONY: cl
cl: heatmap-cl.pdf

.PHONY: tr
tr: heatmap-tr.pdf TR_95_T-diagram.pdf TR_95_T-scatter.pdf TR_95_T-CDK4-MDM2-scatter.pdf

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
	cnvkit.py reference $^ -f $(refgenome_ucsc) -y -o $@

reference-tr-thin.cnn: $(tr_thin_ref_cnns)
	cnvkit.py reference $^ -f $(refgenome_ucsc) -y -o $@

reference-tr.cnn: $(tr_ref_cnns)
	cnvkit.py reference $^ -f $(refgenome_ucsc) -y -o $@

reference-exome.cnn: $(ex_ref_cnns)
	cnvkit.py reference $^ -f $(refgenome_ucsc) -y -o $@


# == Build components

build/CL_seq.cnr: build/%.cnr: cell/%.targetcoverage.cnn cell/%.antitargetcoverage.cnn reference-cell.cnn
	cnvkit.py fix $^ -o $@

build/CL_acgh.cnr: cell/CL_acgh.cnr
	cp $< $@

$(cl_cnrs_flat): build/%_flat.cnr: cell/%.targetcoverage.cnn cell/%.antitargetcoverage.cnn intervals/reference-cl-flat.cnn
	cnvkit.py fix $^ -o $@

$(tr_thin_cnrs): build/%.cnr: tr-thin/%.targetcoverage.cnn tr-thin/%.antitargetcoverage.cnn reference-tr-thin.cnn
	cnvkit.py fix $^ -o $@

$(tr_cnrs): build/%.cnr: targeted/%.targetcoverage.cnn targeted/%.antitargetcoverage.cnn reference-tr.cnn
	cnvkit.py fix $^ -o $@

$(ex_cnrs): build/%.cnr: exome/%.targetcoverage.cnn exome/%.antitargetcoverage.cnn reference-exome.cnn
	cnvkit.py fix $^ -o $@

$(cl_segs) $(tr_thin_segs) $(tr_segs) $(ex_segs): %.cns: %.cnr
	cnvkit.py segment $< -o $@


# == Results

heatmap-cl.pdf: $(cl_segs)
	cnvkit.py heatmap -d -o $@ $^

heatmap-tr-thin.pdf: $(tr_thin_segs)
	cnvkit.py heatmap -d -o $@ $(filter %_T_thin.cns,$^)

heatmap-tr.pdf: $(tr_segs)
	cnvkit.py heatmap -d -o $@ $(filter %_T.cns,$^)

heatmap-tr-nod.pdf: $(tr_segs)
	cnvkit.py heatmap -o $@ $(filter %_T.cns,$^)

heatmap-exome.pdf: $(ex_segs)
	cnvkit.py heatmap -d -o $@ $(filter %_T.cns,$^)


cl-metrics.csv: $(cl_segs)
	cnvkit.py metrics $(cl_segs:.cns=.cnr) -s $^ -o $@

tr-thin-metrics.csv: $(tr_thin_segs)
	cnvkit.py metrics $(tr_thin_cnrs) -s $^ -o $@

tr-metrics.csv: $(tr_segs)
	cnvkit.py metrics $(tr_cnrs) -s $^ -o $@

ex-metrics.csv: $(ex_segs)
	cnvkit.py metrics $(ex_cnrs) -s $^ -o $@

# Example figures
TR_95_T-diagram.pdf: build/TR_95_T.cns build/TR_95_T.cnr
	cnvkit.py diagram -s $^ -y -o $@

TR_95_T-scatter.pdf: build/TR_95_T.cns build/TR_95_T.cnr
	cnvkit.py scatter -s $^ -o $@

TR_95_T-CDK4-MDM2-scatter.pdf: build/TR_95_T.cns build/TR_95_T.cnr
	cnvkit.py scatter -s $^ -o $@ -c chr12:50000000-80000000 -g CDK4,MDM2

