# Test CNVkit with bigger datasets.
#
# Dependency: cnvkit.py must be in $PATH

cnvkit=cnvkit.py
refgenome_ucsc=/db/ucsc.hg19.fasta
refgenome_nochr=/db/hg19-nochr.fa

# ------------------------------------------------------------------------------
# Targeted resequencing samples ("TR")

tr_tcnn=$(wildcard targeted/TR*.targetcoverage.cnn)
tr_cnrs=$(patsubst targeted/%.targetcoverage.cnn,build/%.cnr,$(tr_tcnn))
tr_segs=$(tr_cnrs:.cnr=.cns)

# Pooled reference (normal samples)
tr_ref_cnns=$(wildcard targeted/TR_*_N.*.cnn)


# ------------------------------------------------------------------------------
# Exome samples ("EX")

ex_tcnn=$(wildcard exome/EX*.targetcoverage.cnn)
ex_cnrs=$(patsubst exome/%.targetcoverage.cnn,build/%.cnr,$(ex_tcnn))
ex_segs=$(ex_cnrs:.cnr=.cns)

# Reference: normal male samples
ex_ref_cnns=$(wildcard exome/EX_*_N.*.cnn)


# ------------------------------------------------------------------------------
#  Action!

.PHONY: tr
tr: heatmap-tr.pdf tr-metrics.csv

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

reference-tr.cnn: $(tr_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_ucsc) -y -o $@

reference-exome.cnn: $(ex_ref_cnns)
	$(cnvkit) reference $^ -f $(refgenome_nochr) -y -o $@


# == Build components

$(tr_cnrs): build/%.cnr: targeted/%.targetcoverage.cnn targeted/%.antitargetcoverage.cnn reference-tr.cnn
	$(cnvkit) fix $^ -o $@

$(ex_cnrs): build/%.cnr: exome/%.targetcoverage.cnn exome/%.antitargetcoverage.cnn reference-exome.cnn
	$(cnvkit) fix $^ -o $@

$(tr_segs) $(ex_segs): %.cns: %.cnr
	$(cnvkit) segment $< -o $@


# == Results

heatmap-tr.pdf: $(tr_segs)
	$(cnvkit) heatmap -d -o $@ $(filter %_T.cns,$^)

heatmap-exome.pdf: $(ex_segs)
	$(cnvkit) heatmap -d -o $@ $(filter %_T.cns,$^)


tr-metrics.csv: $(tr_segs)
	$(cnvkit) metrics $(tr_cnrs) -s $^ -o $@

ex-metrics.csv: $(ex_segs)
	$(cnvkit) metrics $(ex_cnrs) -s $^ -o $@


