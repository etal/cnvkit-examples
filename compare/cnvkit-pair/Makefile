# Paired references

refgenome_ucsc=~/db/ucsc.hg19.fasta
# tr_tgt_bed=../../intervals/tr.target-267.bed
# tr_anti_bed=../../intervals/tr.antitarget-10-150kb.bed
# ex_tgt_bed=../../intervals/ex.target-267.bed
# ex_anti_bed=../../intervals/ex.antitarget-6-90kb.bed

tr_tcnn=$(wildcard ../../targeted/TR*.targetcoverage.cnn)
tr_cnrs=$(patsubst ../../targeted/%.targetcoverage.cnn,%_pair.cnr,$(filter %_T.targetcoverage.cnn,$(tr_tcnn)))
tr_segs=$(tr_cnrs:.cnr=.cns)
tr_references=$(patsubst %_T_pair.cnr,reference-%.cnn,$(tr_cnrs))

ex_tcnn=$(wildcard ../../exome/EX*.targetcoverage.cnn)
ex_cnrs=$(patsubst ../../exome/%.targetcoverage.cnn,%_pair.cnr,$(filter %_T.targetcoverage.cnn,$(ex_tcnn)))
ex_segs=$(ex_cnrs:.cnr=.cns)
ex_references=$(patsubst %_T_pair.cnr,reference-%.cnn,$(ex_cnrs))

all: tr ex

.PHONY: tr
tr: heatmap-tr-paired.pdf

.PHONY: ex
ex: $(ex_segs) heatmap-ex-paired.pdf

.PHONY: metrics
metrics:  metrics-tr-paired.csv metrics-ex-paired.csv

.PHONY: clean
clean:
	# Targeted
	rm -vf TR* heatmap-tr*.pdf
	# Exome
	rm -vf EX* heatmap-ex-paired.pdf


$(tr_references): reference-%.cnn: ../../targeted/%_N.targetcoverage.cnn ../../targeted/%_N.antitargetcoverage.cnn
	cnvkit.py reference $^ -f $(refgenome_ucsc) -y -o $@

$(ex_references): reference-%.cnn: ../../exome/%_N.targetcoverage.cnn ../../exome/%_N.antitargetcoverage.cnn
	cnvkit.py reference $^ -f $(refgenome_ucsc) -y -o $@


$(tr_cnrs): %_T_pair.cnr: ../../targeted/%_T.targetcoverage.cnn ../../targeted/%_T.antitargetcoverage.cnn reference-%.cnn
	cnvkit.py fix $^ -o $@

$(ex_cnrs): %_T_pair.cnr: ../../exome/%_T.targetcoverage.cnn ../../exome/%_T.antitargetcoverage.cnn reference-%.cnn
	cnvkit.py fix $^ -o $@


$(tr_segs) $(ex_segs): %.cns: %.cnr
	cnvkit.py segment $< -o $@


metrics-tr-paired.csv: $(tr_segs)
	cnvkit.py metrics $(tr_cnrs) -s $^ -o $@

metrics-ex-paired.csv: $(ex_segs)
	cnvkit.py metrics $(ex_cnrs) -s $^ -o $@


heatmap-tr-paired.pdf: $(tr_segs)
	cnvkit.py heatmap -d -o $@ $^

heatmap-ex-paired.pdf: $(ex_segs)
	cnvkit.py heatmap -d -o $@ $^

