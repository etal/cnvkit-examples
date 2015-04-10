# Flat reference

refgenome_ucsc=~/db/ucsc.hg19.fasta
tr_tgt_bed=../../intervals/tr.target-267.bed
tr_anti_bed=../../intervals/tr.antitarget-10-150kb.bed
ex_tgt_bed=../../intervals/ex.target-267.bed
ex_anti_bed=../../intervals/ex.antitarget-6-90kb.bed

tr_tcnn=$(wildcard ../../targeted/TR*.targetcoverage.cnn)
tr_cnrs_flat=$(patsubst ../../targeted/%.targetcoverage.cnn,%_flat.cnr,$(tr_tcnn))
tr_segs_flat=$(tr_cnrs_flat:.cnr=.cns)

ex_tcnn=$(wildcard ../../exome/EX*.targetcoverage.cnn)
ex_cnrs_flat=$(patsubst ../../exome/%.targetcoverage.cnn,%_flat.cnr,$(ex_tcnn))
ex_segs_flat=$(ex_cnrs_flat:.cnr=.cns)


all: tr ex

.PHONY: tr
tr: heatmap-tr-flat.pdf metrics-tr-flat.csv

.PHONY: ex
ex: $(ex_segs) metrics-ex-flat.csv heatmap-ex-flat.pdf

.PHONY: clean
clean:
	# Targeted
	rm -vf TR* heatmap-tr*.pdf
	# Exome
	rm -vf EX* heatmap-ex-flat.pdf


reference-tr-flat.cnn: $(tr_tgt_bed) $(tr_anti_bed)
	cnvkit.py reference -t $(tr_tgt_bed) -a $(tr_anti_bed) -f $(refgenome_ucsc) -y -o $@

reference-ex-flat.cnn: $(ex_tgt_bed) $(ex_anti_bed)
	cnvkit.py reference -t $(ex_tgt_bed) -a $(ex_anti_bed) -f $(refgenome_ucsc) -y -o $@


$(tr_cnrs_flat): %_flat.cnr: ../../targeted/%.targetcoverage.cnn ../../targeted/%.antitargetcoverage.cnn reference-tr-flat.cnn
	cnvkit.py fix $^ -o $@

$(ex_cnrs_flat): %_flat.cnr: ../../exome/%.targetcoverage.cnn ../../exome/%.antitargetcoverage.cnn reference-ex-flat.cnn
	cnvkit.py fix $^ -o $@

$(tr_segs_flat) $(ex_segs_flat): %.cns: %.cnr
	cnvkit.py segment $< -o $@


metrics-tr-flat.csv: $(tr_segs_flat)
	cnvkit.py metrics $(tr_cnrs_flat) -s $^ -o $@

metrics-ex-flat.csv: $(ex_segs_flat)
	cnvkit.py metrics $(ex_cnrs_flat) -s $^ -o $@


heatmap-tr-flat.pdf: $(tr_segs_flat)
	cnvkit.py heatmap -d -o $@ $(filter %_T_flat.cns,$^)

heatmap-ex-flat.pdf: $(ex_segs_flat)
	cnvkit.py heatmap -d -o $@ $(filter %_T_flat.cns,$^)

