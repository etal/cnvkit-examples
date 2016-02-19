# Flat reference

refgenome_ucsc=~/db/ucsc.hg19.fasta
ref_tr=../../intervals/reference-tr-flat.cnn
ref_ex=../../intervals/reference-ex-flat.cnn
ref_cl=../../intervals/reference-cl-flat.cnn

tr_tcnn=$(wildcard ../../targeted/TR*.targetcoverage.cnn)
tr_cnrs_flat=$(patsubst ../../targeted/%.targetcoverage.cnn,%_flat.cnr,$(filter %_T.targetcoverage.cnn,$(tr_tcnn)))
tr_segs_flat=$(tr_cnrs_flat:.cnr=.cns)

ex_tcnn=$(wildcard ../../exome/EX*.targetcoverage.cnn)
ex_cnrs_flat=$(patsubst ../../exome/%.targetcoverage.cnn,%_flat.cnr,$(filter %_T.targetcoverage.cnn,$(ex_tcnn)))
ex_segs_flat=$(ex_cnrs_flat:.cnr=.cns)


all: tr ex cl

.PHONY: tr
tr: heatmap-tr-flat.pdf

.PHONY: ex
ex: heatmap-ex-flat.pdf

.PHONY: cl
cl: CL_flat.cns

.PHONY: metrics
metrics: metrics-tr-flat.csv metrics-ex-flat.csv metrics-cl-flat.csv

.PHONY: clean
clean:
	# Targeted
	rm -vf TR* heatmap-tr*.pdf
	# Exome
	rm -vf EX* heatmap-ex-flat.pdf


$(tr_cnrs_flat): %_flat.cnr: ../../targeted/%.targetcoverage.cnn ../../targeted/%.antitargetcoverage.cnn $(ref_tr)
	cnvkit.py fix $^ -o $@

$(ex_cnrs_flat): %_flat.cnr: ../../exome/%.targetcoverage.cnn ../../exome/%.antitargetcoverage.cnn $(ref_ex)
	cnvkit.py fix $^ -o $@

CL_flat.cnr: ../../cell/CL_seq.targetcoverage.cnn ../../cell/CL_seq.antitargetcoverage.cnn $(ref_cl)
	cnvkit.py fix $^ -o $@

CL_flat.cns $(tr_segs_flat) $(ex_segs_flat): %.cns: %.cnr
	cnvkit.py segment --drop-low $< -o $@


metrics-tr-flat.csv: $(tr_segs_flat)
	cnvkit.py metrics $(tr_cnrs_flat) -s $^ -o $@

metrics-ex-flat.csv: $(ex_segs_flat)
	cnvkit.py metrics $(ex_cnrs_flat) -s $^ -o $@

metrics-cl-flat.csv: CL_flat.cns
	cnvkit.py metrics CL_flat.cnr -s $^ -o $@


heatmap-tr-flat.pdf: $(tr_segs_flat)
	cnvkit.py heatmap -d -o $@ $^

heatmap-ex-flat.pdf: $(ex_segs_flat)
	cnvkit.py heatmap -d -o $@ $^

