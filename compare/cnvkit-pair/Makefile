# Paired references

# refgenome_ucsc=~/db/ucsc.hg19.fasta
# tr_tgt_bed=../../intervals/tr.target-267.bed
# tr_anti_bed=../../intervals/tr.antitarget-10-150kb.bed
# ex_tgt_bed=../../intervals/ex.target-267.bed
# ex_anti_bed=../../intervals/ex.antitarget-6-90kb.bed

tr_tcnn=$(wildcard ../../targeted/TR*.targetcoverage.cnn)
tr_cnrs=$(patsubst ../../targeted/%.targetcoverage.cnn,%_pair.cnr,$(filter %_T.targetcoverage.cnn,$(tr_tcnn)))
tr_segs=$(tr_cnrs:.cnr=.cns)
tr_orig_ref=../../reference-tr.cnn
tr_references=$(patsubst %_T_pair.cnr,reference-%.cnn,$(tr_cnrs))

ex_tcnn=$(wildcard ../../exome/EX*.targetcoverage.cnn)
ex_cnrs=$(patsubst ../../exome/%.targetcoverage.cnn,%_pair.cnr,$(filter %_T.targetcoverage.cnn,$(ex_tcnn)))
ex_segs=$(ex_cnrs:.cnr=.cns)
ex_orig_ref=../../reference-exome.cnn
ex_references=$(patsubst %_T_pair.cnr,reference-%.cnn,$(ex_cnrs))

cl_cnrs=CL_pair.cnr
cl_segs=$(cl_cnrs:.cnr=.cns)
cl_orig_ref=../../reference-cell.cnn
cl_references=reference-CL.cnn


all: tr ex cl

.PHONY: tr
tr: heatmap-tr-paired.pdf

.PHONY: ex
ex: $(ex_segs) heatmap-ex-paired.pdf

.PHONY: cl
cl: $(cl_segs) $(cl_segs:.cns=.bed)

.PHONY: metrics
metrics:  metrics-tr-paired.csv metrics-ex-paired.csv metrics-cl-paired.csv

.PHONY: clean
clean:
	# Targeted
	rm -vf TR* heatmap-tr*.pdf
	# Exome
	rm -vf EX* heatmap-ex-paired.pdf


$(tr_references): reference-%.cnn: ../../targeted/%_N.targetcoverage.cnn ../../targeted/%_N.antitargetcoverage.cnn
	cnvkit.py reference $^ -y -o $@.tmp.cnn
	python ../add_gc_rmask.py $(tr_orig_ref) $@.tmp.cnn -o $@
	rm $@.tmp.cnn

$(ex_references): reference-%.cnn: ../../exome/%_N.targetcoverage.cnn ../../exome/%_N.antitargetcoverage.cnn
	cnvkit.py reference $^ -y -o $@.tmp.cnn
	python ../add_gc_rmask.py $(ex_orig_ref) $@.tmp.cnn -o $@
	rm $@.tmp.cnn

$(cl_references): ../../cell/normal/MB_1454_N.targetcoverage.cnn ../../cell/normal/MB_1454_N.antitargetcoverage.cnn
	cnvkit.py reference $^ -y -o $@.tmp.cnn
	python ../add_gc_rmask.py $(cl_orig_ref) $@.tmp.cnn -o $@
	rm $@.tmp.cnn


$(tr_cnrs): %_T_pair.cnr: ../../targeted/%_T.targetcoverage.cnn ../../targeted/%_T.antitargetcoverage.cnn reference-%.cnn
	cnvkit.py fix $^ -o $@

$(ex_cnrs): %_T_pair.cnr: ../../exome/%_T.targetcoverage.cnn ../../exome/%_T.antitargetcoverage.cnn reference-%.cnn
	cnvkit.py fix $^ -o $@

$(cl_cnrs): ../../cell/CL_seq.targetcoverage.cnn ../../cell/CL_seq.antitargetcoverage.cnn reference-CL.cnn
	cnvkit.py fix $^ -o $@


$(cl_segs) $(tr_segs) $(ex_segs): %.cns: %.cnr
	cnvkit.py segment --drop-low $< -o $@

$(cl_segs:.cns=.bed): %.bed: %.cns
	cnvkit.py export bed $< --ploidy 6 --show ploidy -g f -y -o $@


metrics-tr-paired.csv: $(tr_segs)
	cnvkit.py metrics $(tr_cnrs) -s $^ -o $@

metrics-ex-paired.csv: $(ex_segs)
	cnvkit.py metrics $(ex_cnrs) -s $^ -o $@

metrics-cl-paired.csv: $(cl_segs)
	cnvkit.py metrics $(cl_cnrs) -s $^ -o $@


heatmap-tr-paired.pdf: $(tr_segs)
	cnvkit.py heatmap -d -o $@ $^

heatmap-ex-paired.pdf: $(ex_segs)
	cnvkit.py heatmap -d -o $@ $^

