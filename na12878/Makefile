# The NA12878 exome from EdgeBio
# http://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html#exome-with-validation-against-reference-materials
# Coverage obtained with several parameter set permutations for comparison

all: giab.cns


# Intervals

exome_bed := ngv3.hg19.bed
duke_exclude_bed := ../intervals/wgEncodeDukeMapabilityRegionsExcludable.hg19.bed
lumpy_exclude_bed := ../intervals/lumpy-exclude.hg19.bed
hg19_fa := ~/db/ucsc.hg19.fasta
refflat_db := ~/db/annot/refFlat.txt


$(exome_bed): NGv3.bed
	sed 's/.*/chr&/g' $< > $@

ngv3-access.hg19.bed: 
	cnvkit.py access $(hg19_fa) -x $(duke_exclude_bed) -x $(lumpy_exclude_bed) -o $@

ngv3.target-100.bed: $(exome_bed)
	cnvkit.py target $< --annotate $(refflat_db) --split -a 100 -o $@

ngv3.target-267.bed: $(exome_bed)
	cnvkit.py target $< --annotate $(refflat_db) --split -o $@

ngv3.antitarget-2-30kb.bed: ngv3-access.hg19.bed $(exome_bed)
	cnvkit.py antitarget -g $^ -a 30000 -m 2000 -o $@

ngv3.antitarget-6-90kb.bed: ngv3-access.hg19.bed $(exome_bed)
	cnvkit.py antitarget -g $^ -a 90000 -m 6000 -o $@


# Generic female references

reference-giab-thin.cnn: ngv3.target-100.bed ngv3.antitarget-2-30kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(hg19_fa) -o $@

reference-giab.cnn: ngv3.target-267.bed ngv3.antitarget-6-90kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(hg19_fa) -o $@


# Pipeline, starting from coverages

giab-thin.cnr giab.cnr: %.cnr: \
	%.targetcoverage.cnn %.antitargetcoverage.cnn reference-%.cnn
	cnvkit.py fix $^ -o $@

giab-thin.cns giab.cns: %.cns: %.cnr
	cnvkit.py segment $< -o $@


# Benchmark
# giab-svclassify-deletions-2015-05-22.bed
# giab-svclassify-insertions-2015-05-22.bed
