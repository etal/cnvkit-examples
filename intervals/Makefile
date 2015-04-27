# Build the targets, antitargets and flat references for CL, TR and EX datasets

cell_intervals=cell-baits.bed
target_intervals=targeted-baits.interval_list
exome_intervals=exome-baits.interval_list
access_intervals=access-100-mappable.hg19.bed
genome_fasta=~/db/ucsc.hg19.fasta
exclude_bed=wgEncodeDacMapabilityConsensusExcludable.hg19.bed
# wgEncodeDukeMapabilityRegionsExcludable.hg19.bed


# reference-ex-flat-thin.cnn
all: reference-cl-flat.cnn reference-tr-flat.cnn reference-ex-flat.cnn \
	reference-tr-flat-thin.cnn

clean:
	rm -vf reference-*.cnn \
		cl.target-*.bed cl.antitarget-*.bed \
		tr.target-*.bed tr.antitarget-*.bed \
		ex.target-*.bed ex.antitarget-*.bed


$(access_intervals): $(exclude_bed)
	genome2access.py $(genome_fasta) -x $^ -s 100 -o $@


# Cell

cl.target-267.bed: $(cell_intervals)
	cnvkit.py target $< --split -o $@

cl.antitarget-10-150kb.bed: $(access_intervals) $(cell_intervals)
	cnvkit.py antitarget -a 150000 -m 10000 -g $^ -o $@


# TR

tr.target-267.bed: $(target_intervals)
	cnvkit.py target $< --split -o $@

tr.antitarget-10-150kb.bed: $(access_intervals) $(target_intervals)
	cnvkit.py antitarget -a 150000 -m 10000 -g $^ -o $@

reference-cl-flat.cnn reference-tr-flat.cnn: \
	reference-%-flat.cnn: %.target-267.bed %.antitarget-10-150kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@


tr.target-100.bed: $(target_intervals)
	cnvkit.py target $< --split -a 100 -o  $@

tr.antitarget-5-75kb.bed: $(access_intervals) $(target_intervals)
	cnvkit.py antitarget -a 75000 -m 5000 -g $^ -o $@

reference-tr-flat-thin.cnn: tr.target-100.bed tr.antitarget-5-75kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@


# EX

ex.target-267.bed: $(exome_intervals)
	cnvkit.py target $< --split -o $@

ex.antitarget-6-90kb.bed: $(access_intervals) $(exome_intervals)
	cnvkit.py antitarget -a 90000 -m 6000 -g $^ -o $@

reference-ex-flat.cnn: ex.target-267.bed ex.antitarget-6-90kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@


ex.target-100.bed: $(exome_intervals)
	cnvkit.py target $< --split -a 100 -o $@

ex.antitarget-3-45kb.bed: $(access_intervals) $(exome_intervals)
	cnvkit.py antitarget -a 45000 -m 3000 -g $^ -o $@

reference-ex-flat-thin.cnn: ex.target-100.bed ex.antitarget-3-45kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@
