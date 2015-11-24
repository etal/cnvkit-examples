# Build the targets, antitargets and flat references for CL, TR and EX datasets

cell_intervals := cell-baits.bed
target_intervals := targeted-baits.interval_list
exome_intervals := exome-baits.interval_list
access_intervals := access-5k-mappable.hg19.bed
genome_fasta := ~/db/ucsc.hg19.fasta
exclude_bed := wgEncodeDacMapabilityConsensusExcludable.hg19.bed wgEncodeDukeMapabilityRegionsExcludable.hg19.bed


all: reference-cl-flat.cnn reference-tr-flat.cnn reference-ex-flat.cnn

clean:
	rm -vf reference-*.cnn \
		cl.target-*.bed cl.antitarget-*.bed \
		tr.target-*.bed tr.antitarget-*.bed \
		ex.target-*.bed ex.antitarget-*.bed


$(access_intervals): $(exclude_bed)
	cnvkit.py access $(genome_fasta) -x $< -x $(lastword $^) -s 5000 -o $@


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


# EX

ex.target-267.bed: $(exome_intervals)
	cnvkit.py target $< --split -o $@

ex.antitarget-6-90kb.bed: $(access_intervals) $(exome_intervals)
	cnvkit.py antitarget -a 90000 -m 6000 -g $^ -o $@

reference-ex-flat.cnn: ex.target-267.bed ex.antitarget-6-90kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@
