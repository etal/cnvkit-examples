# Build the targets, antitargets and flat references for TR and EX datasets

target_intervals=targeted-baits.interval_list
exome_intervals=exome-baits.interval_list
access_intervals=access-10kb.hg19.bed
genome_fasta=~/db/ucsc.hg19.fasta

all: reference-tr-flat-thin.cnn reference-tr-flat-wide.cnn reference-ex-flat-thin.cnn reference-ex-flat-wide.cnn

clean:
	rm -vf tr.antitarget-*kb.bed tr.target-*.bed ex.antitarget-*kb.bed ex.target-*.bed reference-{ex,tr}-flat-{thin,wide}.cnn


# TR

tr.target-100.bed: $(target_intervals)
	cnvkit.py target $< --split -a 100 -o  $@

tr.target-267.bed: $(target_intervals)
	cnvkit.py target $< --split -o $@

tr.antitarget-5-75kb.bed: $(target_intervals)
	cnvkit.py antitarget $< -g $(access_intervals) -a 75000 -m 5000 -o $@

tr.antitarget-10-150kb.bed: $(target_intervals)
	cnvkit.py antitarget $< -g $(access_intervals) -a 150000 -m 10000 -o $@

reference-tr-flat-thin.cnn: tr.target-100.bed tr.antitarget-5-75kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@

reference-tr-flat-wide.cnn: tr.target-267.bed tr.antitarget-10-150kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@


# EX

ex.target-100.bed: $(exome_intervals)
	cnvkit.py target $< --split -a 100 -o $@

ex.target-267.bed: $(exome_intervals) 
	cnvkit.py target $< --split -o $@

ex.antitarget-3-45kb.bed: $(exome_intervals)
	cnvkit.py antitarget $< -a 45000 -m 3000 -g $(access_intervals) -o $@

ex.antitarget-6-90kb.bed: $(exome_intervals)
	cnvkit.py antitarget $< -a 90000 -m 6000 -g $(access_intervals) -o $@

reference-ex-flat-thin.cnn: ex.target-100.bed ex.antitarget-3-45kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@

reference-ex-flat-wide.cnn: ex.target-267.bed ex.antitarget-6-90kb.bed
	cnvkit.py reference -t $< -a $(lastword $^) -f $(genome_fasta) -y -o $@
