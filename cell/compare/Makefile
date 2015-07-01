# CNV profiles for the C0902 cell line (NimV1 platform)
#
# Genes of interest:
# 	
# NTRK1
# ALK
# BRAF
# MET
# ROS1
# RET
# CDKN2A (not FISH'd, but focally deleted)

baits=../../intervals/cell-baits.bed
cnvbuild=../../build
clseqpair=$(cnvbuild)/CL_seq.cnr $(cnvbuild)/CL_seq.cns
clacghpair=$(cnvbuild)/CL_acgh.cnr $(cnvbuild)/CL_acgh.cns


all: heatmap.pdf scatters.pdf

clean:
	rm -vf *.pdf

heatmap.pdf: $(cnvbuild)/CL_acgh.cns $(cnvbuild)/CL_seq.cns $(cnvbuild)/CL_seq_flat.cns
	cnvkit.py heatmap $^ -d -o $@

scatters.pdf: cell-correlation.pdf cell-correlation-flat.pdf cell-paired-genome.pdf \
	CL_acgh-scatter.pdf CL_seq-scatter.pdf CL_seq_flat-scatter.pdf \
	focal-NTRK1.pdf focal-ALK.pdf focal-ROS1.pdf focal-MET-BRAF.pdf focal-CDKN2A.pdf focal-RET.pdf
	pdfunite $^ $@

# Whole-genome scatter plots

cell-paired-genome.pdf: plot_genome_pair.py $(clseqpair) $(clacghpair)
	python $^ -o $@

CL_acgh-scatter.pdf CL_seq-scatter.pdf CL_seq_flat-scatter.pdf: %-scatter.pdf: $(cnvbuild)/%.cns $(cnvbuild)/%.cnr
	cnvkit.py scatter -s $^ -o $@

# Focal scatter plots -- FISH'd genes

focal-NTRK1.pdf: plot_focal_pair.py $(clseqpair) $(clacghpair)
	python $^ NTRK1 -o $@

focal-ALK.pdf: plot_focal_pair.py $(clseqpair) $(clacghpair)
	python $^ ALK -o $@

focal-ROS1.pdf: plot_focal_pair.py $(clseqpair) $(clacghpair)
	python $^ ROS1 -o $@

focal-MET-BRAF.pdf: plot_focal_pair.py $(clseqpair) $(clacghpair)
	python $^ MET,BRAF -w 1e7 -o $@

focal-RET.pdf: plot_focal_pair.py $(clseqpair) $(clacghpair)
	python $^ RET -o $@

focal-CDKN2A.pdf: plot_focal_pair.py $(clseqpair) $(clacghpair)
	python $^ CDKN2A -w 2e6 -o $@

focal-chr9.pdf: $(clseqpair) $(clacghpair)
	cnvkit.py scatter $< -s $(cnvbuild)/CL_seq.cns -c chr9 -o tmp1.pdf
	cnvkit.py scatter $(cnvbuild)/CL_acgh.cnr -s $(cnvbuild)/CL_acgh.cns -c chr9 -o tmp2.pdf
	pdfnup --nup 2x1 tmp1.pdf tmp2.pdf
	mv tmp2-nup.pdf $@
	rm tmp1.pdf tmp2.pdf

# Global correlation

cell-correlation.pdf: plot_paired_segments.py cell-correlation.csv
	python $^ -o $@

cell-correlation.csv: pair_segments.py $(cnvbuild)/CL_acgh.cns $(cnvbuild)/CL_seq.cns
	python $^ -i $(baits) -o $@

cell-correlation-flat.pdf: plot_paired_segments.py cell-correlation-flat.csv
	python $^ -o $@

cell-correlation-flat.csv: pair_segments.py $(cnvbuild)/CL_acgh.cns $(cnvbuild)/CL_seq_flat.cns
	python $^ -i $(baits) -o $@

