# Calculate precision/recall of each tool (CNVkit, CopywriteR, CONTRA) versus
# CBS-segmented aCGH at large and small gains and losses.
# Overlap of CNV called regions is already calculated with `bedtools intersect`,
# with results in '.wao' and '.wbo' files (see ../compare/Makefile).

# wao = acgh as reference (left side), wbo = seq as reference

# tool x size x cnvtype -> precision, recall

tools := CL.contra-pair CL.contra-pool CL.cw-noref CL.cw-pair CL_flat CL_pair CL_seq 
tool_tables := $(addsuffix .tsv,$(tools))

all: alltools.pdf

alltools.pdf: plot_tools.py alltools.tsv
	python $^ -o $@

alltools.tsv: $(tool_tables)
	cat $^ | sort -ur > $@

$(tool_tables): %.tsv: pr_beds.py %.wao.bed %.wbo.bed
	python $^ -o $@

