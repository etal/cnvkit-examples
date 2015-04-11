# Compare segmented CN values by aCGH, DNAseq at targeted genes in tumor samples

int_tr=../intervals/targeted-baits.interval_list
int_ex=../intervals/exome-baits.interval_list

# Not all TR samples have matching aCGH
tr_180k=01 02 05 06 08 11 12 13 14 26 28 29 31 34 36 37 42 43 45 48 49 55 58 59 64 77 79
tr_244k=93 94 95 96 98 99
tr_num=$(tr_180k) $(tr_244k)
ex_num=11 12 19 35 37 40 46 49 53 71


all: tr ex

.PHONY: tr
tr: tr-cnvkit-pool.stats.csv tr-cnvkit-pair.stats.csv tr-cnvkit-flat.stats.csv tr-contra-pool.stats.csv tr-contra-pair.stats.csv

.PHONY: ex
ex: ex-cnvkit-pool.stats.csv ex-cnvkit-pair.stats.csv ex-cnvkit-flat.stats.csv ex-contra-pool.stats.csv ex-contra-pair.stats.csv


.PHONY: clean
clean:
	rm -vf tables/??_*_T.*.csv *.stats.csv *.alt.pdf *.diffs.dat


comparison.pdf: compare_methods.py \
	TR-CNVkit-pool.diffs.dat EX-CNVkit-pool.diffs.dat \
	TR-CNVkit-pair.diffs.dat EX-CNVkit-pair.diffs.dat \
	TR-CNVkit-flat.diffs.dat EX-CNVkit-flat.diffs.dat \
	TR-CONTRA-pool.diffs.dat EX-CONTRA-pool.diffs.dat \
	TR-CONTRA-pair.diffs.dat EX-CONTRA-pair.diffs.dat \
	TR-CopywriteR-pair.diffs.dat EX-CopywriteR-pair.diffs.dat
	python $^ -o $@


# CNVkit, pooled normal reference

cnvbuild=../build
tr_ckpool=$(patsubst %,tables/TR_%_T.cnvkit-pool.csv,$(tr_num))
ex_ckpool=$(patsubst %,tables/EX_%_T.cnvkit-pool.csv,$(ex_num))

tr-cnvkit-pool.stats.csv: alt.py $(tr_ckpool)
	python $^ -n TR-CNVkit-pool > $@

ex-cnvkit-pool.stats.csv: alt.py $(ex_ckpool)
	python $^ -n EX-CNVkit-pool > $@

$(tr_ckpool): tables/%.cnvkit-pool.csv: pair_segments.py acgh/%.cns $(cnvbuild)/%.cns
	python $^ -i $(int_tr) -o $@

$(ex_ckpool): tables/%.cnvkit-pool.csv: pair_segments.py acgh/%.cns $(cnvbuild)/%.cns
	python $^ -i $(int_ex) -o $@

# CNVkit, paired normal references

tr_ckpair=$(patsubst %,tables/TR_%_T.cnvkit-pair.csv,$(tr_num))
ex_ckpair=$(patsubst %,tables/EX_%_T.cnvkit-pair.csv,$(ex_num))

tr-cnvkit-pair.stats.csv: alt.py $(tr_ckpair)
	python $^ -n TR-CNVkit-pair > $@

ex-cnvkit-pair.stats.csv: alt.py $(ex_ckpair)
	python $^ -n EX-CNVkit-pair > $@

$(tr_ckpair): tables/%.cnvkit-pair.csv: pair_segments.py acgh/%.cns cnvkit-pair/%_pair.cns
	python $^ -i $(int_tr) -o $@

$(ex_ckpair): tables/%.cnvkit-pair.csv: pair_segments.py acgh/%.cns cnvkit-pair/%_pair.cns
	python $^ -i $(int_ex) -o $@


# CNVkit, flat reference

tr_ckflat=$(patsubst %,tables/TR_%_T.cnvkit-flat.csv,$(tr_num))
ex_ckflat=$(patsubst %,tables/EX_%_T.cnvkit-flat.csv,$(ex_num))

tr-cnvkit-flat.stats.csv: alt.py $(tr_ckflat)
	python $^ -n TR-CNVkit-flat > $@

ex-cnvkit-flat.stats.csv: alt.py $(ex_ckflat)
	python $^ -n EX-CNVkit-flat > $@

$(tr_ckflat): tables/%.cnvkit-flat.csv: pair_segments.py acgh/%.cns cnvkit-flat/%_flat.cns
	python $^ -i $(int_tr) -o $@

$(ex_ckflat): tables/%.cnvkit-flat.csv: pair_segments.py acgh/%.cns cnvkit-flat/%_flat.cns
	python $^ -i $(int_ex) -o $@


# CONTRA, pooled normal references

tr_conpool=$(patsubst %,tables/TR_%_T.contra-pool.csv,$(tr_num))
ex_conpool=$(patsubst %,tables/EX_%_T.contra-pool.csv,$(ex_num))

tr-contra-pool.stats.csv: alt.py $(tr_conpool)
	python $^ -n TR-CONTRA-pool > $@

ex-contra-pool.stats.csv: alt.py $(ex_conpool)
	python $^ -n EX-CONTRA-pool > $@

$(tr_conpool): tables/%.contra-pool.csv: pair_segments.py acgh/%.cns contra-pool/%.contra-pool.cns
	python $^ -i $(int_tr) -o $@

$(ex_conpool): tables/%.contra-pool.csv: pair_segments.py acgh/%.cns contra-pool/%.contra-pool.cns
	python $^ -i $(int_ex) -o $@


# CONTRA, paired normal references

tr_conpair=$(patsubst %,tables/TR_%_T.contra-pair.csv,$(tr_num))
ex_conpair=$(patsubst %,tables/EX_%_T.contra-pair.csv,$(ex_num))

tr-contra-pair.stats.csv: alt.py $(tr_conpair)
	python $^ -n TR-CONTRA-pair > $@

ex-contra-pair.stats.csv: alt.py $(ex_conpair)
	python $^ -n EX-CONTRA-pair > $@

$(tr_conpair): tables/%.contra-pair.csv: pair_segments.py acgh/%.cns contra-pair/%.contra-pair.cns
	python $^ -i $(int_tr) -o $@

$(ex_conpair): tables/%.contra-pair.csv: pair_segments.py acgh/%.cns contra-pair/%.contra-pair.cns
	python $^ -i $(int_ex) -o $@

