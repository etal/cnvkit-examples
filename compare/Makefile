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
tr: tr-cnvkit-pool.diffs.dat tr-cnvkit-pair.diffs.dat tr-cnvkit-flat.diffs.dat \
	tr-contra-pool.diffs.dat tr-contra-pair.diffs.dat \
	tr-copywriter-pair.diffs.dat tr-copywriter-noref.diffs.dat

.PHONY: ex
ex: ex-cnvkit-pool.diffs.dat ex-cnvkit-pair.diffs.dat ex-cnvkit-flat.diffs.dat \
	ex-contra-pool.diffs.dat ex-contra-pair.diffs.dat \
	ex-copywriter-pair.diffs.dat ex-copywriter-noref.diffs.dat

.PHONY: clean
clean:
	rm -vf tables/??_*_T.*.csv *.stats.csv *.alt.pdf *.diffs.dat


comparison.pdf: compare_methods.py \
	tr-cnvkit-pool.diffs.dat ex-cnvkit-pool.diffs.dat \
	tr-cnvkit-pair.diffs.dat ex-cnvkit-pair.diffs.dat \
	tr-cnvkit-flat.diffs.dat ex-cnvkit-flat.diffs.dat \
	tr-contra-pool.diffs.dat ex-contra-pool.diffs.dat \
	tr-contra-pair.diffs.dat ex-contra-pair.diffs.dat \
	tr-copywriter-pair.diffs.dat ex-copywriter-pair.diffs.dat \
	tr-copywriter-noref.diffs.dat ex-copywriter-noref.diffs.dat
	python $^ -o $@ > comparison.stats.csv


# CNVkit, pooled normal reference

cnvbuild=../build
tr_ckpool=$(patsubst %,tables/TR_%_T.cnvkit-pool.csv,$(tr_num))
ex_ckpool=$(patsubst %,tables/EX_%_T.cnvkit-pool.csv,$(ex_num))

tr-cnvkit-pool.diffs.dat: alt.py $(tr_ckpool)
	python $^ -n tr-cnvkit-pool

ex-cnvkit-pool.diffs.dat: alt.py $(ex_ckpool)
	python $^ -n ex-cnvkit-pool

$(tr_ckpool): tables/%.cnvkit-pool.csv: pair_segments.py acgh/%.cns $(cnvbuild)/%.cns
	python $^ -i $(int_tr) -o $@

$(ex_ckpool): tables/%.cnvkit-pool.csv: pair_segments.py acgh/%.cns $(cnvbuild)/%.cns
	python $^ -i $(int_ex) -o $@

# CNVkit, paired normal references

tr_ckpair=$(patsubst %,tables/TR_%_T.cnvkit-pair.csv,$(tr_num))
ex_ckpair=$(patsubst %,tables/EX_%_T.cnvkit-pair.csv,$(ex_num))

tr-cnvkit-pair.diffs.dat: alt.py $(tr_ckpair)
	python $^ -n tr-cnvkit-pair

ex-cnvkit-pair.diffs.dat: alt.py $(ex_ckpair)
	python $^ -n ex-cnvkit-pair

$(tr_ckpair): tables/%.cnvkit-pair.csv: pair_segments.py acgh/%.cns cnvkit-pair/%_pair.cns
	python $^ -i $(int_tr) -o $@

$(ex_ckpair): tables/%.cnvkit-pair.csv: pair_segments.py acgh/%.cns cnvkit-pair/%_pair.cns
	python $^ -i $(int_ex) -o $@

# CNVkit, flat reference

tr_ckflat=$(patsubst %,tables/TR_%_T.cnvkit-flat.csv,$(tr_num))
ex_ckflat=$(patsubst %,tables/EX_%_T.cnvkit-flat.csv,$(ex_num))

tr-cnvkit-flat.diffs.dat: alt.py $(tr_ckflat)
	python $^ -n tr-cnvkit-flat

ex-cnvkit-flat.diffs.dat: alt.py $(ex_ckflat)
	python $^ -n ex-cnvkit-flat

$(tr_ckflat): tables/%.cnvkit-flat.csv: pair_segments.py acgh/%.cns cnvkit-flat/%_flat.cns
	python $^ -i $(int_tr) -o $@

$(ex_ckflat): tables/%.cnvkit-flat.csv: pair_segments.py acgh/%.cns cnvkit-flat/%_flat.cns
	python $^ -i $(int_ex) -o $@


# CONTRA, pooled normal references

tr_conpool=$(patsubst %,tables/TR_%_T.contra-pool.csv,$(tr_num))
ex_conpool=$(patsubst %,tables/EX_%_T.contra-pool.csv,$(ex_num))

tr-contra-pool.diffs.dat: alt.py $(tr_conpool)
	python $^ -n tr-contra-pool

ex-contra-pool.diffs.dat: alt.py $(ex_conpool)
	python $^ -n ex-contra-pool

$(tr_conpool): tables/%.contra-pool.csv: pair_segments.py acgh/%.cns contra-pool/%.contra-pool.cns
	python $^ -i $(int_tr) -o $@

$(ex_conpool): tables/%.contra-pool.csv: pair_segments.py acgh/%.cns contra-pool/%.contra-pool.cns
	python $^ -i $(int_ex) -o $@

# CONTRA, paired normal references

tr_conpair=$(patsubst %,tables/TR_%_T.contra-pair.csv,$(tr_num))
ex_conpair=$(patsubst %,tables/EX_%_T.contra-pair.csv,$(ex_num))

tr-contra-pair.diffs.dat: alt.py $(tr_conpair)
	python $^ -n tr-contra-pair

ex-contra-pair.diffs.dat: alt.py $(ex_conpair)
	python $^ -n ex-contra-pair

$(tr_conpair): tables/%.contra-pair.csv: pair_segments.py acgh/%.cns contra-pair/%.contra-pair.cns
	python $^ -i $(int_tr) -o $@

$(ex_conpair): tables/%.contra-pair.csv: pair_segments.py acgh/%.cns contra-pair/%.contra-pair.cns
	python $^ -i $(int_ex) -o $@


# CopywriteR, paired normal references

tr_cwpair=$(patsubst %,tables/TR_%_T.cw-pair.csv,$(tr_num))
ex_cwpair=$(patsubst %,tables/EX_%_T.cw-pair.csv,$(ex_num))

tr-copywriter-pair.diffs.dat: alt.py $(tr_cwpair)
	python $^ -n tr-copywriter-pair

ex-copywriter-pair.diffs.dat: alt.py $(ex_cwpair)
	python $^ -n ex-copywriter-pair

$(tr_cwpair): tables/%.cw-pair.csv: pair_segments.py acgh/%.cns copywriter/%.cw-pair.cns
	python $^ -i $(int_tr) -o $@

$(ex_cwpair): tables/%.cw-pair.csv: pair_segments.py acgh/%.cns copywriter/%.cw-pair.cns
	python $^ -i $(int_ex) -o $@

# CopywriteR, no reference

tr_cwnoref=$(patsubst %,tables/TR_%_T.cw-noref.csv,$(tr_num))
ex_cwnoref=$(patsubst %,tables/EX_%_T.cw-noref.csv,$(ex_num))

tr-copywriter-noref.diffs.dat: alt.py $(tr_cwnoref)
	python $^ -n tr-copywriter-noref

ex-copywriter-noref.diffs.dat: alt.py $(ex_cwnoref)
	python $^ -n ex-copywriter-noref

$(tr_cwnoref): tables/%.cw-noref.csv: pair_segments.py acgh/%.cns copywriter/%.cw-noref.cns
	python $^ -i $(int_tr) -o $@

$(ex_cwnoref): tables/%.cw-noref.csv: pair_segments.py acgh/%.cns copywriter/%.cw-noref.cns
	python $^ -i $(int_ex) -o $@

