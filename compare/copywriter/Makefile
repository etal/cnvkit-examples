# NB: Don't run this in parallel (-j).

ex_seg=$(patsubst %-segment.Rdata,%_T.seg,$(wildcard EX_*-segment.Rdata))
ex_cns=$(ex_seg:.seg=.cw-pair.cns)

all: did-import-tr $(ex_cns)

clean:
	rm -vf *.cns did-import-tr

tr-segment.seg: %.seg: ../copywriter2seg.R %.Rdata
	Rscript $^ $@

$(ex_seg): %_T.seg: ../copywriter2seg.R %-segment.Rdata
	Rscript $^ $@

did-import-tr: tr-segment.seg
	sed 's/DM_/TR_/g' $< > tmp.seg
	cnvkit.py import-seg -p chr -c human tmp.seg
	rm tmp.seg
	# Remove normal sample calls
	rm -vf *_N.bam.vs.*
	# Short names for paired calls
	bash rename-pair.sh log2*.vs.log2.TR*.cns
	mv -vf TR_101_.cw-pair.cns TR_101_T.cw-pair.cns
	mv -vf TR_07_T.cw-pair.cns TR_08_T.cw-pair.cns
	# Short names for reference-free calls
	bash rename-noref.sh log2*.vs.none.cns
	mv -vf TR_101_.cw-noref.cns TR_101_T.cw-noref.cns
	mv -vf TR_07_T.cw-noref.cns TR_08_T.cw-noref.cns
	touch $@

$(ex_cns): %.cw-pair.cns: %.seg
	cnvkit.py import-seg -p chr -c human $<
	# Remove normal sample calls
	rm -vf *_N.bam.vs.*
	# Short names for paired calls
	bash rename-pair.sh log2*.vs.log2.EX*.cns
	# Short names for reference-free calls
	bash rename-noref.sh log2*.vs.none.cns
