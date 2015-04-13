
all: did-import-tr

clean:
	rm -vf *.cns did-import-tr

tr-segment.seg ex-segment.seg: %.seg: ../copywriter2seg.R %.Rdata
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
	# Short names for reference-free calls
	bash rename-noref.sh log2*.vs.none.cns
	mv -vf TR_101_.cw-noref.cns TR_101_T.cw-noref.cns
	touch $@

