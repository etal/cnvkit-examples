
tr-segment.seg ex-segment.seg: %.seg: ../copywriter2seg.R %.Rdata
	Rscript $^ $@

did-import-tr: tr-segment.seg
	cnvkit.py import-seg -p chr -c human $<
	rm -vf *_N.bam.vs.* 
	touch $@

