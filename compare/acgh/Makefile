# Import segmented array CGH datasets to CNVkit .cns format

all: 180k 244k 1mil

.PHONY: 180k 244k 1mil
180k 244k 1mil: %: aCGH_GenePattern_%.seg
	cnvkit.py import-seg --from-log10 -p chr -c human $<

.PHONY: clean
clean:
	rm -vf *.cns
