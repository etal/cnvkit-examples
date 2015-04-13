# Paired reference

#  .CBS_1.csv -> .cns
#  .CNATable.csv -> .cnr

cnr=$(patsubst %.CNATable.csv,%.contra-pair.cnr,$(wildcard *.CNATable.csv))
cns=$(patsubst %.CBS_1.csv,%.contra-pair.cns,$(wildcard *.CBS_1.csv))


all: $(cnr) $(cns)

clean:
	rm -vf *.cns *.cnr

$(cns): %.contra-pair.cns: ../contra2cns.py %.CBS_1.csv
	python $^ -o $@

$(cnr): %.contra-pair.cnr: ../contra2cnr.py %.CNATable.csv
	python $^ -o $@

