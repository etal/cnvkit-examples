# Pooled reference

#  .CBS_1.csv -> .cns
#  .CNATable.csv -> .cnr

# cnr=$(patsubst %.CNATable.csv,%.contra-pool.cnr,$(wildcard *.CNATable.csv))
cns=$(patsubst %.CBS_1.csv,%.contra-pool.cns,$(wildcard *.CBS_1.csv))


all: $(cns)

clean:
	rm -vf *.cns *.cnr

$(cns): %.contra-pool.cns: ../contra2cns.py %.CBS_1.csv
	python $^ -o $@

# $(cnr): %.contra-pool.cnr: ../contra2cnr.py %.CNATable.csv
#         python $^ -o $@

