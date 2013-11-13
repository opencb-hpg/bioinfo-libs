COMMONS_DIR = commons/

.PHONY: all clean dbwt csalib util

all: dbwt csalib util

dbwt:
	make -C dbwt/ all

csalib:
	make -C csalib/ all

util:
	make -C util/ all

clean:
	make -C $(COMMONS_DIR) clean
	make -C dbwt/ clean
	make -C csalib/ clean
	make -C search/ clean
	make -C gpu/ clean
	make -C util/ clean
