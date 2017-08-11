INSTALL_PATH?=/usr/local/CTBangBang/
SCRIPT_PATH?=/usr/bin/
SRC_DIR=$(shell pwd)/

all: rebuild

rebuild:
	mkdir -p src/obj
	$(MAKE) -C src ../icd

.PHONY: all clean

clean:
	$(MAKE) -C src clean
