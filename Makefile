SRC_DIR = src
DOC_DIR = doc
DEPS_DIR = deps

.PHONY: doc test test-all

doc:
	$(MAKE) -C $(SRC_DIR)
	$(MAKE) -C $(DOC_DIR) html 

test:
	julia test/runtests.jl

test-all:
	julia test/runtests.jl all