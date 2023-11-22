SUBDIRS = mole_C++ examples_C++ tests_C++

all clean: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: all $(SUBDIRS)
.NOTPARALLEL:
