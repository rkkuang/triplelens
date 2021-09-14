GCC	= /usr/bin/g++
CFLAGS	= -std=c++11
#CFLAGS	= -O3 -std=c++11 -Wall -Wextra -pedantic -DNDEBUG
BUILDDIR = build
BINDIR = bin
SRCDIR = src
TESTDIR = test
SRCEXT = cpp
DATADIR = data

SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
TESTS = $(shell find $(TESTDIR) -type f -name *.$(SRCEXT))


OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
TESTOBJS = $(patsubst $(TESTDIR)/%,$(BUILDDIR)/%,$(TESTS:.$(SRCEXT)=.o))


TARGET = testtriple
#TEST1 = triLightCurve

all:	$(BINDIR)/$(TARGET)
#all:	$(BINDIR)/$(TARGET) $(BINDIR)/$(TEST1)

$(BINDIR)/$(TARGET): $(OBJECTS) $(TESTOBJS)
	$(GCC) $(CFLAGS) $(OBJECTS) $(BUILDDIR)/$(TARGET).o  -o $(BINDIR)/$(TARGET)

#$(BINDIR)/$(TEST1): $(OBJECTS) $(TESTOBJS)
#	$(GCC) $(CFLAGS) $(OBJECTS) $(BUILDDIR)/$(TEST1).o  -o $(BINDIR)/$(TEST1)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(BUILDDIR) $(BINDIR) $(DATADIR)
	$(GCC) $(CFLAGS)  -c -o $@ $<

$(BUILDDIR)/%.o: $(TESTDIR)/%.$(SRCEXT)
	$(GCC) $(CFLAGS)  -c -o $@ $<

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(DATADIR):
	mkdir -p $(DATADIR)

clean:
	$(RM) -r $(BUILDDIR)/*
	$(RM) -r $(BINDIR)/*

.PHONY: clean
