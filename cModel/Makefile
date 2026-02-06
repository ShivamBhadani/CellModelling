CC = gcc
CFLAGS = -Wall -Wextra -O2
LIBS = -lm

# Source files
SOURCES = OCV_SoC.c cellESR.c cellUKF.c lfs.c lfs_util.c
HEADERS = OCV_SoC.h cellUKF.h cellESR.h
LIBS = -lm
CFLAGS += -Wall -Wextra -g
LDFLAGS += -shared -fPIC
DEFINES ?= 

# Targets
.PHONY: all clean test

all:
	$(CC) $(CFLAGS) -o libcell.so $(SOURCES) $(CFLAGS) $(LDFLAGS) $(LIBS)
test:
	$(CC) $(CFLAGS) -o test_cellESR.out $(SOURCES) $(LIBS) -DUKFTEST $(DEFINES)
	@echo "run ./test_cellESR.out to execute the test program"
clean:
	rm -f libcell.so test_cellESR.out