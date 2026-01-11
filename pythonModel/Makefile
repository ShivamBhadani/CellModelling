CC = gcc
CFLAGS = -Wall -Wextra -O2
LIBS = -lm

# Source files
SOURCES = OCV_SoC.c cellESR.c cellUKF.c
HEADERS = OCV_SoC.h cellUKF.h cellESR.h
LIBS = -lm
CFLAGS += -Wall -Wextra -O2
LDFLAGS += -shared -fPIC
DEFINES = -DUKFTEST

# Targets
.PHONY: all clean test

all:
	$(CC) $(CFLAGS) -o libcell.so $(SOURCES) $(CFLAGS) $(LDFLAGS) $(LIBS)
test:
	$(CC) $(CFLAGS) -o test_cellESR.out $(SOURCES) $(LIBS) $(DEFINES)
	@echo "run ./test_cellESR.out to execute the test program"
clean:
	rm -f libcell.so test_cellESR.out