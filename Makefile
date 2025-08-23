# Makefile
# macOS-focused build for Pure Data externals
# Outputs: wg_line_tilde~.pd_darwin and wg_scatter_tilde~.pd_darwin
# Honors ARCHFLAGS (e.g., -arch x86_64 or -arch arm64)

PD_INCLUDE ?= ./pure-data/src
CC ?= clang

CFLAGS = -DPD -O3 -ffast-math -fPIC -Wall -Wextra -Wno-unused-parameter $(ARCHFLAGS) -I"$(PD_INCLUDE)"
LDFLAGS = -bundle -undefined dynamic_lookup $(ARCHFLAGS)

SOURCES = wg_line_tilde.c wg_scatter_tilde.c
TARGETS = wg_line_tilde~.pd_darwin wg_scatter_tilde~.pd_darwin

all: $(TARGETS)

wg_line_tilde~.pd_darwin: wg_line_tilde.o
	$(CC) $< -o $@ $(LDFLAGS)

wg_scatter_tilde~.pd_darwin: wg_scatter_tilde.o
	$(CC) $< -o $@ $(LDFLAGS)

wg_line_tilde.o: wg_line_tilde.c
	$(CC) $(CFLAGS) -c $< -o $@

wg_scatter_tilde.o: wg_scatter_tilde.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.pd_darwin

.PHONY: all clean

