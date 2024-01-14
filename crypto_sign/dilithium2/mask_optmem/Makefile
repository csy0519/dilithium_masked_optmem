# This Makefile can be used with GNU Make or BSD Make

LIB=libdilithium2_mask_optmem.a
HEADERS=api.h masked_api.h ntt.h packing.h params.h poly.h polyvec.h reduce.h rounding.h sign.h symmetric.h gadgets.h masked_poly.h masked_representations.h masked.h masked_utils.h masked_smallpoly.h
OBJECTS=ntt.o packing.o poly.o polyvec.o reduce.o rounding.o sign.o symmetric-shake.o gadgets.o masked_poly.o masked_representations.o masked_utils.o masked_smallpoly.o

CFLAGS=-O3 -Wall -Wextra -Wpedantic -Werror -Wmissing-prototypes -Wredundant-decls -std=c99 -I../../../common $(EXTRAFLAGS)

all: $(LIB)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(LIB): $(OBJECTS)
	$(AR) -r $@ $(OBJECTS)

clean:
	$(RM) $(OBJECTS)
	$(RM) $(LIB)
