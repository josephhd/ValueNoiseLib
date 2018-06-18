CC=gcc
CFLAGS= 

DEPS=valuenoise.h vectors.h
OUT=valuenoise

SRCS=valuenoisetest.c valuenoise.c
OBJS=$(SRCS:.c=.o)

#$@ refers to the name of the target
#$^ refers to the right of the target, or the dependencies
#$< refers to the first item in the dependencies

$(OUT): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c $< $(CFLAGS)

.PHONY: clean

clean:
	rm *.o $(OUT)

run:
	make
	./$(OUT)