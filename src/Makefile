LDFLAGS = -lm -pthread
CFLAGS = -g -Wall -Wextra -std=c99
CC = gcc

OBJS = swCost.o mem.o swCalc.o swOut.o swGetSeq.o

etu: swDnaLinear

aff : swDnaAffineLongGaps

prot : swProtAffine

all: swDnaLinear swDnaAffine swDnaAffineLongGaps swProtAffine

swDnaLinear: swDnaLinear.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

swDnaAffine: swDnaAffine.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

swDnaAffineLongGaps: swDnaAffineLongGaps.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

swProtAffine: swProtAffine.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -MMD -c -o $@ $<

# deps made thanks to -MMD above, just include them
-include $(SOURCES:.c=.d)

clean:
	rm -f *.o *.d *~
