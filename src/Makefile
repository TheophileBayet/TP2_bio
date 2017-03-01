LDFLAGS = -lm
CFLAGS = -g -Wall -Wextra -std=c99
CC = gcc

OBJS = swCost.o mem.o swCalc.o swOut.o swGetSeq.o

etu: swDnaLinear

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

