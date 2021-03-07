CC=gcc
CFLAGS=-O0 -Wall -Wno-deprecated-declarations
ifeq ($(shell uname -s), Darwin)
	LDFLAGS=-lpthread -lm -framework OpenGL -framework GLUT
else
	LDFLAGS=-lpthread -lm -lGL -lGLU -lglut
endif
RM=rm -f

export CC CFLAGS LDFLAGS RM

all: gl2Body

gl2Body: gl2Body.o
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) gl2Body.o gl2Body
