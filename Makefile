CC=gcc

all: run

spline: spline.c
	$(CC) -o spline spline.c -lGL -lGLU -lglut

run: spline
	./spline

tags: spline.c 
	ctags spline.c

clean:
	rm -f spline tags
