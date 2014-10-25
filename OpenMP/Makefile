# ----------------------- #
# MAKEFILE by Tiago Brito #
# ----------------------- #

CFLAGS=-fopenmp -Wall -ansi -pedantic -lm
CXX=g++

# ------------------------------------------------------------ #
# patsubst - substitui para todos os *.cpp - $(wildcard *.cpp) #
# a extensão %.cpp por %.o                                     #
# ------------------------------------------------------------ #

OBJECTS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))

# ------------------------------------------------------------ #
# wildcard - como não se pode colocar *.cpp numa destas flags  #
# pois isso seria interpretado como uma string usa-se o        #
# wildcard para fazer uma lista dos ficheiros .cpp que         #
# queremos						#
# ------------------------------------------------------------ #

PROGRAM=lcs_par

# ----------------------------------------- #
# $@ nome da regra que está antes dos :     #
# $^ tudo o que está depois dos :           #
# $< o primeiro nome que está depois dos :  #
# $> o último nome que está depois dos :    #
# ----------------------------------------- #

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CXX) $(CFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(PROGRAM)
