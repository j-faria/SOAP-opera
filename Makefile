# Put the directory of 'DNest4/code' into this variable
DNEST4_PATH = DNest4/code

EIGEN_PATH = eigen

GSL_LIB := $(shell ./gsl/bin/gsl-config --libs)
GSL_INCLUDE = $(shell ./gsl/bin/gsl-config --cflags)

includes = -I$(DNEST4_PATH) -I$(EIGEN_PATH) $(GSL_INCLUDE)

CXX = g++-4.9
CXXFLAGS = -pthread -std=c++11 -O3 -DNDEBUG -w -DEIGEN_MPL2_ONLY
CC = gcc
CCFLAGS = -pthread -O3 -DNDEBUG -w 
# -DEIGEN_MPL2_ONLY

LIBS = -L$(DNEST4_PATH) -ldnest4 $(GSL_LIB)
LDFLAGS = -Wl,-rpath,$(shell ./gsl/bin/gsl-config --prefix)/lib

SRCDIR = ./src
SRCS =\
$(SRCDIR)/starspot.c \
$(SRCDIR)/Data.cpp \
$(SRCDIR)/ConditionalPrior.cpp \
$(SRCDIR)/Model.cpp \
$(SRCDIR)/main.cpp

OBJS1=$(subst .cpp,.o,$(SRCS))
OBJS=$(subst .c,.o,$(OBJS1))

HEADERS1=$(subst .cpp,.h,$(SRCS))
HEADERS=$(subst .c,.h,$(HEADERS1))

all: main
#all:
#	@echo $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(includes) $(CXXFLAGS)

%.o: %.c
	$(CC) -c -o $@ $< $(includes) $(CCFLAGS)


main: $(DNEST4_PATH)/libdnest4.a $(OBJS)
	$(CXX) -o main $(OBJS) $(CXXFLAGS) $(LIBS) $(LDFLAGS)

$(DNEST4_PATH)/libdnest4.a:
	make noexamples -C $(DNEST4_PATH)

clean:
	rm -f main $(OBJS)

cleanout:
	rm -f sample.txt sample_info.txt levels.txt weights.txt posterior_sample.txt sampler_state.txt

cleanall: clean
	rm -f sample.txt sample_info.txt \
		levels.txt weights.txt posterior_sample.txt sampler_state.txt \
		posterior_sample_lnlikelihoods.txt
