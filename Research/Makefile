SRCS := write.cpp
OBJS := $(SRCS:%.cpp=%.o)
EXEC := $(SRCS:%.cpp=%)
all: compile link
compile:
%.o: %.cpp
	g++ -I/home/mango/Software/eigen/ -I/home/mango/Software/pythia8309/include -I/home/mango/Software/root-6.26.10/montecarlo/pythia8/inc `root-config --cflags` -c $<

link: $(OBJS)
	g++ $(OBJS) -o $(EXEC) -lpythia8 -L/home/mango/Software/pythia8309/lib `root-config --glibs`
