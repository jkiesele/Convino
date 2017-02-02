CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
LD_FLAGS := `root-config --cflags --glibs` -lMinuit -lMathMore -lMinuit2 
CC_FLAGS := -fPIC -Wall `root-config --cflags`
CC_FLAGS += -I./include -O2 


ifeq ($(USE_MP),true)
	CC_FLAGS +=-fopenmp -DUSE_MP
else
	
endif

all: $(patsubst bin/%.cpp, %, $(wildcard bin/*.cpp)) libconvino.so




%: bin/%.cpp Makefile $(OBJ_FILES)
	g++ $(CC_FLAGS) $(LD_FLAGS) $(OBJ_FILES) $< -o $@ 

libconvino.so: $(OBJ_FILES)
	g++ -shared $(LD_FLAGS) -o $@ $^

obj/%.o: src/%.cpp
	g++ $(CC_FLAGS) -c -o $@ $<


clean: 
	rm -f obj/*.o obj/*.d
	touch bin/*
