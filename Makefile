# copyright: Michael Safyan


program_NAME := rtk2
program_C_SRCS := $(wildcard *.c)
program_CXX_SRCS := $(wildcard *.cpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS :=
program_LIBRARY_DIRS :=
program_LIBRARIES :=

CPPFLAGS += -std=c++0x -Wall -O3 -DnotRpackage=1
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += -pthread$(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))
LDLIBS += -lz

.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	echo "LINK: $(LINK.cc) "
	$(LINK.cc) $(program_OBJS) -o $(program_NAME)  $(LDLIBS) 


clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)
	if [ -d data/out/ ] && [ `ls data/out/ | wc -l` -ne 0 ]; then \
		echo "Removing data/out/*"; \
		rm -rf data/out/*; \
	fi

test: $(program_NAME)
	# check if data/table.tsv exists
	@if [ ! -f data/table.tsv ]; then \
		echo "Test failed: data/table.tsv does not exist"; \
		exit 1; \
	fi
	mkdir -p data/out/
	./$(program_NAME) -h 
	./$(program_NAME) memory -i data/table.tsv -o data/out/table. -ns
	# Check that data/out/table.lobal_diversity.tsv has 5 lines
	@if [ `cat data/out/table.global_diversity.tsv | wc -l` -eq 5 ]; then \
		echo "Test passed: 5 lines found in tableglobal_diversity.tsv"; \
	else \
		echo "Test failed: table.global_diversity.tsv does not have 5 lines"; \
		exit 1; \
	fi

	# Check if "F99" is in data/out/table.alpha_chao1.tsv 
	@if [ `grep -c "F99" data/out/table.alpha_chao1.tsv` -eq 1 ]; then \
		echo "Test passed: F99 found in table.alpha_chao1.tsv"; \
	else \
		echo "Test failed: F99 not found in table.alpha_chao1.tsv"; \
		exit 1; \
	fi
	

distclean: clean
