INCLUDES=-Iincludes/
DEFINES=

#CXX =  g++-4.2.1
FLAGS=-MMD -g -O0
CFLAGS=$(FLAGS) $(INCLUDES) $(DEFINES)
CXXFLAGS=$(CFLAGS)

LDFLAGS=$(FLAGS)
LDLIBS=

TARGET=amgpp

SRCS=$(wildcard src/*.cpp)
OBJS=$(SRCS:%.cpp=%.o)
DEPS=$(SRCS:%.cpp=%.d)
DEPS_EXIST=$(wildcard $(DEPS))

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

clean:
	rm -f $(TARGET) $(OBJS) $(DEPS_EXIST)

depend: $(DEPS)

%.d: %.cpp
	g++ -E $(CXXFLAGS) -o /dev/null -MF $@ $<
# Actually -E should force gcc (in case of the -MD or -MMD switch) to
# interpret -o as the dependency output file and prevent the output of
# the preprocessed code. However this is not the case for -MMD which
# seems to be a bug.

release:
ifndef BUILDNAME
	$(error Error: BUILDNAME not given)
else ifneq ($(wildcard release/$(BUILDNAME)/),)
	$(error Error: Release with buildname $(BUILDNAME) already present)
else
	@if [ "$(RELEASEDIR)" == "" ]; then         \
		releasedir="./release";                 \
	else                                        \
		releasedir="$(RELEASEDIR)";             \
	fi;                                         \
	svn export . "$${releasedir}/$(BUILDNAME)"; \
	cd "$${releasedir}";                        \
	rm -rf "$(BUILDNAME)/data";                 \
	tar -cjf "$(BUILDNAME).tar.bz2" "$(BUILDNAME)/";
endif

.PHONY: all clean depend release

-include $(DEPS_EXIST)
# Including possibly non-existent dependency files would cause make to build
# them due to the implicit rule above. For efficiency reasons we prevent this
# by only including existing dependencies.
