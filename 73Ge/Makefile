CXX = $(shell root-config --cxx)
CXXFLAGS = $(shell root-config --cflags) -std=c++17 -Wall -Wextra -O2 -fPIC
ROOTLIBS = $(shell root-config --libs)

SRCDIR = src
INCDIR = include
OBJDIR = obj
LIBDIR = lib

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
HEADERS = $(wildcard $(INCDIR)/*.hpp)

LIBNAME = lib-nm-toolkit
SHAREDLIB = $(LIBDIR)/$(LIBNAME).so
STATICLIB = $(LIBDIR)/$(LIBNAME).a

all: $(SHAREDLIB) $(STATICLIB)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(SHAREDLIB): $(OBJECTS) | $(LIBDIR)
	$(CXX) -shared -o $@ $^ $(ROOTLIBS)
	@echo "Built shared library: $(SHAREDLIB)"

$(STATICLIB): $(OBJECTS) | $(LIBDIR)
	ar rcs $@ $^
	@echo "Built static library: $(STATICLIB)"


clean:
	rm -rf $(OBJDIR) $(LIBDIR)

.PHONY: all clean
