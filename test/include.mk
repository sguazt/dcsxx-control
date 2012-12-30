.PHONY: test-build

test_SOURCES := $(wildcard $(addsuffix /*.cpp,$(test_srcdirs)))
#test_OBJS := $(patsubst $(test_srcdirs)/%,$(test_buildtmpdir)/%,$(patsubst %.cpp,%.$(obj_ext),$(test_SOURCES)))
test_OBJS := $(patsubst $(test_srcdir)/%,$(test_buildtmpdir)/%,$(patsubst %.cpp,%.$(obj_ext),$(test_SOURCES)))
test_TARGETS := $(addprefix $(test_bindir)/,$(patsubst %.cpp,%,$(patsubst $(test_srcdir)/%,%,$(test_SOURCES))))


test-build: override CC=$(CXX)
test-build: $(test_OBJS) $(test_TARGETS)

#$(info TEST BUILDTMPDIR ==> $(test_buildtmpdir))
#$(info TEST BINDIR ==> $(test_bindir))
#$(info TEST SRCDIRS ==> $(test_srcdirs))
#$(info TEST SOURCES ==> $(test_SOURCES))
#$(info TEST OBJS ==> $(test_OBJS))
#$(info TEST TARGETS ==> $(test_TARGETS))

$(test_bindir)/%: $(test_buildtmpdir)/%.$(obj_ext)
	mkdir -p $(dir $@)
	#$(CXX) $(LDFLAGS) -o $@ $<
	$(CXX) -o $@ $< $(LDFLAGS)


## Source to Object rules

$(test_buildtmpdir)/%.$(obj_ext): $(test_srcdir)/%.cpp
	@echo "=== (Test) Compiling: $@ ==="
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ -c $<


## Header to Precompiled header rules

ifeq ($(use_pch),true)
$(test_buildtmpdir)/%.$(pch_ext): %.hpp
	@echo "=== (Test) Pre-compiling header: $@ ==="
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@
else
$(test_buildtmpdir)/%.$(pch_ext): ;
endif


## Source to Dependency rules

$(test_buildtmpdir)/%.d: %.cpp
	@echo "=== (Test) Creating dependencies file: $@ ==="
	@set -e; rm -f $@; \
		$(CXX) $(CPPFLAGS) $< > $@.$$$$; \
		sed ’s,\($*\)\.$(obj_ext)[ :]*,\1.$(obj_ext) $@ : ,g’ < $@.$$$$ > $@; \
		rm -f $@.$$$$
