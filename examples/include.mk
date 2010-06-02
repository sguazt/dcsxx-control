.PHONY: xmp-build

xmp_SOURCES := $(wildcard $(addsuffix /*.cpp,$(xmp_srcdirs)))
#xmp_OBJS := $(patsubst $(xmp_srcdirs)/%,$(xmp_buildtmpdir)/%,$(patsubst %.cpp,%.$(obj_ext),$(xmp_SOURCES)))
xmp_OBJS := $(patsubst $(xmp_srcdir)/%,$(xmp_buildtmpdir)/%,$(patsubst %.cpp,%.$(obj_ext),$(xmp_SOURCES)))
xmp_TARGETS := $(addprefix $(xmp_bindir)/,$(patsubst %.cpp,%,$(patsubst $(xmp_srcdir)/%,%,$(xmp_SOURCES))))


xmp-build: override CC=$(CXX)
xmp-build: $(xmp_OBJS) $(xmp_TARGETS)

#$(info TEST BUILDTMPDIR ==> $(xmp_buildtmpdir))
#$(info TEST BINDIR ==> $(xmp_bindir))

$(xmp_bindir)/%: $(xmp_buildtmpdir)/%.$(obj_ext)
	mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) -o $@ $<


## Source to Object rules

$(xmp_buildtmpdir)/%.$(obj_ext): $(xmp_srcdir)/%.cpp
	@echo "=== (Test) Compiling: $@ ==="
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ -c $<


## Header to Precompiled header rules

ifeq ($(use_pch),true)
$(xmp_buildtmpdir)/%.$(pch_ext): %.hpp
	@echo "=== (Test) Pre-compiling header: $@ ==="
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@
else
$(xmp_buildtmpdir)/%.$(pch_ext): ;
endif


## Source to Dependency rules

$(xmp_buildtmpdir)/%.d: %.cpp
	@echo "=== (Test) Creating dependencies file: $@ ==="
	@set -e; rm -f $@; \
		$(CXX) $(CPPFLAGS) $< > $@.$$$$; \
		sed ’s,\($*\)\.$(obj_ext)[ :]*,\1.$(obj_ext) $@ : ,g’ < $@.$$$$ > $@; \
		rm -f $@.$$$$
