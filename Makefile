CXX := g++
CXXFLAGS := -std=c++17 -O2

# ROOT configuration (required for calibration and data quality codes)
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

# Directories
CALIB_DIR := calibration/src
DATAQ_DIR := dataquality/src
BSD_DIR   := dataquality_withBSD/src
DIV_DIR   := divideevent/src

# Targets
CALIB_TARGETS := $(CALIB_DIR)/calibration $(CALIB_DIR)/convertlightyield
DATAQ_TARGET  := $(DATAQ_DIR)/dataqualityplot
BSD_TARGET    := $(BSD_DIR)/dataqualityplot_withBSD
DIV_TARGET    := $(DIV_DIR)/extract_events

ALL_TARGETS := $(CALIB_TARGETS) $(DATAQ_TARGET) $(BSD_TARGET) $(DIV_TARGET)

.PHONY: all clean

all: $(ALL_TARGETS)

# --- calibration ---

$(CALIB_DIR)/calibration: $(CALIB_DIR)/calibration.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $< -o $@ $(ROOTLIBS) -lSpectrum

$(CALIB_DIR)/convertlightyield: $(CALIB_DIR)/convertlightyield.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $< -o $@ $(ROOTLIBS) -lSpectrum

# --- dataquality ---

$(DATAQ_DIR)/dataqualityplot: $(DATAQ_DIR)/dataqualityplot.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $< -o $@ $(ROOTLIBS)

$(BSD_DIR)/dataqualityplot_withBSD: $(BSD_DIR)/dataqualityplot_withBSD.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $< -o $@ $(ROOTLIBS)

# --- divideevent ---

$(DIV_DIR)/extract_events: $(DIV_DIR)/extract_events.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(ALL_TARGETS)
