CXX := clang++

BUILD_DIR = build
SRC_FILE = src/main.cpp
EXE_FILE = $(BUILD_DIR)/main.exe

COMPILER_FLAGS := -o $(EXE_FILE) -Wall -pedantic-errors -Wextra -Wconversion -Wsign-conversion -std=c++17
COMPILER_FLAGS_DEBUG :=-g

all: debug

build_dir:
	@mkdir -p $(BUILD_DIR)/

debug: build_dir
	$(CXX) $(COMPILER_FLAGS) $(COMPILER_FLAGS_DEBUG) $(SRC_FILE)

run:
	$(EXE_FILE)

clean:
	@rm -rf $(BUILD_DIR)/
