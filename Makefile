CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -Wpedantic -fopenmp
CPPFLAGS = -Iinclude -Iexamples
IMPL = src/sparse.cpp
.PHONY: all test submission clean
all: build/demo build/bench build/check
build:
	mkdir -p build
build/sparse.o: $(IMPL) include/sparse.h | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(IMPL) -o $@
build/libsparse.a: build/sparse.o
	ar rcs $@ $<
build/demo: examples/demo.cpp examples/problems.h build/libsparse.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< build/libsparse.a -o $@
build/bench: examples/bench.cpp examples/problems.h build/libsparse.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< build/libsparse.a -o $@
build/check: tests/check.cpp build/libsparse.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< build/libsparse.a -o $@
test: build/check
	./build/check
submission: | build
	python3 -c 'import zipfile; z=zipfile.ZipFile("build/submission.zip","w",zipfile.ZIP_DEFLATED); z.write("src/sparse.cpp","sparse.cpp"); z.close()'
clean:
	rm -rf build
