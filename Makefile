build_folder_check:
	@if [ ! -d "cpp/build" ]; then \
		mkdir -p cpp/build; \
		echo "Created build folder."; \
	fi
	cd cpp/build && cmake ..

build: build_folder_check
	cd build && make

clean:
	@rm -rf build

.PHONY: build_folder_check build

# Set the default target to build
.DEFAULT_GOAL := build
