check_conda_env:
	@if [ "$$CONDA_DEFAULT_ENV" != "ModSims" ]; then \
		echo "Error: The 'ModSims' environment is not active. Activate the 'ModSims' environment and try again."; \
		exit 1; \
	fi

build_folder_check: check_conda_env
	@if [ ! -d "cpp/build" ]; then \
		mkdir -p cpp/build; \
	fi

configure: build_folder_check
	cd cpp/build && cmake -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DEIGEN_TEST_NOQT=ON ..

build: build_folder_check
	cd cpp/build && make

help:
	cd cpp/build && make help

all_options:
	cd cpp/build && cmake -LAH ..

clean:
	@rm -rf bin
	@rm -rf cpp/build

run:
	./bin/kernel

export_env:
	conda env export --no-builds | grep -v "prefix" > environment.yml

all: build_folder_check configure build run

.PHONY: build_folder_check build configure clean run

# Set the default target to build
.DEFAULT_GOAL := all
