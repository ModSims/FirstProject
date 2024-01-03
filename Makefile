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
	cd cpp/build && make -j

test: build_folder_check
	cd cpp/build && make test

help:
	cd cpp/build && make help

all_options:
	cd cpp/build && cmake -LAH ..

clean:
	@rm -rf out
	@rm -rf cpp/build

export_env:
	conda env export --no-builds | grep -v "prefix" > environment.yml

install_python_library:
	conda develop -u $(shell pwd)/out/lib/ || true
	conda develop $(shell pwd)/out/lib/

clone_vtk:
	@if [ ! -d "cpp/extern/vtk" ]; then \
		git clone --depth 1 --branch v9.1.0 git@github.com:Kitware/VTK.git cpp/extern/vtk; \
	fi

all: build_folder_check configure build install_python_library test

.PHONY: build_folder_check build configure clean

# Set the default target to build
.DEFAULT_GOAL := all
