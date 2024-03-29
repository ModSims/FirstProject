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

run_benchmarks:
	cd benchmarks && bash run_benchmarks.sh

benchmark_4_screens:
	cd ../tests && tmux new-session \; split-window -h \; split-window -v \; split-window -v \; select-layout even-horizontal \; send-keys -t 0 'time lid_driven_cavity_2d --imax 64 --jmax 64 --solver jacobi' C-m \; send-keys -t 1 'time lid_driven_cavity_2d --imax 64 --jmax 64 --solver conjugated_gradient' C-m \; send-keys -t 2 'time lid_driven_cavity_2d --imax 64 --jmax 64 --solver multigrid_jacobi' C-m \; send-keys -t 3 'time lid_driven_cavity_2d --imax 64 --jmax 64 --solver multigrid_pcg' C-m

gource:
	gource -s 0.1 -f --hide-usernames --auto-skip-seconds 1 --key --highlight-users --hide progress --file-idle-time 0 --max-files 0 --background-colour 000000 --font-size 24 --title "ModSims"

all: build_folder_check configure build install_python_library test

.PHONY: build_folder_check build configure clean

# Set the default target to build
.DEFAULT_GOAL := all
