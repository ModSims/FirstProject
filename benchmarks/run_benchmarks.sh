#!/bin/bash

declare -A solvers=( ["jacobi"]="jacobi/"
                    ["conjugated_gradient"]="conjugated_gradient/"
                    ["multigrid_jacobi"]="multigrid_jacobi/"
                    ["multigrid_pcg"]="multigrid_pcg/" )

declare -a resolutions=("64" "128" "256" "512")

run_solver() {
  local solver_name="$1"
  local resolution="$2"
  local solver_dir="${resolution}x${resolution}/lid_driven_cavity_2d/${solvers[$solver_name]}"

  (
    cd "$solver_dir" || exit 1
    yes | rm -rf * && lid_driven_cavity_2d --imax $resolution --jmax $resolution --solver "$solver_name" --eps 0.000001 --itermax $resolution --omg 1
    echo "$solver_name" > "${solver_name}_${resolution}.completed"
  )
}

# Run solvers for each resolution in parallel
for resolution in "${resolutions[@]}"; do
  for solver_name in "${!solvers[@]}"; do
    run_solver "$solver_name" "$resolution" &
  done
done

# Wait for all background processes to finish
wait

# Print a final message
echo "All solvers done"
