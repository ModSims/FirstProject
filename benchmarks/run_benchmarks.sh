#!/bin/bash

declare -A solvers=( ["jacobi"]="jacobi/"
                    ["conjugated_gradient"]="conjugated_gradient/"
                    ["multigrid_jacobi"]="multigrid_jacobi/"
                    ["multigrid_pcg"]="multigrid_pcg/" )

declare -a resolutions=("64" "128" "256")
declare -a REs=("100" "400" "1000")

run_solver() {
  local solver_name="$1"
  local resolution="$2"
  local RE="$3"
  local solver_dir="${RE}/${resolution}x${resolution}/lid_driven_cavity_2d/${solvers[$solver_name]}"
  local csv_file="stats.csv"

  rm -rf "$solver_dir"
  mkdir -p "$solver_dir"

  (
    cd "$solver_dir" || exit 1

    # Capture the start time before running the program
    start_time=$(date +%s.%N)

    # Run the solver and capture memory usage
    lid_driven_cavity_2d --imax "$resolution" --jmax "$resolution" --solver "$solver_name" --eps 0.000001 --Re "$RE" --omg 1.7 --t_end 50 &

    # Capture the process ID of the background job
    pid=$!

    # Monitor memory usage until the process finishes
    while ps -p $pid > /dev/null; do
      # Extract memory usage from /proc filesystem (in KB)
      memory_consumed=$(awk '/VmRSS/ {print $2}' /proc/$pid/status)
      sleep 1
    done

    # Capture the end time after the program finishes
    end_time=$(date +%s.%N)

    # Write CSV header if the file doesn't exist
    if [ ! -f "$csv_file" ]; then
      echo "Solver,Resolution,RE,Elapsed Time,Memory Consumed" > "$csv_file"
    fi

    # Append solver information to the CSV file
    elapsed_time=$(echo "$end_time - $start_time" | bc)
    echo "$solver_name,$resolution,$RE,$elapsed_time,$memory_consumed" >> "$csv_file"

    echo "$solver_name" > "${solver_name}_${resolution}.completed"
  )
}

# Run solvers for each resolution in parallel
for resolution in "${resolutions[@]}"; do
  for solver_name in "${!solvers[@]}"; do
    for RE in "${REs[@]}"; do
      echo "Running ${solvers[$solver_name]} for $resolution x $resolution with Re=$RE"
      run_solver "$solver_name" "$resolution" "$RE" &
    done
  done
  wait
done

# Wait for all background processes to finish
wait

# Print a final message
echo "All solvers done"
