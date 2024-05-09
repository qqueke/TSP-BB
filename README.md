# Traveling Salesman Problem using Branch and Bound

This project is an implementation of the Traveling Salesman Problem (TSP) using the Branch and Bound algorithm. It provides several implementations:

- **Serial:** Sequential execution of the Branch and Bound algorithm.
- **OpenMP:** Parallel implementation using OpenMP.
- **MPI:** Distributed memory parallel implementation using MPI.
- **Hybrid:** Combined parallel implementation using both OpenMP and MPI.

## Requirements

To compile and run this project, you will need:

- MPI compiler
- OpenMP  

## Compilation

To compile the project, use `make` followed by the implementation you want to compile:

- For Serial: `make` inside `serial` directory
- For OpenMP: `make` inside `omp` directory
- For MPI: `make` inside `mpi` directory
- For Hybrid: `make` inside `hybrid` directory

## Execution

After compiling, you can run the program by executing the generated binary:
```bash
./main <name_of_pub_instance> <max_bound>
```

- `<name_of_pub_instance>`: The name of the TSP instance you want to solve.
- `<max_bound>`: The maximum bound for the Branch and Bound algorithm.

## Output

The solutions for the TSP instances are stored in `<name_of_pub_instance>.out` file.

## Project Description

You can find the project description in the `Description.pdf` file.

For any additional information or issues, please refer to the project documentation or contact the project maintainers.
