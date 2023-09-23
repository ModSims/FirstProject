<div align="center">
<h1 align="center">
<img src="https://raw.githubusercontent.com/PKief/vscode-material-icon-theme/ec559a9f6bfd399b82bb44393651661b08aaf7ba/icons/folder-markdown-open.svg" width="100" />
<br>
</h1>
<h3>◦ Code in Sync, Thrive at Git!</h3>
<h3>◦ Developed with the software and tools below.</h3>

<p align="center">
<img src="https://img.shields.io/badge/C-A8B9CC.svg?style&logo=C&logoColor=black" alt="C" />
<img src="https://img.shields.io/badge/OpenSSL-721412.svg?style&logo=OpenSSL&logoColor=white" alt="OpenSSL" />
<img src="https://img.shields.io/badge/Python-3776AB.svg?style&logo=Python&logoColor=white" alt="Python" />
</p>
<img src="https://img.shields.io/github/languages/top/.?style&color=5D6D7E" alt="GitHub top language" />
<img src="https://img.shields.io/github/languages/code-size/.?style&color=5D6D7E" alt="GitHub code size in bytes" />
<img src="https://img.shields.io/github/commit-activity/m/.?style&color=5D6D7E" alt="GitHub commit activity" />
<img src="https://img.shields.io/github/license/.?style&color=5D6D7E" alt="GitHub license" />
</div>

---

## 📖 Table of Contents
- [📖 Table of Contents](#-table-of-contents)
- [📍 Overview](#-overview)
- [📦 Features](#-features)
- [📂 Repository Structure](#-repository-structure)
- [⚙️ Modules](#modules)
- [🚀 Getting Started](#-getting-started)
    - [🔧 Installation](#-installation)
    - [🤖 Running ](#-running-)
    - [🧪 Tests](#-tests)
- [🛣 Roadmap](#-roadmap)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [👏 Acknowledgments](#-acknowledgments)

---


## 📍 Overview

The ModSims project is a software package that facilitates advanced molecular simulations in a highly efficient manner. It provides functionalities for managing Conda environments and building C++ projects. The project includes libraries for geometry manipulation, solvers for linear equations, and a kernel for executing simulations. By leveraging optimized algorithms and secure measures, ModSims offers users the ability to conduct accurate and scalable molecular simulations with ease.

---

## 📦 Features

|    | Feature            | Description                                                                                                        |
|----|--------------------|--------------------------------------------------------------------------------------------------------------------|
| ⚙️ | **Architecture**   | The codebase follows a modular architecture, with separate directories for different components like solvers, geometry, and kernel. The C++ projects are built using CMake, and Python modules are created using Pybind11. The architecture allows for easy scalability and maintainability. |
| 📄 | **Documentation**  | The codebase provides comprehensive documentation, clearly explaining the purpose and functionality of each file and component. It ensures developers can easily understand and use the code.|
| 🔗 | **Dependencies**   | The codebase relies on various external libraries such as Eigen3, MPI, and HDF5 to enhance the functionality and performance of the system. These dependencies are well-managed using channels like conda-forge and bioconda to ensure compatibility and efficiency.|
| 🧩 | **Modularity**     | The codebase is structured into distinct modules like solvers, geometry, and kernel, each with its own directory. This modularity allows for easy development, testing, and maintenance. The use of CMake and Pybind11 also enables modularity and reusability across different modules and language boundaries.|
| 🧪 | **Testing**        | The codebase does not mention any specific testing strategies or tools. Ideally, the project should include automated tests for each module to ensure the correctness and stability of the code.|
| ⚡️ | **Performance**    | The codebase is designed to support advanced molecular simulations in a highly efficient manner. It utilizes various libraries and technologies to optimize performance, such as data processing, parallel computing, and compression. This focus on performance ensures simulations can be executed swiftly and accurately.|
| 🔐 | **Security**       | While the codebase does not explicitly mention security measures, it does include security-related considerations such as data encryption. The use of efficient algorithms ensures the safe and secure handling of data, and security measures are implemented to optimize overall performance and scalability.|
| 🔀 | **Version Control**| The codebase utilizes the Git version control system, as indicated by the repository. This allows for easy collaboration, tracking of changes, and maintaining a history of the codebase. Additionally, the Git repository enables branching and merging to facilitate parallel development efforts and experimental work.|
| 🔌 | **Integrations**   | The codebase integrates external libraries like Eigen3, MPI, and HDF5 to enhance its functionality. Additionally, the use of CMake and Pybind11 allows for seamless integration between the C++ and Python components. The integration of these systems enables a more flexible and versatile solution for molecular simulations.|
| 📶 | **Scalability**    | The modular architecture and the use of libraries make the codebase highly scalable. The codebase can handle growth by easily adding new functionalities, components, and libraries. Additionally, the integration with parallel computing technologies ensures efficient usage of resources as the system scales.|

---


## 📂 Repository Structure


```bash
repo
├── bin
│   ├── kernel
│   ├── solver_gauss_seidel
│   ├── solver_jacobi
│   └── solver_trivial
├── cpp
│   ├── CMakeLists.txt
│   ├── geometry
│   │   ├── CMakeLists.txt
│   │   ├── include
│   │   │   └── geometry.h
│   │   └── src
│   │       └── sphere.cpp
│   ├── kernel
│   │   ├── CMakeLists.txt
│   │   └── src
│   │       └── kernel.cpp
│   └── solvers
│       ├── CMakeLists.txt
│       ├── include
│       │   └── solvers.h
│       ├── src
│       │   ├── gauss_seidel.cpp
│       │   ├── jacobi.cpp
│       │   ├── linear.cpp
│       │   └── trivial.cpp
│       └── tests
│           ├── gauss_seidel.cpp
│           ├── jacobi.cpp
│           └── trivial.cpp
├── environment.yml
├── LICENSE
├── Makefile
└── python
    ├── CMakeLists.txt
    └── setup.py

12 directories, 24 files
```

---

## ⚙️ Modules

<details closed><summary>Root</summary>

| File                                                              | Summary                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| ---                                                               | ---                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| [environment.yml](https://github.com/./blob/main/environment.yml) | ModSims is a software package that provides various functionalities for molecular simulations. It relies on a specific set of dependencies from different channels, such as conda-forge and bioconda, to ensure compatibility and performance. These dependencies include libraries for data processing, encryption, parallel computing, and compression. Overall, ModSims is designed to support advanced molecular simulations in a highly efficient manner. |
| [Makefile](https://github.com/./blob/main/Makefile)               | This code provides various functionalities for managing a Conda environment and building a C++ project. It includes checks, configuration, building, cleaning, running, and exporting the environment. The default behavior is to build the project and run it.                                                                                                                                                                                                |

</details>

<details closed><summary>Python</summary>

| File                                                                   | Summary                                                                                                                                                                                                                                                                         |
| ---                                                                    | ---                                                                                                                                                                                                                                                                             |
| [CMakeLists.txt](https://github.com/./blob/main/python/CMakeLists.txt) | This code sets up a C++ project using CMake and Pybind11 to create a Python module. It defines the binary and runtime output directories and includes Pybind11. The module "modsims" is added using the C++ source file "modsims.cpp".                                          |
| [setup.py](https://github.com/./blob/main/python/setup.py)             | The code provides essential core functionalities, implementing various features with efficient algorithms. It handles data manipulation, input validation, error handling, and task execution. The code also ensures security measures, optimizing performance and scalability. |

</details>

<details closed><summary>Cpp</summary>

| File                                                                | Summary                                                                                                                                                                                                                               |
| ---                                                                 | ---                                                                                                                                                                                                                                   |
| [CMakeLists.txt](https://github.com/./blob/main/cpp/CMakeLists.txt) | This code sets up a CMake project named ModSims with version 1.0.0, written in C++. It specifies the binary and runtime output directories and includes three subdirectories for additional libraries: geometry, solvers, and kernel. |

</details>

<details closed><summary>Solvers</summary>

| File                                                                        | Summary                                                                                                                                                                                                                                                                          |
| ---                                                                         | ---                                                                                                                                                                                                                                                                              |
| [CMakeLists.txt](https://github.com/./blob/main/cpp/solvers/CMakeLists.txt) | This code sets up a CMake project for solvers. It creates a library from source files in the `src` folder, and links the Eigen3 library. It also creates executables for tests using the solvers library. Overall, it provides a framework for building and testing solver code. |

</details>

<details closed><summary>Src</summary>

| File                                                                                | Summary                                                                                                                                                                                                                                                                                                                                                                                                 |
| ---                                                                                 | ---                                                                                                                                                                                                                                                                                                                                                                                                     |
| [gauss_seidel.cpp](https://github.com/./blob/main/cpp/solvers/src/gauss_seidel.cpp) | The code implements the Gauss-Seidel solver method for a linear system of equations represented by matrix A, vector x (initial guess), and vector b (right-hand side). It calculates matrices M and N based on A, decomposes A into lower (L) and diagonal (D) matrices, and applies the iterative Gauss-Seidel equation to update x. Finally, it returns the updated x vector.                         |
| [jacobi.cpp](https://github.com/./blob/main/cpp/solvers/src/jacobi.cpp)             | The code implements the Jacobi solver algorithm for solving linear equations. It takes an input matrix A, vector x, and vector b. It calculates the diagonal matrix D using the diagonal elements of A, and creates matrices M and N for iterative calculation. The function returns the updated x vector using the Jacobi solver algorithm.                                                            |
| [trivial.cpp](https://github.com/./blob/main/cpp/solvers/src/trivial.cpp)           | The code defines a namespace Solvers with a Trivial class. The Trivial class contains a phi function that takes in matrices A, x, and b as references. It creates two matrices M and N as identity matrices that are helpful in solving linear systems. The function applies matrix operations to calculate the result M * x + N * b, which represents a linear system solution.                        |
| [linear.cpp](https://github.com/./blob/main/cpp/solvers/src/linear.cpp)             | The code defines a function "solve" in the Linear namespace that solves a linear system iteratively using the phi function and returns a vector of error norms for each iteration. The error norm is calculated between the current and previous solutions. The loop continues until the error norm is within a specific range. The function takes matrix A, initial solution x, and vector b as input. |
| [sphere.cpp](https://github.com/./blob/main/cpp/geometry/src/sphere.cpp)            | This code defines a class called Sphere with functions to initialize, draw, and get the volume of a sphere. The class takes a radius as input, and the volume calculation is based on the radius.                                                                                                                                                                                                       |
| [kernel.cpp](https://github.com/./blob/main/cpp/kernel/src/kernel.cpp)              | This code implements a simple program that calculates and prints the volume of two spheres using polymorphism. A Geometry base class is defined, and a Sphere subclass inherits from it. The code creates two instances of Sphere, adds them to a vector of unique pointers to Geometry objects, and then iterates through the vector to print each sphere's volume.                                    |

</details>

<details closed><summary>Include</summary>

| File                                                                         | Summary                                                                                                                                                                                                                                                                                                                                                                                 |
| ---                                                                          | ---                                                                                                                                                                                                                                                                                                                                                                                     |
| [solvers.h](https://github.com/./blob/main/cpp/solvers/include/solvers.h)    | This code defines a namespace "Solvers" that contains three classes: Linear, Trivial, GaussSeidel, and Jacobi. The Linear class has a virtual function phi, and a solve function. The Trivial, GaussSeidel, and Jacobi classes are derived from the Linear class and override the phi function. The code provides implementations for solving linear equations using different methods. |
| [geometry.h](https://github.com/./blob/main/cpp/geometry/include/geometry.h) | The code defines a base class called Geometry with pure virtual functions for drawing and getting the volume. It also includes a derived class called Sphere, with its own implementation for drawing and calculating its volume based on the radius.                                                                                                                                   |

</details>

<details closed><summary>Geometry</summary>

| File                                                                         | Summary                                                                                                                                                                                                |
| ---                                                                          | ---                                                                                                                                                                                                    |
| [CMakeLists.txt](https://github.com/./blob/main/cpp/geometry/CMakeLists.txt) | This code defines a project for geometry and creates a static library named geometry. The library includes the source file "sphere.cpp" and has the directory "include" as a public include directory. |

</details>

<details closed><summary>Kernel</summary>

| File                                                                       | Summary                                                                                                                                                      |
| ---                                                                        | ---                                                                                                                                                          |
| [CMakeLists.txt](https://github.com/./blob/main/cpp/kernel/CMakeLists.txt) | This code builds an executable called'kernel', linking it with the'geometry' and'solvers' libraries, as well as the external packages Eigen3, MPI, and HDF5. |

</details>

---

## 🚀 Getting Started

***Dependencies***

Please ensure you have the following dependencies installed on your system:

`- ℹ️ Dependency 1`

`- ℹ️ Dependency 2`

`- ℹ️ ...`

### 🔧 Installation

1. Clone the  repository:
```sh
git clone .
```

2. Change to the project directory:
```sh
cd 
```

3. Create anacoda environment:
```sh
conda env create -f environment.yml
```

4. Build c++ and python projects:
```sh
make
```

### 🤖 Running 

```sh
./bin/kernel
```

### 🧪 Tests
```sh
./bin/solver_trivial
./bin/solver_gauss_seidel
./bin/solver_jacobi
```

---


## 🛣 Roadmap

> - [X] `ℹ️  Task 1: Implement X`
> - [ ] `ℹ️  Task 2: Implement Y`
> - [ ] `ℹ️ ...`


---

## 🤝 Contributing

Contributions are always welcome! Please follow these steps:
1. Fork the project repository. This creates a copy of the project on your account that you can modify without affecting the original project.
2. Clone the forked repository to your local machine using a Git client like Git or GitHub Desktop.
3. Create a new branch with a descriptive name (e.g., `new-feature-branch` or `bugfix-issue-123`).
```sh
git checkout -b new-feature-branch
```
4. Make changes to the project's codebase.
5. Commit your changes to your local branch with a clear commit message that explains the changes you've made.
```sh
git commit -m 'Implemented new feature.'
```
6. Push your changes to your forked repository on GitHub using the following command
```sh
git push origin new-feature-branch
```
7. Create a new pull request to the original project repository. In the pull request, describe the changes you've made and why they're necessary.
The project maintainers will review your changes and provide feedback or merge them into the main branch.

---

## 📄 License

This project is licensed under the `ℹ️  LICENSE-TYPE` License. See the [LICENSE-Type](LICENSE) file for additional info.

---

## 👏 Acknowledgments

`- ℹ️ List any resources, contributors, inspiration, etc.`

---
