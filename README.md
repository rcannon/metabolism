
## Contents
[Installing and Running](#installation-and-running))
[Cleaning](#cleaning)
[Bash Scripts](#bash-scripts)
[Notes](#notes)


## Installation and Running

From the top-level directory in deception or marianas:
```
$ source env.sh
$ bash_scripts/setup_prereqs.sh
$ bash_scripts/build_metabolism.sh
$ cd build
$ ./run
```

## Cleaning

Make sure you are out of the build directory. Then run
```
$ bash_scripts/clean_all.sh
```

## Bash Scripts

These are located in the `bash_scripts` folder under the top-level directory.
Feel free to call them from anywhare in the file system after you have set `PROJECT_DIR to be the top level directory of
the project (see the section on [Installing and Running](#installation-and-running)).

- `setup_prereqs.sh` - sets up Eigen and IFOPT from the provided git submodules.
- `build_metabolism.sh` - compiles the metabolism code using cmake.
- `clean_all.sh` - undoes everything from `setup_prereqs.sh` and `build_metabolism.sh`. Make sure you are not in the `build` directory when running this.

## Notes

1. When running the solver, IFOPT will probably say that
3 of the relaxed_regulation_upper_constraints are violated.
This can be fixed by increasing `Mb` (aka `big_M_value` in `maximum_entropy_relaxed.cc`)
to at least 200 (it is currently set to 50). However, in my experience it is not necessary
to do this in order to get convergence to the same value as the python version.

2. Assorted notes on working with Eigen and ifopt can be found in the `notes.txt` file in the same directory as this README.

