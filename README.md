
## Installation and Running

From the top-level directory in deception or marianas:
```
$ source bash_scripts/env.sh
$ bash_scripts/setup_prereqs.sh
$ bash_scripts/build_metabolism.sh
$ cd build
$ ./metabolism
```

## Cleaning

Make sure you are out of the build directory. Then run
```
$ bash_scripts/clean_all.sh
```

## Bash Scripts

These are located in the `bash_scripts` folder under the top-level directory.
Feel free to call them from anywhare in the file system; it shouldn't matter.

- `env.sh` : loads the necessary modules (gcc, cmake).
- `setup_prereqs.sh` - sets up Eigen and IFOPT from the provided git submodules.
- `build_metabolism.sh` - compiles the metabolism code using cmake.
- `clean_all.sh` - undoes everything from `setup_prereqs.sh` and `build_metabolism.sh`. Make sure you are not in the `build` directory when running this.


