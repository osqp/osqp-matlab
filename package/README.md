# Matlab Interface packaging functions

Make sure you have `cmake` on your path, and have cloned the entire source tree (including submodules) locally.

```
git submodule update --init --recursive
```

Then simply run `package_osqp.m` from within MATLAB. This will compile the interface and package it as a `osqp-matlab-<platform>64.tar.gz` file.
This can also be done on the command line:

```
/path/to/matlab -nodisplay -nosplash -nodesktop -r "cd package; package_osqp(); exit;"
```

Additionally, you can pass a version number to the `package_osqp` function. This is done automatically for the Linux
platform by Github actions (by looking at the release tag), but would have to be done manually for Windows/MacOS. For example,

```
/path/to/matlab -nodisplay -nosplash -nodesktop -r "cd package; package_osqp('0.6.2'); exit;"
```

Once the `.tar.gz` files for Windows/MacOS have been generated, upload them manually to the appropriate release as assets (Release -> Edit -> Upload files).
