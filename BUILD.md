# Building STAG

## Building the documentation

In order to build the STAG documentation, you should will need

- Doxygen: install with `apt-get install doxygen` or by downloading from the
  [doxygen website](https://www.doxygen.nl/download.html). Make sure that the doxygen
  executable is available on your path.
- Sphinx: install with `apt-get install python3-sphinx` or by following the instructions
  on the [sphinx website](https://www.sphinx-doc.org/en/master/usage/installation.html).
  Make sure that the `sphinx-build` executable is available on your path.
- The sphinx `breathe` extension. You can install this with `pip install breathe`.
  This should be installed for the default python version installed on your build system.
- The sphinx `furo` theme. You can install this with `pip install furo`.

Then, building the cmake target `Documentation` should build the documentation successfully.
For example, the following cmake commands should do this:

```bash
cmake -S . -B cmake-build-debug
cmake --build cmake-build-debug --target Documentation
```

The built documentation files are then available in `cmake-build-debug/docs/docs/sphinx/`.