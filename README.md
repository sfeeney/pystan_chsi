# pystan_chsi
Cubic Hermite spline interpolation for the PyStan package

Author
 - Stephen Feeney

The repository, for now, consists of a Python module (defining an object that calculates the properties of the interpolant), a test file, and Stan model code which constructs the interpolant and samples the test model. Just type `python chsi_test.py` to test it out! And note the following dependencies:

 - [PyStan](https://pystan.readthedocs.io/en/latest/)
 - [Corner](https://corner.readthedocs.io/en/latest/)
