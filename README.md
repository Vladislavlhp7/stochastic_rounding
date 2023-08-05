# Stochastic Rounding
## Description
This C implementation demonstrates the effect of Stochastic Rounding (SR) in software applied in Floating-Point (FP) arithmetic.

The experiments performed are:

* Repeated rounding: This experiment shows how SR can help to prevent the "stagnation" of rounding errors that can occur with deterministic rounding schemes.
* Harmonic series summation: This experiment shows how SR can improve the accuracy of the harmonic series summation, which is a notoriously difficult problem for floating-point arithmetic.
* Vector inner product simulation: This experiment shows how SR can improve the accuracy of the vector inner product, which is a common operation in many computational applications.

The numerical results provided show that SR can significantly reduce the absolute error of floating-point computations. This makes SR a promising rounding scheme for computationally-intensive tasks where accuracy is important.

## Usage
To run the experiments, simply compile the code and then run the main.c file. The results will be printed to the console.

To compile the code, run the following command in the terminal:

```
gcc -o main main.c
```

To run the code, run the following command in the terminal:

```
./main
```

## Documentation
For more information, please see the accompanying documentation.

I hope this is helpful! Let me know if you have any other questions.
