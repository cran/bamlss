# bamlss 1.1-2

* Bug fixing.

* Removed bit package depends.

* More stable implementation of the batchwise backfitting optimizer.

# bamlss 1.1-1

* Bug fixing.

* Support new features for generating weights for neural networks.

# bamlss 1.1-0

* New features.

* Bug fixing.

* New paper introducing bamlss at <https://arxiv.org/abs/1909.11784>

* In addition, there is now a new website <http://www.bamlss.org/> with technical examples, textbook
  examples, as well as examples from publications. The website will be expanded in the future and
  will serve as the main source for presenting new features of the bamlss package.

* Adding a vignette to CRAN which links to the new website.


# bamlss 1.0-2

* Bug fixing.

* First neural network implementations, the `n()` constructor.

* Support for very large data sets using the ff and ffbase package.

* Experimental batchwise backfitting algorithm `bbfit()`.


# bamlss 1.0-1

* Bug fixing.

* More model fitting engines: `lasso()`, `stabsel()`


# bamlss 0.1-2

* Fixed some issues with `tx()` for BayesX.

* Added `tx3()` for BayesX.

* Fixed environment saving problem when using `light = TRUE` in `bamlss()`.

* Dropped some parts of the `jm.mode()` return value.

* Solved `PROTECT()` problems in the C code.

* Many other small edits ...
