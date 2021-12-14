# gillis.nmf R package
This repository contains a R package for performing Nonnegative Matrix Factorization algorithms. 
This R code is translated from [MATLAB code](https://gitlab.com/ngillis/nmfbook/-/blob/master/algorithms/separable%20NMF/SPA/SPA.m) in the Nonnegative Matrix Factorization book (N. Gillis, 2020).

# Installing the package
To install the `gillis.nmf` package first you need to install devtools:

```ruby
install.packages("devtools")
library(devtools)
install_github("swanjin/gillis.nmf")
```

# Running the package

The one of the functions in the `gillis.nmf` package is spa. To get minimal help:

```ruby
library(gillis.nmf)
?spa
```

# spa example

```ruby
spa(matrix(1:16,4),2,c('normalize'))
```
