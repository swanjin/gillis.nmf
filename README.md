# spa R package
This repository contains a R package for performing Successive Projection Algorithm (SPA). 
This R code is translated from [MATLAB code](https://gitlab.com/ngillis/nmfbook/-/blob/master/algorithms/separable%20NMF/SPA/SPA.m) in the Non-negative Matrix Factorization book (N. Gillis, 2020).

# Installing the package
To install the `spa` package first you need to install devtools:

```ruby
install.packages("devtools")
library(devtools)
install_github("swanjin/spa")
```

# Running the package

The main function in the spa package is SPA. To get minimal help:

```ruby
library(spa)
?SPA
```

# Example

```ruby
SPA(matrix(1:16,4),2,c('normalize'))
```
