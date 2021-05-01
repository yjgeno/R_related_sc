### attributes of R objects:
#attributes()

### 6 data types: typeof() - what is the object's data type (low-level)?
# character: "a", "swc"
# numeric: 2, 15.5
# integer: 2L (the L tells R to store this as an integer)
# logical: TRUE, FALSE
# complex: 1+4i (complex numbers with real and imaginary parts)

###'Everything in R is an object' (base objects)
### data structures: class() - what kind of object is it (high-level)?
# atomic vector
# list
# matrix
# data frame
# factors

### Object-oriented objects:
# S3
# S4

# https://adv-r.hadley.nz/base-types.html


x <- 1L
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- 1
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- letters
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- TRUE
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- cars
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- cars[1]
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- cars[[1]]
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- matrix(cars)
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- new.env()
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- expression(1 + 1)
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- quote(y <- 1 + 1)
print(c(class(x), mode(x), storage.mode(x), typeof(x)))

x <- ls
print(c(class(x), mode(x), storage.mode(x), typeof(x)))
