# A package for dominance analysis

| **Build Status** |
|:----------------:|
| [![Build Status][build-img]][build-url] |
| [![Codecov branch][codecov-img]][codecov-url] |

[build-img]: https://github.com/guberger/Dominance.jl/workflows/CI/badge.svg?branch=main
[build-url]: https://github.com/guberger/Dominance.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/guberger/Dominance.jl/coverage.svg?branch=main
[codecov-url]: http://codecov.io/github/guberger/Dominance.jl?branch=main

## Workflow

### To run examples
```julia
(@v1.6) pkg> activate ./examples
```
Run your favorite example by executing the corresponding file, or running
```julia
julia> include("./examples/example_filename.jl")
```
### To run tests outside of CI
```julia
(@v1.6) pkg> activate ./test
```
Define a variable `TestLocal`
```julia
julia> TestLocal = ""
```
Run your favorite test by executing the corresponding file, or running
```julia
julia> include("./test/test_filename.jl")
```
