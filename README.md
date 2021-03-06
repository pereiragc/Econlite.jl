# Econlite
Lightweight package providing very basic building blocks in computational economics. 

The reason for creating this package was that I wanted to avoid unnecessary repetition of things like: 

- Setting up a CRRA utility function 
- Setting up Cobb-Douglas/CES production function 

along with their derivatives and inverse functions of those derivatives. Plus, allow some more general statements like 

``` julia
income_tax = AdValorem(0.2)
w = 1.

net(w, income_tax)
```

instead of
``` julia
τ = 0.2
w = 1.

w * (1 - τ)
```

That allows you to easily reparametrize for a more complicated tax system. Same goes for aggregators, utility functions, etc.


## Related packages
Check out [QuantEcon.jl](https://github.com/QuantEcon/QuantEcon.jl). 


## Planned features
- [ ] Pretty printing
- [ ] Modified single-valued utilities (e.g., recentering, rescaling, general monotone transformations)
- [x] Basic markov chain types and methods
  - [x] simulate
  - [ ] estimate
- [x] Iteration utility 
- [ ] Derivatives/inverses of tax functions 
- [ ] Documentation



## Some examples

### Fixed point iteration
Let's solve the advanced equation 1 + 0.5x = x in a fancy way:

``` julia
x0 = [0.1]
F!(x, x0)=(x.=1 .+0.5x0)

contract!(x0, copy(x0), (), F!,
          IterationSpec(tol = 1e-12,
                        iterdisp=DisplayMod(2)))

x0 
```

Output: 

``` julia
julia> contract!(x0, copy(x0), (), F!,
                 IterationSpec(tol = 1e-12,
                               iterdisp=DisplayMod(2)))
Iteration 2
Iteration 4
Iteration 6
Iteration 8
Iteration 10
Iteration 12
Iteration 14
Iteration 16
Iteration 18
Iteration 20
Iteration 22
Iteration 24
Iteration 26
Iteration 28
Iteration 30
Iteration 32
Iteration 34
Iteration 36
Iteration 38
Iteration 40
Converged after 41 iterations

julia> x0
1-element Array{Float64,1}:
 1.999999999999136
```


### Markov Chain

``` julia
julia> import Random.seed!

julia> mc = MarkovChain([-0.5, 0.5], [
           0.1 0.9
           0.9 0.1
       ])
2-element MarkovChain{Float64}:
 -0.5
  0.5

julia> mc[1]
-0.5

julia> mc[2]
0.5

julia> Random.seed!(1234);

julia> mc[markov_genpath(mc, 1, 7)]
7-element Array{Float64,1}:
  0.5
 -0.5
  0.5
 -0.5
  0.5
 -0.5
  0.5
```

Also, you can "combine" two Markov chains independently: 

``` julia


julia> c = MarkovChain([true, false], [0.1  0.9;0.5  0.5])
2-element MarkovChain{Bool}:
 1
 0

julia> x = MarkovChain([-1., 1.], [0. 1.; 1. 0.])
2-element MarkovChain{Float64}:
 -1.0
  1.0

julia> combined = indeptensor(x, c)
4-element MarkovChainExtended{Tuple{Float64,Bool}}:
 (-1.0, 1)
 (1.0, 1)
 (-1.0, 0)
 (1.0, 0)

julia> transition(combined)
4×4 Array{Float64,2}:
 0.0  0.1  0.0  0.9
 0.1  0.0  0.9  0.0
 0.0  0.5  0.0  0.5
 0.5  0.0  0.5  0.0

```
