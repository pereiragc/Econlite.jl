using Econlite
using Test
using ForwardDiff

@testset "Econlite.jl" begin


    # utility.jl
    tol = 1e-7
    l = LogUtility()
    crra = IsoElastic(0.42)
    cara = CARA(1.2)

    for ps in (l, crra, cara)
        c = 1.23
        u = util(c, ps)
        cc = invutil(u, ps)

        @test abs(c - cc) < tol

        du=dutil(c, ps)
        ccc=invdutil(du, ps)

        @test abs(c - ccc) < tol


        # check with forward diff
        du_fd = ForwardDiff.derivative(t -> util(t, ps), c)

        @test abs(du_fd - du) < tol
    end


    # aggregators.jl

end
