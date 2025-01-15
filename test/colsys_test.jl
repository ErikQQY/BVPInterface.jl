@testset "Test COLNEW" begin
    using BVPInterface
    @BVPInterface.import_huge
    
    a, b = -pi/2, pi/2
    orders = [1, 1,]
    ζ = [a, b]
    
    global ε = 1.0
    function rhs(x, z, f)
        s² = sin(x)^2
        f[1] = (s²-z[2]*s²*s²/z[1])/ε
        f[2] = 0.0
    end
    
    function Drhs(x, z, df)
        df[:] .= 0.0
        s⁴ = sin(x)^4
        df[1,1] = z[2]*s⁴/(z[1]^2)
        df[1,2] = -s⁴/z[1]
    end
    
    function bc(i, z, bc)
        bc[1] = z[1]-1.0
    end
    
    function Dbc(i, z, dbc)
        dbc[1] = 1.0
        dbc[2] = 0.0
    end
    
    function initial_guess(x, z, dmz)
        z[1] = 0.5
        z[2] = 1.0
        rhs(x, z, dmz)
    end
    
    opt = OptionsODE("COLSYS example",
          OPT_BVPCLASS => 2, OPT_COLLOCATIONPTS => 7,
          OPT_RTOL => [1e-4, 1e-4], OPT_MAXSUBINTERVALS => 200)
    xx = collect(LinRange(a,b, 400));


    sol, retcode, stats = colsys([a,b], orders, ζ, rhs, Drhs, bc, Dbc, initial_guess ,opt)
    @test retcode ≥ 0
end