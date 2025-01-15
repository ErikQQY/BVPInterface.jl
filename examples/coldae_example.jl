using BVPInterface
@BVPInterface.import_huge

a, b = 0., 1.
orders = [1, 1, 1]
ζ = [a, a, b]
ny = 1

function rhs(x, z, y, f)
    f[1] = (1 + z[2] - sin(x)) * y[1] + cos(x)
    f[2] = cos(x)
    f[3] = y[1]
    f[4] = (z[1] - sin(x)) * (y[1] - exp(x))
end

function Drhs(x, z, y, df)
    df[:] .= 0.0
    df[1,2] = y[1]
    df[1,4] = (1 + z[2] - sin(x))
    df[3,4] = 1.0
    df[4,1] = y[1] - exp(x)
    df[4,4] = z[1] - sin(x)
end

function bc(i, z, bc)
    if i == 1
        bc[1] = z[1]
    elseif i == 2
        bc[1] = z[3] - 1.0
    elseif i == 3
        bc[1] = z[2] - sin(1.0)
    end
end

function Dbc(i, z, dbc)
    if i == 1
        dbc[1] = 1.0
        dbc[2] = 0.0
        dbc[3] = 0.0
    elseif i == 2
        dbc[1] = 0.0
        dbc[2] = 0.0
        dbc[3] = 1.0
    else
        dbc[1] = 0.0
        dbc[2] = 1.0
        dbc[3] = 0.0
    end
end

function initial_guess(x, z, y, dmz)
    z[1] = 0.5
    z[2] = 1.0
    z[3] = 1.0
    y[1] = 0.0
    rhs(x, z, y, dmz)
end

opt = OptionsODE("example 6",
      OPT_BVPCLASS => 2, OPT_COLLOCATIONPTS => 7,
      OPT_RTOL => [1e-4, 1e-4, 1e-4], OPT_MAXSUBINTERVALS => 200)
xx = collect(LinRange(a,b, 400));  
sol, retcode, stats = coldae([a,b], orders, ny, ζ, rhs, Drhs, bc, Dbc, nothing ,opt)