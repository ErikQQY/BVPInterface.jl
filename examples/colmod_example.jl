using BVPInterface
@BVPInterface.import_huge

Epsmin = 1e-9
ncomp = 1
orders = [2, ]
aleft = 0.0; aright = 1.0
zeta = [0.0, 1.0]

function rhs(x, z, f, e)
    f[1] = z[1]*(1.0-z[2])/e
end

function Drhs(x, z, df, ncomp, e)
    df[:] .= 0.0
    df[1,1] = (1.0-z[2])/e
    df[1,2] = -z[1]/e
end

function bc(i, z, g, e)
    if i == 1
        g[1] = z[1]-1.0
    elseif i == 2
        g[1] = z[1]+1.0/3.0
    end
end

function Dbc(i, z, dg, e)
    dg[1] = 1.0
    dg[2] = 0.0
end

function guess(x, z, dmval, e)
    z[1] = 1.0
end

opt = OptionsODE("COLMOD example",
      OPT_BVPCLASS => 1, OPT_COLLOCATIONPTS => 7,
      OPT_RTOL => [1e-8, 1e-4], OPT_MAXSUBINTERVALS => 200)
xx = collect(LinRange(aleft,aright,400));
sol, retcode, stats = colmod([aleft, aright], orders, 0.5, Epsmin, zeta, rhs, Drhs, bc, Dbc, guess ,opt)
zz = evalSolution(sol, xx)