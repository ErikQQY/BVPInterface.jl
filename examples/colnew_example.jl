using BVPInterface
@BVPInterface.import_huge

# We look at the boundary value problem
#
#                       sin⁴(x)
#    ε⋅y'(x) = sin²(x)-λ──────  ;   y(-π/2)=y(π/2)=1
#                         y
# 
# Here ε is given (e.g. ε=0.1) and λ is an unknown parameter: λ'=0
#
# We use colsys/colnew and a "homotopy": starting with ε=1.0, using this
# solution as start-guess for ε=0.5, then using this as start guess
# for ε=0.2 and then ε=0.1

a, b = -pi/2, pi/2
orders = [1, 1,]
ζ = [a, b]

global ε = nothing 
global ε_old = nothing
global sol_old = nothing

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

opt = OptionsODE("example 6",
      OPT_BVPCLASS => 2, OPT_COLLOCATIONPTS => 7,
      OPT_RTOL => [1e-4, 1e-4], OPT_MAXSUBINTERVALS => 200)
xx = collect(LinRange(a,b, 400));

sol = nothing
ε = 1.0
guess = sol_old !== nothing ? sol_old : initial_guess    
sol, retcode, stats = colnew([a,b], orders, ζ, rhs, Drhs, bc, Dbc, guess ,opt)
zz = evalSolution(sol, xx)