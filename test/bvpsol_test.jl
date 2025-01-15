@testset "Test BVPSOL" begin
    using BVPInterface
    @BVPInterface.import_huge

    ivpopt = OptionsODE("ivpoptions",
    OPT_RHS_CALLMODE => RHS_CALL_INSITU)

    opt = OptionsODE("BVPSOL example",
    OPT_RHS_CALLMODE => RHS_CALL_INSITU,
    OPT_MAXSTEPS     => 10,
    OPT_RTOL         => 1e-6,
    OPT_BVPCLASS     => 0,
    OPT_SOLMETHOD    => 0,
    OPT_IVPOPT       => ivpopt)

    tNodes = [0,5]
    xInit = [ 5.0  0.45938665265299; 1  1];
    odesolver = nothing

    function f(t,x,dx)
        dx[1] =  x[2]
        dx[2] = -x[1]
        return nothing
    end

    function bc(xa,xb,r)
        r[1] = xa[1] - 5
        r[2] = xb[1] - 0.45938665265299
        return nothing
    end

    (t,x,code,stats) = bvpsol(f,bc,tNodes,xInit,odesolver,opt)
    @test code ≥ 0
end