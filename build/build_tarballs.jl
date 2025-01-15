using BinaryBuilder

name = "BVPInterface"
version = v"0.0.1"

function get_autogenerated_build_script()
  return raw"""
    cd $WORKSPACE/srcdir/ODEInterfac*/src
    if [[ $target == *mingw* ]] ; then
    
    # lapack
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack.o ./dc_lapack.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack_i32.o ./dc_lapack.f
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack.o ./lapack.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack_i32.o ./lapack.f
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc.o ./lapackc.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc_i32.o ./lapackc.f
    
    # slatec
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec.o ./slatec.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec_i32.o ./slatec.f
    
    # bvpsol
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol.o ./bvpsol.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol_i32.o ./bvpsol.f
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol.o ./linalg_bvpsol.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol_i32.o ./linalg_bvpsol.f
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst.o ./zibconst.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst_i32.o ./zibconst.f
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol.o ./ma28_bvpsol.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol_i32.o ./ma28_bvpsol.f
    gfortran -shared -o ./bvpsol.dll ./bvpsol.o ./linalg_bvpsol.o ./zibconst.o ./ma28_bvpsol.o
    gfortran -shared -o ./bvpsol_i32.dll ./bvpsol_i32.o ./linalg_bvpsol_i32.o ./zibconst_i32.o ./ma28_bvpsol_i32.o
    
    # colnew
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew.o ./colnew.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew_i32.o ./colnew.f
    gfortran -shared -o ./colnew.dll ./colnew.o
    gfortran -shared -o ./colnew_i32.dll ./colnew_i32.o

    # coldae
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae.o ./coldae.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae_i32.o ./coldae.f
    gfortran -shared -o ./coldae.dll ./coldae.o
    gfortran -shared -o ./coldae_i32.dll ./coldae_i32.o

    # colmod
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colmod.o ./colmod.f
    gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colmod_i32.o ./colmod.f
    gfortran -shared -o ./colmod.dll ./colmod.o
    gfortran -shared -o ./colmod_i32.dll ./colmod_i32.o
    
    # bvpm2
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvp_la-2.o ./bvp_la-2.f
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -std=f2008 -o ./bvp_m-2.o ./bvp_m-2.f90
    gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wall -Wextra -Wimplicit-interface -std=f2008ts -o ./bvp_m_proxy.o ./bvp_m_proxy.f90
    gfortran -shared -o ./bvp_m_proxy.dll ./bvp_m_proxy.o ./bvp_m-2.o ./bvp_la-2.o
    cp *.dll $libdir
    
    
    
    elif [[ $target == *apple* ]] ; then
    
    # lapack
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack.o ./dc_lapack.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack_i32.o ./dc_lapack.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack.o ./lapack.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack_i32.o ./lapack.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc.o ./lapackc.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc_i32.o ./lapackc.f
    
    # slatec
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec.o ./slatec.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec_i32.o ./slatec.f
    
    # bvpsol
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol.o ./bvpsol.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol_i32.o ./bvpsol.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol.o ./linalg_bvpsol.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol_i32.o ./linalg_bvpsol.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst.o ./zibconst.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst_i32.o ./zibconst.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol.o ./ma28_bvpsol.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol_i32.o ./ma28_bvpsol.f
    gfortran -shared -fPIC -o ./bvpsol.dylib ./bvpsol.o ./linalg_bvpsol.o ./zibconst.o ./ma28_bvpsol.o
    gfortran -shared -fPIC -o ./bvpsol_i32.dylib ./bvpsol_i32.o ./linalg_bvpsol_i32.o ./zibconst_i32.o ./ma28_bvpsol_i32.o
    
    # colnew
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew.o ./colnew.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew_i32.o ./colnew.f
    gfortran -shared -fPIC -o ./colnew.dylib ./colnew.o
    gfortran -shared -fPIC -o ./colnew_i32.dylib ./colnew_i32.o

    # coldae
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae.o ./coldae.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae_i32.o ./coldae.f
    gfortran -shared -fPIC -o ./coldae.dylib ./coldae.o
    gfortran -shared -fPIC -o ./coldae_i32.dylib ./coldae_i32.o

    # colmod
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colmod.o ./colmod.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colmod_i32.o ./colmod.f
    gfortran -shared -fPIC -o ./colmod.dylib ./colmod.o
    gfortran -shared -fPIC -o ./colmod_i32.dylib ./colmod_i32.o
    
    # bvpm2
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvp_la-2.o ./bvp_la-2.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -std=f2008 -o ./bvp_m-2.o ./bvp_m-2.f90
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wall -Wextra -Wimplicit-interface -std=f2008ts -o ./bvp_m_proxy.o ./bvp_m_proxy.f90
    gfortran -shared -fPIC -o ./bvp_m_proxy.dylib ./bvp_m_proxy.o ./bvp_m-2.o ./bvp_la-2.o
    cp *.dylib $libdir
    
    
    
    else
    
    # lapack
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack.o ./dc_lapack.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack_i32.o ./dc_lapack.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack.o ./lapack.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack_i32.o ./lapack.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc.o ./lapackc.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc_i32.o ./lapackc.f
    
    # slatec
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec.o ./slatec.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec_i32.o ./slatec.f
    
    # bvpsol
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol.o ./bvpsol.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol_i32.o ./bvpsol.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol.o ./linalg_bvpsol.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol_i32.o ./linalg_bvpsol.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst.o ./zibconst.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst_i32.o ./zibconst.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol.o ./ma28_bvpsol.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol_i32.o ./ma28_bvpsol.f
    gfortran -shared -fPIC -o ./bvpsol.so ./bvpsol.o ./linalg_bvpsol.o ./zibconst.o ./ma28_bvpsol.o
    gfortran -shared -fPIC -o ./bvpsol_i32.so ./bvpsol_i32.o ./linalg_bvpsol_i32.o ./zibconst_i32.o ./ma28_bvpsol_i32.o
    
    # colnew
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew.o ./colnew.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew_i32.o ./colnew.f
    gfortran -shared -fPIC -o ./colnew.so ./colnew.o
    gfortran -shared -fPIC -o ./colnew_i32.so ./colnew_i32.o

    # coldae
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae.o ./coldae.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae_i32.o ./coldae.f
    gfortran -shared -fPIC -o ./coldae.so ./coldae.o
    gfortran -shared -fPIC -o ./coldae_i32.so ./coldae_i32.o

    # colmod
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colmod.o ./colmod.f
    gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colmod_i32.o ./colmod.f
    gfortran -shared -fPIC -o ./colmod.so ./colmod.o
    gfortran -shared -fPIC -o ./colmod_i32.so ./colmod_i32.o
    
    # bvpm2
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvp_la-2.o ./bvp_la-2.f
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -std=f2008 -o ./bvp_m-2.o ./bvp_m-2.f90
    gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wall -Wextra -Wimplicit-interface -std=f2008ts -o ./bvp_m_proxy.o ./bvp_m_proxy.f90
    gfortran -shared -fPIC -o ./bvp_m_proxy.so ./bvp_m_proxy.o ./bvp_m-2.o ./bvp_la-2.o
    cp *.so $libdir
    

    fi
  """
end

sources = [
  GitSource("https://github.com/ErikQQY/BVPInterface.jl", "27e1087bd70f5de9c78dadf3bd0a09a217af047b"),
]

script = get_autogenerated_build_script()

platforms = expand_gfortran_versions(supported_platforms())
# platforms = filter(p -> !isa(p, MacOS), platforms)
# platforms = [
#  BinaryBuilder.Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(libgfortran_version=v"5.0.0")),
# ]

products = [
  LibraryProduct("bvpsol",       :libbvpsol),
  LibraryProduct("bvpsol_i32",   :libbvpsol_i32),
  LibraryProduct("colnew",       :libcolnew),
  LibraryProduct("colnew_i32",   :libcolnew_i32),
  LibraryProduct("coldae",       :libcoldae),
  LibraryProduct("coldae_i32",   :libcoldae_i32),
  LibraryProduct("colmod",       :libcolmod),
  LibraryProduct("colmod_i32",   :libcolmod_i32),
  LibraryProduct("bvp_m_proxy",  :libbvp_m_proxy),
]

dependencies = [
  Dependency("CompilerSupportLibraries_jll")
]

build_tarballs(ARGS, name, version, sources, script,
               platforms, products, dependencies; preferred_gcc_version=v"9.1.0")