#!/bin/bash
if [[ $target == *mingw* ]] ; then

# dopri
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5.o ./dopri5.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5_i32.o ./dopri5.f
gfortran -shared -o ./dopri5.dll ./dopri5.o
gfortran -shared -o ./dopri5_i32.dll ./dopri5_i32.o
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853.o ./dop853.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853_i32.o ./dop853.f
gfortran -shared -o ./dop853.dll ./dop853.o
gfortran -shared -o ./dop853_i32.dll ./dop853_i32.o

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

# dopri
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5.o ./dopri5.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5_i32.o ./dopri5.f
gfortran -shared -fPIC -o ./dopri5.dylib ./dopri5.o
gfortran -shared -fPIC -o ./dopri5_i32.dylib ./dopri5_i32.o
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853.o ./dop853.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853_i32.o ./dop853.f
gfortran -shared -fPIC -o ./dop853.dylib ./dop853.o
gfortran -shared -fPIC -o ./dop853_i32.dylib ./dop853_i32.o

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

# dopri
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5.o ./dopri5.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5_i32.o ./dopri5.f
gfortran -shared -fPIC -o ./dopri5.so ./dopri5.o
gfortran -shared -fPIC -o ./dopri5_i32.so ./dopri5_i32.o
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853.o ./dop853.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853_i32.o ./dop853.f
gfortran -shared -fPIC -o ./dop853.so ./dop853.o
gfortran -shared -fPIC -o ./dop853_i32.so ./dop853_i32.o

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

