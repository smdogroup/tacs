export TACS_DIR=${SRC_DIR}

if [[ $(uname) == Darwin ]]; then
  export SO_EXT="dylib"
  export SO_LINK_FLAGS="-fPIC -dynamiclib"
  export LIB_SLF="${SO_LINK_FLAGS} -install_name @rpath/libtacs.dylib"
  export F5TOVTK_SLF="${SO_LINK_FLAGS} -install_name @rpath/f5tovtk"
  export LAPACK_LIBS="-framework accelerate"
elif [[ "$target_platform" == linux-* ]]; then
  export SO_EXT="so"
  export SO_LINK_FLAGS="-fPIC -shared"
  export LIB_SLF="${SO_LINK_FLAGS}"
  export F5TOVTK_SLF="${SO_LINK_FLAGS}"
  export LAPACK_LIBS="-L${PREFIX}/lib/ -llapack -lpthread -lblas"
fi

if [[ $scalar == "complex" ]]; then
  alias MAKE="make complex"
  alias PIP="CFLAGS=-DTACS_USE_COMPLEX ${PYTHON} -m pip"
elif [[ $mpi == "mpich" ]]; then
  alias MAKE="make default"
  alias PIP=" ${PYTHON} -m pip"
fi

cp Makefile.in.info Makefile.in;
MAKE TACS_DIR=${TACS_DIR} \
     LAPACK_LIBS="${LAPACK_LIBS}" \
     METIS_INCLUDE=-I${PREFIX}/include/ METIS_LIB="-L${PREFIX}/lib/ -lmetis" \
     SO_LINK_FLAGS="${LIB_SLF}" SO_EXT=${SO_EXT};
mv ${TACS_DIR}/lib/libtacs.${SO_EXT} ${PREFIX}/lib;
PIP install --no-deps --prefix=${PREFIX} . -vv;
cd ${TACS_DIR}/extern/f5tovtk;
# purposely lower-case make, only default supported below
make default TACS_DIR=${TACS_DIR} \
             LAPACK_LIBS="${LAPACK_LIBS}" \
             METIS_INCLUDE=-I${PREFIX}/include/ METIS_LIB="-L${PREFIX}/lib/ -lmetis" \
             SO_LINK_FLAGS="${F5TOVTK_SLF}" SO_EXT=${SO_EXT};
mv ./f5tovtk ${PREFIX}/bin;

