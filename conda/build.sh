export TACS_DIR=${SRC_DIR}
#export CONDA_BUILD_SYSROOT=${HOME}/opt/MacOSX10.10.sdk
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

cp Makefile.in.info Makefile.in;
make default TACS_DIR=${TACS_DIR} \
             LAPACK_LIBS="${LAPACK_LIBS}" \
             METIS_INCLUDE=-I${PREFIX}/include/ METIS_LIB="-L${PREFIX}/lib/ -lmetis" \
             SO_LINK_FLAGS="${LIB_SLF}" SO_EXT=${SO_EXT};
mv ${TACS_DIR}/lib/libtacs.${SO_EXT} ${PREFIX}/lib;
${PYTHON} -m pip install --no-deps . -vv;
cd ${TACS_DIR}/extern/f5tovtk;
make default TACS_DIR=${TACS_DIR} \
             LAPACK_LIBS="${LAPACK_LIBS}" \
             METIS_INCLUDE=-I${PREFIX}/include/ METIS_LIB="-L${PREFIX}/lib/ -lmetis" \
             SO_LINK_FLAGS="${F5TOVTK_SLF}" SO_EXT=${SO_EXT};
mv ./f5tovtk ${PREFIX}/bin;

