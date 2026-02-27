#!/bin/bash
# Find and set DYLD_LIBRARY_PATH for all required JLLs

# Find the directories containing the libraries
GDAL_LIB_DIR=$(dirname $(find ~/.julia/artifacts -name "libgdal.37.dylib" | head -n 1))
OPENSSL_LIB_DIR=$(dirname $(find ~/.julia/artifacts -name "libcrypto.3.dylib" | head -n 1))
ZSTD_LIB_DIR=$(dirname $(find ~/.julia/artifacts -name "libzstd.1.dylib" | head -n 1))
LZ4_LIB_DIR=$(dirname $(find ~/.julia/artifacts -name "liblz4.1.dylib" | head -n 1))

# Construct the library path
export DYLD_LIBRARY_PATH=$GDAL_LIB_DIR:$OPENSSL_LIB_DIR:$ZSTD_LIB_DIR:$LZ4_LIB_DIR:$DYLD_LIBRARY_PATH

# Execute the command passed to the script
exec "$@"
