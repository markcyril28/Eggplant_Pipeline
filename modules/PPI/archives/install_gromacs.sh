#!/usr/bin/env bash
set -euo pipefail
set -e

#######################################################################
# GROMACS 2024 Installation Script for AMD MI210 GPU
#
# Tested with: ROCm 7.0.1 (LLVM 20), AdaptiveCpp (hipSYCL)
# Target GPU: AMD Instinct MI210 (gfx90a)
# 
# Note: This script applies a workaround for LLVM version mismatch:
#   - ROCm 7.x uses LLVM 20
#   - AdaptiveCpp's clang plugin is built against system LLVM 18
#   - The plugin is patched out since direct HIP compilation works
#######################################################################

# Setup logging - capture all output to log file
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs/console_logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/gromacs_install_$(date +%Y%m%d_%H%M%S).log"

# Redirect all output to both console and log file
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=============================================="
echo "GROMACS Installation Log"
echo "Started: $(date)"
echo "Log file: ${LOG_FILE}"
echo "=============================================="

# 1. Set ROCm environment (adjust ROCM_PATH if different)
export ROCM_PATH=${ROCM_PATH:-/opt/rocm}
export PATH="${ROCM_PATH}/bin:${ROCM_PATH}/llvm/bin:$PATH"
export LD_LIBRARY_PATH="${ROCM_PATH}/lib:${ROCM_PATH}/llvm/lib:$LD_LIBRARY_PATH"

# Verify ROCm is available
if [ ! -d "$ROCM_PATH" ]; then
    echo "ERROR: ROCm not found at $ROCM_PATH"
    echo "Please install ROCm first: https://rocm.docs.amd.com/en/latest/deploy/linux/installer/install.html"
    exit 1
fi

# 2. Create and activate conda environment
# shellcheck disable=SC1091

# Initialize conda for shell activation
source "${HOME}/miniconda3/etc/profile.d/conda.sh"

# Also initialize mamba if available
if [ -f "${HOME}/miniconda3/etc/profile.d/mamba.sh" ]; then
    source "${HOME}/miniconda3/etc/profile.d/mamba.sh"
fi

# Prefer mamba, fallback to conda
if command -v mamba &> /dev/null; then
    PKG_MGR="mamba"
    ACTIVATE_CMD="mamba activate"
    echo "Using mamba for package management"
else
    PKG_MGR="conda"
    ACTIVATE_CMD="conda activate"
    echo "Using conda for package management"
fi

# Define environment path explicitly to avoid confusion
ENV_PATH="${HOME}/miniconda3/envs/gromacs"

# Create the environment if it doesn't exist
if [ ! -d "$ENV_PATH" ]; then
    echo "Creating conda environment 'gromacs' at ${ENV_PATH}..."
    $PKG_MGR create -p "$ENV_PATH" python=3.11 -y
fi

# Activate using the explicit path
$ACTIVATE_CMD "$ENV_PATH"

# 3. Install dependencies
$PKG_MGR install -c conda-forge -y \
  cmake \
  fftw \
  openmpi \
  pkg-config \
  python \
  boost

# 4. Install AdaptiveCpp (hipSYCL)
# Check for existing working AdaptiveCpp installations
ACPP_CMAKE_DIR=""

# Check if acpp works in conda prefix (test by running acpp --version)
if [ -x "${CONDA_PREFIX}/bin/acpp" ]; then
    if "${CONDA_PREFIX}/bin/acpp" --version >/dev/null 2>&1; then
        echo "Working AdaptiveCpp found in ${CONDA_PREFIX}"
        ACPP_CMAKE_DIR="${CONDA_PREFIX}/lib/cmake/AdaptiveCpp"
    else
        echo "AdaptiveCpp in conda env is broken, will rebuild..."
        rm -rf "${CONDA_PREFIX}/include/AdaptiveCpp" \
               "${CONDA_PREFIX}/lib/libacpp"* \
               "${CONDA_PREFIX}/lib/cmake/AdaptiveCpp" \
               "${CONDA_PREFIX}/lib/cmake/hipSYCL" \
               "${CONDA_PREFIX}/lib/cmake/OpenSYCL" \
               "${CONDA_PREFIX}/etc/AdaptiveCpp" \
               "${CONDA_PREFIX}/bin/acpp"* \
               "${CONDA_PREFIX}/bin/syclcc"* 2>/dev/null || true
    fi
fi

if [ -z "$ACPP_CMAKE_DIR" ]; then
    echo "Installing AdaptiveCpp (hipSYCL) from source..."
    
    cd "$HOME"
    # Clean up previous failed build
    if [ -d "AdaptiveCpp" ]; then
        echo "Removing previous AdaptiveCpp source directory..."
        rm -rf AdaptiveCpp
    fi
    
    git clone --recurse-submodules https://github.com/AdaptiveCpp/AdaptiveCpp.git
    cd AdaptiveCpp
    
    # Use develop branch for LLVM 20 compatibility (ROCm 6.x uses LLVM 20)
    git checkout develop
    git pull origin develop
    
    mkdir -p build && cd build
    
    # Build AdaptiveCpp with ROCm's LLVM
    # Note: ROCm 7.x doesn't ship LLVM cmake config files, so we disable the clang plugin
    # which would require building against LLVM. The HIP backend works without the plugin.
    cmake .. \
        -DCMAKE_C_COMPILER="${ROCM_PATH}/llvm/bin/clang" \
        -DCMAKE_CXX_COMPILER="${ROCM_PATH}/llvm/bin/clang++" \
        -DCMAKE_CXX_FLAGS="-Wno-deprecated-declarations" \
        -DCMAKE_SHARED_LINKER_FLAGS="-L${ROCM_PATH}/llvm/lib -Wl,-rpath,${ROCM_PATH}/llvm/lib" \
        -DCMAKE_EXE_LINKER_FLAGS="-L${ROCM_PATH}/llvm/lib -Wl,-rpath,${ROCM_PATH}/llvm/lib" \
        -DCMAKE_INSTALL_PREFIX="${CONDA_PREFIX}" \
        -DDEFAULT_TARGETS="hip:gfx90a" \
        -DDEFAULT_CLANG="${ROCM_PATH}/llvm/bin/clang++" \
        -DWITH_ROCM_BACKEND=ON \
        -DWITH_CUDA_BACKEND=OFF \
        -DWITH_SSCP_COMPILER=OFF \
        -DWITH_OPENCL_BACKEND=OFF \
        -DWITH_LEVEL_ZERO_BACKEND=OFF
    
    make -j32
    make install
    
    # Patch acpp to skip loading the clang plugin (incompatible LLVM versions)
    # ROCm 7.x uses LLVM 20 but AdaptiveCpp's plugin builds against system LLVM 18
    echo "Patching acpp script to skip plugin loading (LLVM version mismatch workaround)..."
    sed -i 's/flags += \["-fpass-plugin=" + self._config.acpp_plugin_path\]/pass  # Skip plugin: LLVM version mismatch/g' "${CONDA_PREFIX}/bin/acpp"
    sed -i 's/flags += \["-fplugin=" + self._config.acpp_plugin_path\]/pass  # Skip plugin: LLVM version mismatch/g' "${CONDA_PREFIX}/bin/acpp"
    
    # Update config to use ROCm clang by default
    sed -i 's|"default-clang".*:.*"/usr/lib/llvm-18/bin/clang++"|"default-clang" : "/opt/rocm/llvm/bin/clang++"|' \
        "${CONDA_PREFIX}/etc/AdaptiveCpp/acpp-core.json"
    
    ACPP_CMAKE_DIR="${CONDA_PREFIX}/lib/cmake/AdaptiveCpp"
    echo "AdaptiveCpp installed successfully"
fi

echo "Using AdaptiveCpp from: ${ACPP_CMAKE_DIR}"

# Verify AdaptiveCpp works
echo "Testing AdaptiveCpp..."
if ! "${CONDA_PREFIX}/bin/acpp" --version; then
    echo "ERROR: AdaptiveCpp installation failed - acpp not working"
    exit 1
fi

# 5. Download GROMACS 2024 source (if not already present)
GROMACS_SRC="${HOME}/gromacs"
if [ ! -d "$GROMACS_SRC" ]; then
    cd "$HOME"
    git clone https://gitlab.com/gromacs/gromacs.git
fi

cd "$GROMACS_SRC"
git fetch --all
git checkout release-2024

# 6. Clean previous build and configure for AMD MI210
rm -rf build
mkdir -p build && cd build

# Use ROCm Clang compiler (required for AdaptiveCpp/hipSYCL)
export CC="${ROCM_PATH}/llvm/bin/clang"
export CXX="${ROCM_PATH}/llvm/bin/clang++"

cmake .. \
  -DGMX_GPU=SYCL \
  -DGMX_SYCL=ACPP \
  -DAdaptiveCpp_DIR="${ACPP_CMAKE_DIR}" \
  -Dhipsycl_DIR="${ACPP_CMAKE_DIR}" \
  -DHIPSYCL_TARGETS='hip:gfx90a' \
  -DGMX_GPU_FFT_LIBRARY=VkFFT \
  -DGMX_MPI=ON \
  -DGMX_OPENMP=ON \
  -DGMX_BUILD_OWN_FFTW=OFF \
  -DCMAKE_C_COMPILER="${CC}" \
  -DCMAKE_CXX_COMPILER="${CXX}" \
  -DCMAKE_INSTALL_PREFIX="${CONDA_PREFIX}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DSYCL_CXX_FLAGS_EXTRA="-DHIPSYCL_ALLOW_INSTANT_SUBMISSION=1"

# 7. Compile and install
make -j32
make install

# 8. Source environment and set library paths
# shellcheck disable=SC1091
source "${CONDA_PREFIX}/bin/GMXRC"

# Ensure libraries are in the library path
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${ROCM_PATH}/lib:${LD_LIBRARY_PATH}"

# Create activation script for future sessions
mkdir -p "${CONDA_PREFIX}/etc/conda/activate.d"
cat > "${CONDA_PREFIX}/etc/conda/activate.d/gromacs-env.sh" << 'EOF'
#!/bin/bash
export ROCM_PATH=${ROCM_PATH:-/opt/rocm}
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${ROCM_PATH}/lib:${LD_LIBRARY_PATH}"
if [ -f "${CONDA_PREFIX}/bin/GMXRC" ]; then
    source "${CONDA_PREFIX}/bin/GMXRC"
fi
EOF
chmod +x "${CONDA_PREFIX}/etc/conda/activate.d/gromacs-env.sh"

echo "Created activation script at ${CONDA_PREFIX}/etc/conda/activate.d/gromacs-env.sh"

# 9. Verify installation
echo ""
echo "=== Verifying GROMACS installation ==="
gmx --version

# Check if GPU is detected
echo ""
echo "=== Checking GPU detection ==="
gmx mdrun -nb gpu -h 2>&1 | head -20 || true

echo ""
echo "=============================================="
echo "Installation completed: $(date)"
echo "Log file saved to: ${LOG_FILE}"
echo "=============================================="
