#!/usr/bin/env python3

"""Generate the GNU, CUDA, and Intel ABACUS Dockerfiles.

The Dockerfiles use the ABACUS toolchain as the single source of dependency
installation logic. The generator only selects the base image, toolchain
configuration, and ABACUS CMake options.
"""

import argparse
import io
from pathlib import Path
from typing import Any


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="check that generated Dockerfiles are up to date",
    )
    args = parser.parse_args()

    with OutputFile("Dockerfile.gnu-openmpi", args.check) as f:
        f.write(
            install_deps_toolchain(
                base_image="ubuntu:24.04",
                mpi_mode="openmpi",
                math_mode="openblas",
                with_gcc="system",
                with_openmpi="install",
                with_openblas="install",
                with_libtorch="install",
                with_libnpy="install",
                with_dftd4="install",
                with_nep="install",
            )
        )
        f.write(install_abacus())

    with OutputFile("Dockerfile.gnu-mpich", args.check) as f:
        f.write(
            install_deps_toolchain(
                base_image="ubuntu:24.04",
                mpi_mode="mpich",
                math_mode="openblas",
                with_gcc="system",
                with_mpich="install",
                with_openblas="install",
                with_libtorch="install",
                with_libnpy="install",
                with_dftd4="install",
                with_nep="install",
            )
        )
        f.write(install_abacus())

# TODO: Check and re-generate Intel and CUDA testers
#
#    with OutputFile("Dockerfile.cuda", args.check) as f:
#        f.write(
#            install_deps_toolchain_cuda(
#                gpu_ver="80",
#                mpi_mode="openmpi",
#                math_mode="openblas",
#                with_gcc="system",
#                with_openmpi="install",
#                with_openblas="install",
#                with_dftd4="install",
#            )
#        )
#        f.write(install_abacus())

#    with OutputFile("Dockerfile.intel", args.check) as f:
#        f.write(
#            install_deps_toolchain(
#                base_image="intel/oneapi-toolkit:2026.1.0-devel-ubuntu24.04",
#                mpi_mode="intelmpi",
#                math_mode="mkl",
#                with_gcc="no",
#                with_intel="system",
#                with_intelmpi="system",
#                with_mkl="system",
#                with_ifx="yes",
#                with_libtorch="no",
#                with_dftd4="install",
#            )
#        )
#        f.write(install_abacus())


# ======================================================================================
def install_deps_toolchain(base_image: str, **kwargs: str) -> str:
    return f"FROM {base_image}\n\n" + install_toolchain(base_image, **kwargs)


# ======================================================================================
def install_deps_toolchain_cuda(gpu_ver: str, **kwargs: str) -> str:
    return rf"""
FROM nvidia/cuda:12.9.2-devel-ubuntu24.04

ARG GPU_VER={gpu_ver}

# Setup CUDA environment.
ENV CUDA_PATH=/usr/local/cuda
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64

""".lstrip() + install_toolchain(
        "ubuntu",
        enable_cuda="",
        gpu_ver="${GPU_VER}",
        **kwargs,
    )


# ======================================================================================
def install_toolchain(base_image: str, **kwargs: str) -> str:
    install_args = []
    for key, value in kwargs.items():
        option = key.replace("_", "-")
        if value == "":
            install_args.append(f"    --{option} \\")
        else:
            install_args.append(f"    --{option}={value} \\")
    install_args_str = "\n".join(install_args)

    return rf"""
# Install requirements for the toolchain.
WORKDIR /opt/abacus-toolchain
COPY ./toolchain/root_requirements/install_requirements*.sh ./
RUN ./install_requirements.sh {base_image}

# Configure the toolchain.
RUN mkdir scripts
COPY ./toolchain/scripts/VERSION \
     ./toolchain/scripts/package_versions.sh \
     ./toolchain/scripts/tool_kit.sh \
     ./toolchain/scripts/common_vars.sh \
     ./toolchain/scripts/signal_trap.sh \
     ./toolchain/scripts/get_openblas_arch.sh \
     ./toolchain/scripts/parse_if.py \
     ./scripts/
COPY ./toolchain/scripts/lib/ ./scripts/lib/
COPY ./toolchain/install_abacus_toolchain_new.sh .
RUN ./install_abacus_toolchain_new.sh \
{install_args_str}
    --skip-system-checks \
    --dry-run

# The dry run leaves the configuration files consumed by each stage script.
# Split the installation into stages so Docker can cache completed work.
COPY ./toolchain/scripts/stage0/ ./scripts/stage0/
RUN ./scripts/stage0/install_stage0.sh && rm -rf ./build

COPY ./toolchain/scripts/stage1/ ./scripts/stage1/
RUN ./scripts/stage1/install_stage1.sh && rm -rf ./build

COPY ./toolchain/scripts/stage2/ ./scripts/stage2/
RUN ./scripts/stage2/install_stage2.sh && rm -rf ./build

COPY ./toolchain/scripts/stage3/ ./scripts/stage3/
RUN ./scripts/stage3/install_stage3.sh && rm -rf ./build

COPY ./toolchain/scripts/stage4/ ./scripts/stage4/
RUN ./scripts/stage4/install_stage4.sh && rm -rf ./build
""".lstrip()


# ======================================================================================
def install_abacus(**cmake_options: str) -> str:

    return rf"""
# Additional prerequisites for building and testing
RUN apt-get -qq update && apt-get install -qq ninja-build bc python3-numpy

# Install ABACUS sources.
WORKDIR /opt/abacus
COPY ./source ./source
COPY ./tests ./tests
COPY ./cmake ./cmake
COPY ./CMakeLists.txt .

# Compile and install ABACUS.
COPY ./docker/scripts/build_abacus.sh ./build_abacus.sh
RUN ./build_abacus.sh

# Test ABACUS.
COPY ./docker/scripts/test_abacus.sh ./test_abacus.sh
RUN /bin/bash -o pipefail -c "./test_abacus.sh |& tee report.log"

# Output the report if the image is old and was therefore pulled from the build cache.
CMD cat $(find ./report.log -mmin +10) | sed '/^Summary:/ s/$/ (cached)/'
ENTRYPOINT []

# EOF
"""


# ======================================================================================
class OutputFile:
    def __init__(self, filename: str, check: bool) -> None:
        self.filename = filename
        self.check = check
        self.content = io.StringIO()
        self.content.write("#\n")
        self.content.write("# This file was created by generate_dockerfiles.py.\n")
        self.content.write("# Do not edit this file directly.\n")
        self.content.write(
            f"# Usage: podman build --shm-size=1g -f ./{filename} ..\n"
        )
        self.content.write("#\n\n")

    def __enter__(self) -> "OutputFile":
        return self

    def write(self, text: str) -> None:
        self.content.write(text)

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        if exc_type is not None:
            return

        output_path = Path(__file__).resolve().parent / self.filename
        generated = self.content.getvalue()

        if self.check:
            assert output_path.read_text(encoding="utf8") == generated
            print(f"File {output_path} is consistent with generator script.")
        else:
            output_path.write_text(generated, encoding="utf8")
            print(f"Wrote {output_path}")


# ======================================================================================
if __name__ == "__main__":
    main()

# EOF
