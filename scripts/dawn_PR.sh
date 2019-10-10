#!/bin/sh

exitError()
{
    echo "ERROR $1: $3" 1>&2
    echo "ERROR     LOCATION=$0" 1>&2
    echo "ERROR     LINE=$2" 1>&2
    exit $1
}

test -n "${target}" || exitError 303 ${LINENO} "Option <target> is not set"
export GRIDTOOLS_REPOSITORY=$(pwd)/gridtools #TODO adjust
export PYTUILS_REPOSITORY=$(pwd)/perf-benchmarks #TODO adjust

startdir=$(pwd)
workdir=/dev/shm/tmp_dawn

rm -rf ${workdir}
cp -r ${startdir} ${workdir}

cd ${workdir}/dawn

bash scripts/jenkins/build.sh
bash scripts/jenkins/build.sh

cd ${workdir}/gtclang
# TODO check we pass a target
if [ ${target} == "gpu" ]; then
  gpu_str="-g"
fi
if [ ${target} == "cuda" ]; then
  gpu_str="-c"
fi

dawn_path="-d ${workdir}/dawn/bundle/install/cmake"
bash scripts/jenkins/build.sh ${gpu_str} ${dawn_path}

cd ${workdir}
#git clone ${startdir}/clang-gridtools #TODO fix
cd clang-gridtools

if [ ${target} == "gpu" ]; then
  gpu_str="-g"
fi

if [ ${target} == "cuda" ]; then
  gpu_str="-u"
fi

bash scripts/jenkins/build.sh ${gpu_str} -c ${workdir}/gtclang/bundle/install/cmake

