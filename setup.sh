# Setup environment to run Obit software
# Modify OBIT_ROOT to point to your installed Obit binary distribution
export OBIT_ROOT=/scratch2/tmauch/refpointingJune23/obit-distro-1.1.648
export OBIT=${OBIT_ROOT}/lib/obit
export OBIT_EXEC=${OBIT}
export PATH=${PATH}:${OBIT}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${OBIT}/lib
export PYTHONPATH=${PYTHONPATH}:${OBIT}/share/obittalk/python:${OBIT}/share/python
