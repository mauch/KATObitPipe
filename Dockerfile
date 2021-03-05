# KATObitPipe Dockerfile
# 
# MAINTAINER: Tom Mauch (tmauch@ska.ac.za)

FROM ubuntu:18.04

# Suppress debconf warnings
ENV DEBIAN_FRONTEND noninteractive

ENV PACKAGES \
    build-essential \
    ca-certificates \
    curl \
    libgfortran3 \
    gosu \
    libfreetype6 \
    libxml2 \ 
    pkg-config \
    python3-dev \
    rsync \
    virtualenv

# Update, upgrade and install packages
RUN apt-get -y update && \
    apt-get --no-install-recommends -y install $PACKAGES

# Set up python 3 virtualenv
RUN virtualenv -p /usr/bin/python3 /ve3 && \
	. /ve3/bin/activate && \
	pip install pip --upgrade
ENV PATH="/ve3/bin:${PATH}"
ENV VIRTUAL_ENV="/ve3"

# Retrieve Obit and set up environment
ENV OBIT_BASE_PATH="/Obit"
ENV OBIT_REVISION="630"
ENV OBIT_TARBALL="obit-distro-1.1.${OBIT_REVISION}.tar.gz"

RUN mkdir -p $OBIT_BASE_PATH && cd $OBIT_BASE_PATH && \
    curl https://svn.cv.nrao.edu/obit/linux_distro/${OBIT_TARBALL} | tar xzf -

ENV OBIT_ROOT="${OBIT_BASE_PATH}/obit-distro-1.1.${OBIT_REVISION}"
ENV OBIT="${OBIT_ROOT}/lib/obit"
ENV OBIT_EXEC="${OBIT}"
ENV PATH="${OBIT}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${OBIT_ROOT}/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH="${OBIT_ROOT}/share/obittalk/python:${OBIT_ROOT}/share/python:${PYTHONPATH}"

# Get rid of Obit's gfortran because AIPS doesn't like it
RUN rm ${OBIT_ROOT}/lib/libgfortran*
RUN ln -s /usr/lib/x86_64-linux-gnu/libgfortran.so.3 ${OBIT_ROOT}/lib

# Install the pipeline
COPY . /KATObitPipe
WORKDIR /KATObitPipe

# Copy metadata to Obit share
RUN cp ./FITS/* ${OBIT_ROOT}/share/obit/data
RUN cd ${OBIT_ROOT}/share/obit/data && gunzip *.fits.gz
RUN chmod -R 777 ${OBIT_ROOT}/share/obit/data


RUN pip install -r requirements.txt
RUN pip install .

RUN mkdir /scratch
VOLUME /scratch

WORKDIR /scratch

RUN chmod +x /KATObitPipe/entrypoint.sh
ENTRYPOINT ["/KATObitPipe/entrypoint.sh"]
CMD ["/bin/bash"]
