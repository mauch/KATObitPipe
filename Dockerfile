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
    libfreetype6 \
    libxml2 \ 
    pkg-config \
    python3-dev \
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
ENV OBIT_TARBALL="Obit.AVX-1.1.${OBIT_REVISION}.tar.gz"

RUN mkdir -p $OBIT_BASE_PATH && cd $OBIT_BASE_PATH && \
    curl https://svn.cv.nrao.edu/obit/linux_distro/${OBIT_TARBALL} | tar xzf -

ENV OBIT="${OBIT_BASE_PATH}/obit-distro-1.1.${OBIT_REVISION}"
ENV PATH="${OBIT}/bin:${PATH}"
ENV LD_LIBRARY_PATH="${OBIT}/lib:${LD_LIBRARY_PATH}"
ENV PYTHONPATH="$OBIT/share/obittalk/python:$OBIT/share/python:${PYTHONPATH}"

# Install the pipeline
COPY . /KATObitPipe
WORKDIR /KATObitPipe
RUN pip install -r requirements.txt
RUN pip install .




