ARG PYVER=3.10
FROM python:${PYVER}-slim
LABEL maintainer Christophe Trophime <christophe.trophime@lncmi.cnrs.fr>

ARG PYVER=3.10
ARG PV_VERSION_MAJOR=5.12
ARG PV_VERSION_MINOR=0
ARG PV_VERSION=${PV_VERSION_MAJOR}.${PV_VERSION_MINOR}
ARG PV_FLAVOR=osmesa-MPI
#ARG PV_FLAVOR=egl-MPI
ARG ARCH=x86_64

USER root

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/opt/paraview/bin:$PATH \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

# This Dockerfile adds a non-root user with sudo access. Use the "remoteUser"
# property in devcontainer.json to use it. On Linux, the container user's GID/UIDs
# will be updated to match your local UID/GID (when using the dockerFile property).
# See https://aka.ms/vscode-remote/containers/non-root-user for details.
ARG USERNAME=feelpp
ARG USER_UID=1001
ARG USER_GID=$USER_UID

# add python package for python dev
RUN pip install --upgrade pip && \
    pip install autopep8 black bandit flake8 mypy pycodestyle pydocstyle pylint pytest

# install Feelpp from BinTray Debian repository
RUN apt update \
    && apt install -y lsb-release wget curl sudo bash-completion openssh-client git git-lfs \
    && lsb_release -cs

# install paraview headless
ARG PV_URL=https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v${PV_VERSION_MAJOR}&type=binary&os=Linux&downloadFile=ParaView-${PV_VERSION}-${PV_FLAVOR}-Linux-Python${PYVER}-${ARCH}.tar.gz
# https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.12&type=binary&os=Linux&downloadFile=ParaView-5.12.0-egl-MPI-Linux-Python3.10-x86_64.tar.gz
RUN mkdir -p /opt/paraview && \
    cd /opt/paraview && \
    wget -qO- ${PV_URL} | tar --strip-components=1 -xzv

# Make ssh dir
# Create known_hosts
# Add bitbuckets key
RUN mkdir -p /root/.ssh/ \
    && touch /root/.ssh/known_hosts \
    && ssh-keyscan -T 60 github.com >> /root/.ssh/known_hosts

# install Paraview pre-requisites
RUN apt update \
    && apt install -y libgomp1 libgl1

# create a virtual env
# install meshlib, pandas, pint, ...

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

#set up user so that we do not run as root
RUN useradd -m -s /bin/bash -G sudo,video feelpp && \
    mkdir -p  /etc/sudoers.d/ && \
    echo "feelpp ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/feelpp

USER $USERNAME
ENV HOME /home/$USERNAME
WORKDIR /home/feelpp

