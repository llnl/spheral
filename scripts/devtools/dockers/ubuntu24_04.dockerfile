FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y bzip2 ca-certificates g++ gcc gfortran git gzip lsb-release patch python3 tar unzip xz-utils zstd libtool curl wget libcurl4-openssl-dev tk-dev autotools-dev build-essential python3-dev python3-pip python3-venv cmake autoconf automake libopenmpi-dev libreadline-dev
ARG username=testuser
RUN useradd -m -s /bin/bash $username
USER $username
RUN git clone --recursive https://github.com/LLNL/spheral /home/$username/spheral
WORKDIR /home/$username/spheral
RUN ./scripts/devtools/tpl-manager.py --spec "spheral+mpi%gcc"
RUN mv *.cmake build.cmake
RUN python3 scripts/devtools/host-config-build.py --host-config build.cmake --build-dir almabuild --install-dir /home/$username/install
WORKDIR /home/$username/spheral/almabuild/build
RUN make -j 36
RUN make -j 36 install
WORKDIR /home/$username/install
RUN mpirun -n 2 ./spheral -c "import Spheral"
