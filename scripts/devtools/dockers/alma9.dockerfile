FROM almalinux:9
RUN dnf update -y
RUN dnf install -y epel-release
RUN dnf group install -y "Development Tools"
RUN dnf install -y gcc-fortran gcc-c++ unzip python3-devel environment-modules cmake autoconf automake mpich-devel ncurses
ARG username=testuser
RUN useradd -m -s /bin/bash $username
USER $username
RUN git clone --recursive https://github.com/LLNL/spheral /home/$username/spheral
WORKDIR /home/$username/spheral
RUN source /home/$username/.bashrc && module load mpi && ./scripts/devtools/tpl-manager.py --spec "spheral+mpi%gcc"
RUN mv *.cmake build.cmake
RUN python3 scripts/devtools/host-config-build.py --host-config build.cmake --build-dir almabuild --install-dir /home/$username/install
WORKDIR /home/$username/spheral/almabuild/build
RUN make -j 36
RUN make -j 36 install
WORKDIR /home/$username/install
RUN source /home/$username/.bashrc && module load mpi && mpirun -n 2 ./spheral -c "import Spheral"
