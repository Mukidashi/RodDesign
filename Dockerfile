FROM --platform=linux/x86-64 ubuntu

ENV DEBIAN_FRONTEND=noninteractive
ENV NVIDIA_VISIBLE_DEVICES ${NVIDIA_VISIBLE_DEVICES:-all}
ENV NVIDIA_DRIVER_CAPABILITIES ${NVIDIA_DRIVER_CAPABILITIES:+$NVIDIA_DRIVER_CAPABILITIES,}graphics

RUN echo "Asia/Tokyo" > /etc/timezone

RUN apt update -y && apt upgrade -y && apt install -y build-essential cmake
RUN apt install -y libeigen3-dev  libgl1-mesa-dev libglew-dev freeglut3 freeglut3-dev \
        libglu1-mesa-dev mesa-common-dev qtbase5-dev qttools5-dev-tools qt5-default


RUN mkdir /work
COPY ./src /work
