FROM nvidia/cuda:11.8.0-devel-ubuntu20.04

# Update default packages
RUN apt-get update

# Get Ubuntu packages
RUN apt-get install -y \
 build-essential \
 curl 
 
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y 

ENV RUSTUP_TOOLCHAIN=1.64.0
