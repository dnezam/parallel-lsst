FROM ubuntu:20.04

RUN apt-get update
RUN apt-get install -y python3.8
RUN apt-get install -y gdb
RUN apt-get install -y gcc-8
RUN apt-get install -y g++-8

RUN apt-get install -y wget
RUN wget https://github.com/bazelbuild/bazelisk/releases/download/v1.15.0/bazelisk-linux-amd64
RUN mv bazelisk-linux-amd64 /usr/local/bin/bazel
RUN chmod +x /usr/local/bin/bazel

RUN apt-get remove -y gcc-9
RUN ln -s /usr/bin/gcc-8 /usr/bin/gcc
RUN ln -s /usr/bin/g++-8 /usr/bin/g++
RUN ln -s /usr/bin/python3.8 /usr/bin/python

RUN apt-get install -y clangd-12
RUN apt-get install -y git

RUN ln -s /usr/bin/clangd-12 /usr/bin/clangd

RUN echo "HELLO!"