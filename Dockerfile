# base image
FROM continuumio/miniconda3

WORKDIR /app

# create the environment
COPY environment.yml .
RUN conda env create -f environment.yml
RUN echo "conda activate compact" >> ~/.bashrc

# install build-essential, required for rust compilation of fastrbo
RUN apt-get update && apt-get install build-essential -y

# make RUN commands use the newly created environment
SHELL ["conda", "run", "-n", "compact", "/bin/bash", "-c"]

# clone and install fastrbo
RUN git clone https://github.com/joerivstrien/fastrbo.git
WORKDIR /app/fastrbo
RUN maturin develop --release

WORKDIR /app
