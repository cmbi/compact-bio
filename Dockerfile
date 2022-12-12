# base image
FROM continuumio/miniconda3

WORKDIR /app

# create the environment
COPY environment.yml .
RUN conda env create -f environment.yml
RUN echo "conda activate compact" >> ~/.bashrc

# make RUN commands use the newly created environment
SHELL ["conda", "run", "-n", "compact", "/bin/bash", "-c"]

# check that the environment is activated
RUN echo "Make sure rbo is installed"
RUN python -c "import rbo"


