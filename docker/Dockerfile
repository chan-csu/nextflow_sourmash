FROM conda/miniconda3
RUN conda install -y -c conda-forge -c bioconda sourmash khmer conda 
RUN sed -i 's/tkagg/agg/g' /usr/local/lib/python*/*-packages/matplotlib/mpl-data/matplotlibrc 

RUN apt-get update --fix-missing && apt-get install -y \
    wget \
    gcc \
    g++ \
    git \
    unzip \
    make \
    curl \
    bzip2 && \
    apt-get clean

#installs scythe
RUN cd / && git clone https://github.com/vsbuffalo/scythe.git && \
    cd scythe && \
    make all && \
    cd ..

#installs sickle
RUN cd / && git clone https://github.com/najoshi/sickle.git && \
    cd sickle && \
    make && \
    cd ..
    
RUN mv /scythe/scythe /sickle/sickle /usr/local/bin/ && \
    rm -Rf /scythe /sickle

# ps is needed for nextflow monitering
RUN apt-get install -y procps