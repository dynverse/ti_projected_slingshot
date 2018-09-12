FROM dynverse/dynwrap:bioc

LABEL version 0.1.3

RUN apt-get install -y libcgal-dev libglu1-mesa-dev

RUN R -e 'devtools::install_github("kstreet13/slingshot")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
