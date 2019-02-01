FROM dynverse/dynwrap:bioc

RUN apt-get update && apt-get install -y libcgal-dev libglu1-mesa-dev libgsl-dev

RUN R -e 'devtools::install_github("kstreet13/slingshot")'

LABEL version 0.1.7.2

ADD . /code

ENTRYPOINT Rscript /code/run.R
