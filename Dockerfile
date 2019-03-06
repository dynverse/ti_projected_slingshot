FROM dynverse/dynwrapr:v0.1.0

ARG GITHUB_PAT

RUN apt-get update && apt-get install -y libcgal-dev libglu1-mesa-dev libgsl-dev

RUN R -e 'devtools::install_github("kstreet13/slingshot")'

COPY definition.yml example.R run.R /code/

ENTRYPOINT ["/code/run.R"]
