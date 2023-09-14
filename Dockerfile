# docker build . -t deep-dili
# docker run --name deep-dili -d -v $(pwd)/data:/data --network=host deep-dili

# First stage: compile mold2.
FROM gcc:4.9 as compile

COPY mold2-src/*.cpp mold2-src/*.h /deep-dili/mold2-src/
RUN g++ -std=c++11 /deep-dili/mold2-src/*.cpp -o /deep-dili/mold2

# Second stage: install Python scripts and dependencies, setup app data area.
FROM python:3.10-slim
COPY --from=compile /deep-dili/mold2 /deep-dili/

ENV PIP_ROOT_USER_ACTION=ignore

RUN pip install --upgrade pip && pip install poetry==1.1.8

WORKDIR /deep-dili

COPY *.py pyproject.toml poetry.lock ./

RUN poetry install

RUN mkdir -p /deep-dili/webapp
COPY webapp/* ./webapp/

ADD data.tgz /

ENTRYPOINT ["poetry", "run", "python"]
CMD ["/deep-dili/main_cli.py", "/data/training/online-objs.pkl", "/dev/stdin"]
