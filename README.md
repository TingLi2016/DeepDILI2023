# Running Deep DILI

## Prerequisites

- Due to the limitation of file size, data.tgz file is provided by request.
- Some variant of Docker or Moby should be installed. For Windows and MacOS,
  [Docker Desktop](https://www.docker.com/products/docker-desktop/) is
  recommended. For Linux, you should probably install via your package manager.
- Git source control manager will be useful, though is not required, for
  downloading the source code repository for study, and for the example data
  which is required to run the analysis.
- It is not necessary to have Python or Anaconda installed on the system where
  the analysis will be run. Instead, controlled versions of these will be provided
  within a Docker image which have no dependence on (and will not interfere with)
  any software installed on your own system.

## Pull Image from DockerHub

Unless you build your own `deep-dili` image as described below, you will need to pull
the `deep-dili` container image from DockerHub:

```sh
docker pull scharris/deep-dili
```

## Fetch Model Data

Download the model [data](https://ncsvmgitlab.fda.gov/deep-dili-devs/deep-dili/-/raw/main/data.tgz),
and unzip it somewhere. This directory contains input files that are necessary for model training
at the top level, and its `training` subdirectory will serve as an output location for files
generated during training.

## Run via command line

First perform training to generate the `data/training/online-objs.pkl`
file needed by the command line interface.

```sh
./run-training.sh scharris/deep-dili:1.1
```

The `run-cli.sh` accepts the image name as argument and sdf data via standard
input:
```sh
cat data/4mols.sdf | ./run-cli.sh scharris/deep-dili:1.1
```

If you've built your own local image (described below), you can use that name
as the first script argument instead in any of the above commands. For example:
```sh
# use local/custom image name
cat data/4mols.sdf | ./run-cli.sh deep-dili
```

### Run Web Service

To start a web service which can be used to make predictions over a network:

```sh
# linux, using host networking
docker run --name deep-dili -d -v $(pwd)/data:/data --network=host deep-dili

# MacOS
docker run --name deep-dili -d -v $(pwd)/data:/data -p127.0.0.1:80:5000 deep-dili

# Windows using portmapping
docker run --name deep-dili -d -v $pwd/data:/data -p127.0.0.1:80:5000 deep-dili
```

The server takes several minutes to start because it has to do training prior to accepting requests.
You can follow the log output to wait for the startup message:
```
docker logs deep-dili -f
```

When the service is ready for requests, you should see messages like below:
```
Starting web application.
 * Serving Flask app 'deep-dili'
 * Debug mode: off
WARNING: This is a development server. Do not use it in a production deployment. Use a production WSGI server instead.
 * Running on all addresses (0.0.0.0)
 * Running on http://127.0.0.1:5000
 * Running on http://192.168.0.136:5000
Press CTRL+C to quit
```

Then open a browser to http://localhost:5000/ to try the app.

The predicition service can also be used from the command line as shown in the examples below.


```sh
# These examples are for the service running via portmapping (second run command above).
# If running via host networking instead (first command above), add ":5000" after localhost.

# via curl:
curl -XPOST --data-binary "@data/4mols.sdf" http://localhost/deep-dili/predict-sdf

# or httpie:
http http://localhost/deep-dili/predict-sdf @data/4mols.sdf

# via SMILES code (experimental)
http 'http://localhost/deep-dili/predict-smiles?s=CO[C@H]1\C=C\O[C@@]2(C)Oc3c(C)c(O)c4C(=C(NC(=O)\C(=C/C=C/[C@H](C)[C@H](O)[C@@H](C)[C@@H](O)[C@@H](C)[C@H](OC(=O)C)[C@@H]1C)\C)\C(=C/NN5CCN(CC5)C6CCCC6)\C(=O)c4c3C2=O)O'
```

## Container Development

The below sections are only for developers wanting to build their own container
images after making modifications to the project source code and/or container
configuration.

Start by cloning the repository:

```sh
git clone https://ncsvmgitlab.fda.gov/deep-dili-devs/deep-dili.git

cd deep-dili-svc
```

### Build Local Container Image

The container image build definition is provided in file `Dockerfile`. To build
a new container image after this file or any project source code has been
changed:

```sh
# (in project directory)
docker build . -t deep-dili
```

The name `deep-dili` here can be replaced with another name, so long as
publication and run commands elsewhere in this documentation are modified
accordingly.


### Development Workflow

```
docker build . -t deep-dili
./run.sh deep-dili
docker logs deep-dili -f
# (ctrl-c when log shows that web server has started)

http http://127.0.0.1/deep-dili/predict-sdf @data/4mols.sdf
```

### Start a Shell in the Container

If you want to look around or experiment inside the container used for analysis,
you can start the container with a shell instead of automatically running the analysis:

```sh
docker run -it --rm -v $(pwd)/data:/data --entrypoint /bin/bash deep-dili
```

This should put you at the `#` prompt of a root-shell within the container.

Type `[CTRL]-D` or `exit` to exit the container shell, which will also dispose
of the container itself (because of the `--rm` argument to our `docker run`
command above).

### Push Docker Image to DockerHub

After you've built an image, you can publish it to DockerHub as shown below.
You'll first need to register as a user on DockerHub if you haven't already.

```sh
user="<your dockerhub username here>" # or replace all "$user" occurrences below with your DockerHub username
# Login to DockerHub. Enter your DockerHub username and password when prompted.
docker login

# (replace the version number at the end with whatever version number)
docker tag deep-dili $user/deep-dili:x.y
docker push $user/deep-dili:x.y

# It's good to also provide a "latest" tag if this is latest version of the image. If a user doesn't
# provide a tag on an image name within a docker run/pull command, "latest" is assumed.
docker tag deep-dili $user/deep-dili:latest
docker push $user/deep-dili:latest
```

## Deploy to server
Create a directory for the application files on the server, such as `/opt/deep-dili`.

Download the
[data](https://ncsvmgitlab.fda.gov/deep-dili-devs/deep-dili/-/raw/main/data.tgz),
and unzip it in the application directory created above. There should now be a directory `/opt/deep-dili/data/`
with file `QSAR_year_338_pearson_0.9.csv` (and some others) and a subdirectory `training`.

Create unit file `deep-dili.service` somewhere on the server with contents like below:

```
[Unit]
Description=Deep DILI
After=docker.service
Requires=docker.service

[Service]
TimeoutStartSec=0
Restart=always
ExecStartPre=-/bin/docker stop deep-dili
ExecStartPre=-/bin/docker rm deep-dili
ExecStartPre=/bin/docker pull scharris/deep-dili:0.9
ExecStart=/bin/docker run --name deep-dili -v /opt/deep-dili/data:/data --network=host scharris/deep-dili:1.0

[Install]
WantedBy=multi-user.target
```

Then copy the Systemd unit file to the system directory:

```console
 cp deep-dili.service /etc/systemd/system/
```

and try starting the service:

```console
systemctl start deep-dili
# check status
systemctl status deep-dili
```

The service should now be running on port 5000, so try accessing it in a browser.

If everything is working, enable the service to start at boot:

```console
systemctl enable deep-dili
```

## Change Notes vs TingLi2016 / DeepDILI

- Separate the test data prediction process from the training and validation process, and use only in-memory
  structures for the test data predictions.
- Unify treatment of base classifiers to simplify the code a bit.
- Select significant features in a simpler way.
- Prefer to get molecule ids from SDF header.
- main.py
  - Allow passing in data directories and number of significant features as command line arguments.
  - A new `data` directory was added in the repository which is ready to be used with (mapped into) the Docker container.
  - Allow controlling the external validation set. If file 'external_mold2.csv' exists in the data directory, then its contents will be used as an external validation set, else a subset of file 'QSAR_year_338_pearson_0.9.csv' is used.
- mold2\_\*.py
  Allow directories in results to already exist without causing error, to make it convenient to re-run analyses without having to clean the data directory of some results directories between runs.
- mold2_DeepDILI.py
  This file had a couple of incompatibilities with Tensorflow 2, which became the default version in Anaconda3 in Jan 2019.
  - Around line 40, the way that a random seed was being set would not work for TF2.
  - Around line 157, in dili_prediction "tf.reset_default_graph()" is no longer valid in TF2.
