#!/bin/bash

script_name="$0"
docker_img="docker.io/jupyter/datascience-notebook"

usage="
USAGE:

${script_name} runner mode src_dir [-p port] [-k key] [-b base_url] [-t tag] [-c container]

Required argument:
- runner: 'docker' | 'podman'
- mode: 'start' or 'stop'
- src_dir: Host path to container's src folder

Optional arguments:
- p:    Jupyter service port
- k:    Jupyter service custom hashed password
- b:    Base URL
- t:    ${docker_img} Docker image tag
- c:    Container name

"

# Set runner
if [ $# -eq 0 ] || { [ "$1" != "docker" ] && [ "$1" != "podman" ]; }; then
    echo "$usage"
    exit 1
fi

runner="$1"
shift

# Set running mode
if [ $# -eq 0 ] || { [ "$1" != "start" ] && [ "$1" != "stop" ]; }; then
    echo "$usage"
    exit 1
fi

mode="$1"
shift

# Set host source path
if [ $# -eq 0 ]; then
    echo "$usage"
    exit 1
fi

host_src="$1"
shift

# Set custom parameter values
port="8888"
key="argon2:\$argon2id\$v=19\$m=10240,t=10,p=8\$01wrTKfayzlVtEHrHl8lJA\$E2vIbqk/3gr4WtGJKFSLTTqqKBqSplC37puuI3gp644"
base_url="small_rna_metavir"
tag="2022-10-09"
container_name="small_rna_metavir_jupyter"

while getopts p:k:b:t:c: arg
do
    case "${arg}" in
        p) port=${OPTARG};;
        k) key=${OPTARG};;
        b) base_url=${OPTARG};;
        t) tag=${OPTARG};;
        c) container_name=${OPTARG};;
    esac
done

user=$base_url

# Stop docker jupyter service
echo "Stoping Jupyter service..."
$runner stop $container_name && $runner rm $container_name
echo "Success!"
echo "Container ${container_name} stopped..."

# Start docker jupyter service
if [ $mode == "start" ]; then

    echo "Starting Jupyter service..."
    
    $runner run -d -t --rm \
            --name $container_name \
            -p ${port}:8888 \
            --user root \
            -e NB_USER="$user" \
            -e NB_UID=1002 \
            -e CHOWN_HOME=yes \
            -w "/home/${NB_USER}" \
            -v $host_src:/home/$user/src \
        ${docker_img}:${tag} \
            start-notebook.sh \
                --NotebookApp.base_url=/$base_url \
                --NotebookApp.password=\'$key\'

    echo "Success!"
    echo "Container ${container_name} started..."
    echo "Server running at http://localhost:${port}/${base_url}"
fi
