#!/bin/bash
imgname=mosaic

echo "[d2s] Starting a local docker registry ..."
docker run -d -p 5000:5000 --restart=always --name registry registry:2
echo -e "\n"


echo "[d2s] Build Docker image from Dockerfile ..."
docker build -t ${imgname}:dev --network host .
echo -e "\n"


echo "[d2s] Tag newly build image to local registry ..."
docker tag ${imgname}:dev localhost:5000/${imgname}
echo -e "\n"


echo "[d2s] Push newly build image to local registry ..."
docker push localhost:5000/${imgname}
echo -e "\n"


echo "[d2s] Build singularity image file from Docker image ..."
SINGULARITY_NOHTTPS=1 singularity build ${imgname}.simg docker://localhost:5000/${imgname}
echo -e "\n"


echo "[d2s] Stop local registry ..."
docker container stop registry
echo -e "\n"


echo "[d2s] Remove registry container ..."
docker container ls -aqf name=registry | xargs docker rm

