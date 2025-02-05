FROM index.docker.io/python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive

COPY ../dist/asap_o_matic-2.3.0-cp310-abi3-linux_x86_64.whl /opt/
RUN pip install /opt/asap_o_matic-2.3.0-cp310-abi3-linux_x86_64.whl

ENTRYPOINT [ "asap-o-matic" ]