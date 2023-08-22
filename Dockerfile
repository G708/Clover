FROM python:3.8
USER root

WORKDIR /opt
COPY requirements.txt ./

RUN python -m pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install --no-cache-dir -r requirements.txt
