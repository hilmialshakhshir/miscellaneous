#!/bin/sh

./build.sh
./dcm_pipeline.sh -id /data/macsbio_DCM/input/WES -od /data/macsbio_DCM/output -rd /data/macsbio_DCM/reference -td /data/macsbio_DCM/tempdir
