#!/bin/bash
# $1:userDirectoryName; $2:RScriptName; $3:parameters
source activate BPCells;
cd ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/temp;
mkdir -p $1;
chmod 777 ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/$2;
Rscript ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/$2 ${@:3} $1
