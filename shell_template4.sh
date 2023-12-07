#!/bin/bash
# $1:RScriptName; $2:parameters
source activate giotto;
cd ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/temp;
chmod 777 ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/$1;
Rscript ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/$1 ${@:2}
