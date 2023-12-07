#!/bin/bash
# $1: RScriptName; $2: conda env name; $3: parameters

cd ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/temp
chmod 777 ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/$1
if [ "$2" == "KEGG" ] || [ "$2" == "global_feature" ]; then
    echo "Regular Modules"
    source activate BPCells
else
    source activate "$2"
fi

Rscript ~/apache-tomcat-9.0.68/webapps/PulmoScope/tool/$1 "${@:3}"



