#!/usr/bin/env bash
java -server -XX:-UseParallelGC -Dcom.sun.management.jmxremote -Dcom.sun.management.jmxremote.port=9000 -Dcom.sun.management.jmxremote.authenticate=false -Dcom.sun.management.jmxremote.ssl=false -jar dist/Genend.jar
