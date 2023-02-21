#!/bin/bash

# Simulation executable
EXEC=main.x

# Path to simulation executable
EXEC_PATH=/home/spaceghoul/Desktop/Spring2023/AE6042/MHD_MUSCL_n_WENO

# Path to output file
OUTPUT=/home/spaceghoul/Desktop/Spring2023/AE6042/MHD_MUSCL_n_WENO/stdout

# Set time interval (s) between system resource checks
t_interval=10

# Start simulation in background
cd $EXEC_PATH
nohup ./$EXEC <<< "1 2 0" > $OUTPUT &

# Get PID of simulation
PID=$!

# Check system resources
while kill -0 $PID > /dev/null 2>&1; do
	echo "Simulation running (PID $PID)"
	# Check more info about simulation
	ps aux | grep $EXEC
	# Check system resources
	echo "------------------------------"
	echo "CPU usage: $(top -b -n1 | grep "Cpu(s)" | awk '{print $2+$4 "%"}')"
	echo "RAM usage: $(free -m | awk 'NR==2{printf "%sMB/%sMB (%.2f%%)\n", $3,$2,$3*100/$2 }')"
	echo "Disk usage: $(df -h | awk '$NF=="/"{printf "%dGB/%dGB (%s)\n", $3,$2,$5}')"
	echo "------------------------------"
	sleep $t_interval
done
