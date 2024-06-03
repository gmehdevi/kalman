#!/bin/bash

make speed > /dev/null

echo "Testing 10 runs for 90 minutes each."

for i in {1..10}
do

  echo -n "run : ${i} - "
  # ./imu-sensor-stream-linux -d 90 > /dev/null &
  ./imu-sensor-stream-linux -d 90 -g 0.1 > /dev/null &
  IMU_PID=$!

  # ./kalman  > /dev/null &
  ./kalman --gps 0.01 > /dev/null &
  KAL_PID=$!

  wait $IMU_PID
  IMU_RET=$?

  wait $KAL_PIDs

  if [[ IMU_RET -eq 0 ]]; then
    echo "OK"
  fi

  if [[ $i -eq 3 ]]; then
    echo "Do you want to continue? (y/n)"
    read -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      exit 1
    fi
  fi
done



