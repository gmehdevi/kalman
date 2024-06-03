#! /bin/bash

echo "Testing filter mean speed."

make speed > /dev/null

./imu-sensor-stream-linux --duration 10  -g 0.1 --filterspeed &
IMU_PID=$!

./kalman -a 0.01 --gps 0.01 > /dev/null &
KAL_PID=$!

wait $IMU_PID
IMU_RET=$?

wait $KAL_PID

echo "Expected to be under 0.0001"

echo "Press any key to conitnue or Ctrl+C to exit."
read

echo "Graphing first 10 minutes of random seed."
./graph_seed.sh $RANDOM

echo "Press any key to conitnue or Ctrl+C to exit."
read

echo "Testing 5 runs for 60 minutes each with increasing amounts of noise."
for i in {1..5}; do
  IMU_RET=1
  echo -n "run : noise: 1.$i - "
  while [[ $IMU_RET -ne 0 ]]; do
    ./imu-sensor-stream-linux --duration 60  -n 1.$i 2>/dev/null >/dev/null &
    IMU_PID=$!

    ./kalman -a 0.01 --gps  2500 -n 1.$i > /dev/null &
    KAL_PID=$!

    wait $IMU_PID
    IMU_RET=$?

    wait $KAL_PID
  done
  echo "OK"
done
