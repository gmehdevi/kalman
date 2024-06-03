seed=$@

make speed > /dev/null
make graph > /dev/null

./imu-sensor-stream-linux --seed $seed --duration 10 -g 0.1 --debug > /dev/null &
IMU_PID=$!

./kalman -p > truth.txt
KAL_PID=$!

wait $IMU_PID
wait $KAL_PID


./imu-sensor-stream-linux --seed $seed -g 0.1 --duration 10  > /dev/null &
IMU_PID=$!

./kalman -p --gps 0.01 > predict.txt
KAL_PID=$!

wait $IMU_PID
wait $KAL_PID


./graph predict.txt truth.txt variances.txt

