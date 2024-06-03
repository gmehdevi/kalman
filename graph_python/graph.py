#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.style as mplstyle
mplstyle.use('fast')

def parse_trajectory_file(file_path, n):
    with open(file_path, 'r') as file:
        if n <= 0:
            lines = file.readlines()
        else:
            lines = [file.readline() for _ in range(n)]

    trajectory = []
    lines = [line for line in lines if re.match(r'x : -?\d+(\.\d+)?, y : -?\d+(\.\d+)?, z : -?\d+(\.\d+)?', line) is not None]
    for line in lines:
        line = re.sub(r'[xyz:,\n]', ' ', line)
        x, y, z = line.split()
        trajectory.append((float(x), float(y), float(z)))

    return np.array(trajectory)

def plot_3d_trajectories(estimated_trajectory, real_trajectory, variances, percent = 10):
    print("Estimated trajectory length: ", len(estimated_trajectory))
    print("Real trajectory length: ", len(real_trajectory))

    if variances is not None and len(variances) != len(estimated_trajectory):
        print("Variance and trajectory lengths do not match.")
        return

    indices = np.arange(0, len(real_trajectory), int(100 / percent))
    indices_est = indices[indices < len(estimated_trajectory)]

    real_trajectory = real_trajectory[indices]
    estimated_trajectory = estimated_trajectory[indices_est]

    x_est, y_est, z_est = zip(*estimated_trajectory)
    x_real, y_real, z_real = zip(*real_trajectory)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(x_est, y_est, z_est, label="Estimated Trajectory", color='blue')
    ax.plot(x_real, y_real, z_real, label="Real Trajectory", color='green')

    if variances is not None:
        variances = variances[indices_est]
        num_points = 5

        for (x, y, z), variance in zip(estimated_trajectory, variances):
            radius = np.linalg.norm(np.sqrt(variance))

            phi = np.random.uniform(0, 2 * np.pi, num_points)
            theta = np.random.uniform(0, np.pi, num_points)

            gx = x + radius * np.sin(theta) * np.cos(phi)
            gy = y + radius * np.sin(theta) * np.sin(phi)
            gz = z + radius * np.cos(theta)

            ax.scatter(gx, gy, gz, color='b', alpha=0.05)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title("Trajectory comparison")

    plt.legend()
    plt.show()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Visualize estimated and real trajectories in 3D.")
    parser.add_argument("estimated_file", help="Path to the estimated trajectory file", default=None)
    parser.add_argument("real_file", help="Path to the real trajectory file", default=None)
    parser.add_argument("variance_file", help="Path to the variance file", nargs='?', default=None)
    parser.add_argument("-d", "--delta", help="Substracts second trajectory from the first one to get the real trajectory")
    parser.add_argument("-n", "--number", help="The number of points to plot over the trajectory",  default=0, type=int)
    parser.add_argument("-p", "--percent", help="The percentage of selected points to plot",  default=10, type=float)
    parser.add_argument("-s", "--start", help="The start index of the trajectory",  default=0, type=int)
    parser.add_argument("-e", "--end", help="The end index of the trajectory",  default=-1, type=int)
    args = parser.parse_args()

    variances = None
    estimated_trajectory = None
    real_trajectory = None

    estimated_trajectory = parse_trajectory_file(args.estimated_file, args.number)
    real_trajectory = parse_trajectory_file(args.real_file, args.number)

    if args.delta:
        real_trajectory = estimated_trajectory - real_trajectory[:len(estimated_trajectory)]

    if args.variance_file != None:
        variances = parse_trajectory_file(args.variance_file, args.number)[:len(estimated_trajectory)]
        variances = variances[args.start:args.end]


    estimated_trajectory = estimated_trajectory[args.start:args.end]
    real_trajectory = real_trajectory[args.start:args.end]

    plot_3d_trajectories(estimated_trajectory, real_trajectory, variances, args.percent)