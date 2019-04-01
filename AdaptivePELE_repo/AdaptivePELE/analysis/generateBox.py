import os
import argparse
import numpy as np
from AdaptivePELE.utilities import utilities


def parse_arguments():
    """
        Create command-line interface
    """
    desc = "Create a visualization for a cylindrical box"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-t", "--top", type=float, nargs="*", required=True, help="Coordinates for the bottom base")
    parser.add_argument("-b", "--bottom", type=float, nargs="*", required=True, help="Coordinates for the top base")
    parser.add_argument("--radius", type=float, default=10, required=True, help="Radius of the cylinder base")
    parser.add_argument("--n_points", type=int, default=10, required=True, help="Number of points to write, the larger it is the finer the discretization")
    parser.add_argument("--title", type=str, default="box.pdb", help="Name to store the output")
    parser.add_argument("--output", type=str, default="", help="Path to store the output")
    args = parser.parse_args()
    return args.top, args.bottom, args.radius, args.title, args.output, args.n_points


def main(top, bot, radius, title="box.pdb", output="", n_points=10):
    if output:
        if not os.path.exists(output):
            os.makedirs(output)
        title = os.path.join(output, title)
    top = np.array(top)
    bot = np.array(bot)
    l = np.linalg.norm(top-bot)
    alpha = top-bot
    alpha /= np.linalg.norm(alpha)
    xp, yp = np.random.rand(2)
    alpha_rot = np.array([xp, yp, -(alpha[0]*xp+alpha[1]*yp)/alpha[2]])
    alpha_rot /= np.linalg.norm(alpha_rot)
    alpha_rot *= radius
    angle = 2*np.pi/n_points
    rot_mat = utilities.generateRotationMatrixAroundAxis(alpha, angle)
    rotated_vectors = []
    new_alpha = alpha_rot
    for i in range(n_points):
        new_alpha = np.dot(rot_mat, new_alpha)
        assert np.abs(np.linalg.norm(new_alpha) - radius) < 1e-4, (np.linalg.norm(new_alpha), radius)
        rotated_vectors.append(new_alpha)

    coords = [top, bot, bot+alpha_rot, bot-alpha_rot, top+alpha_rot, top-alpha_rot]
    points = [bot+i*alpha for i in np.linspace(0, l, num=n_points)]
    for point in points:
        coords += [point+alpha_rot, point-alpha_rot]+[point+alpha_2 for alpha_2 in rotated_vectors]

    coords = np.array(coords)
    utilities.write_PDB_clusters(coords, title)


if __name__ == "__main__":
    t, b, rad, file_name, out, n = parse_arguments()
    main(t, b, rad, file_name, out, n)
