# systems/lorenz.py

def system(X, param):
    """Правая часть уравнений движения Аттрактора Лоренца."""
    dx = param[0] * (X[1] - X[0])
    dy = X[0] * (param[1] - X[2]) - X[1]
    dz = X[0] * X[1] - (8.0 / 3.0) * X[2]
    return [dx, dy, dz]


def runge_kutta(points, param, step):
    """Метод Рунге-Кутты 4-го порядка на чистых списках Python."""
    k1 = system(points, param)

    points1 = []
    for i in range(len(points)):
        points1.append(points[i] + step * 0.5 * k1[i])

    k2 = system(points1, param)

    points2 = []
    for i in range(len(points)):
        points2.append(points[i] + step * 0.5 * k2[i])

    k3 = system(points2, param)

    points3 = []
    for i in range(len(points)):
        points3.append(points[i] + step * k3[i])

    k4 = system(points3, param)

    next_point = []
    for i in range(len(points)):
        next_point.append(
            points[i] + step * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6
        )

    return next_point
