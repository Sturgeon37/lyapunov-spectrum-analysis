# methods/Lyapunov-computation.py
import math
from math import log, inf
import numpy as np
from ..systems.lorenz import runge_kutta


def create_eps(dim_sys, eps):
    E = []
    for i in range(dim_sys):
        stroka = []
        for j in range(dim_sys):
            if i == j:
                stroka.append(eps)
            else:
                stroka.append(0.0)
        E.append(stroka)
    return E


def create_zero_matrix(dim_sys):
    E = []
    for i in range(dim_sys):
        stroka = []
        for j in range(dim_sys):
            stroka.append(0.0)
        E.append(stroka)
    return E


def scalar(a, b):
    res = 0.0
    for i in range(len(a)):
        res += a[i] * b[i]
    return res


def norma(a):
    sum_of_squares = sum([x ** 2 for x in a])
    res = np.sqrt(sum_of_squares)
    return res


def project(a, b):
    scalar_ab = scalar(a, b)
    scalar_bb = scalar(b, b)
    if scalar_bb == 0:
        return [0.0] * len(b)
    scale = scalar_ab / scalar_bb
    return [scale * bi for bi in b]


def ort_gram_shmidt(vectors):
    orto = [vectors[0].copy()] 
    for i in range(1, len(vectors)):
        vector = vectors[i].copy()
        for j in range(i):
            proj = project(vectors[i], orto[j])
            for k in range(len(vector)):
                vector[k] = vector[k] - proj[k]
        orto.append(vector)
    return orto


def count_of_lyap_var2(point, param, amount_of_lyap, step, iterations, eps):
    """Расчет спектра показателей Ляпунова методом shadow-траекторий (Бенеттин)."""
    dim_sys = len(point)
    perturb = create_eps(dim_sys, eps)
    perturbed_point = create_zero_matrix(dim_sys)
    start_point_for_perturb = create_zero_matrix(dim_sys)
    deltas = create_zero_matrix(dim_sys)
    norms = [0.0] * dim_sys
    summa = [0.0] * amount_of_lyap
    lam = [0.0] * amount_of_lyap

    # Прогрев системы на аттрактор
    for i in range(200000):
        point = runge_kutta(point, param, step)
        if any(x > 100 for x in point):
            lam = [inf] * amount_of_lyap
            return lam

    # Основной рабочий цикл
    for i in range(iterations):
        usual_point = runge_kutta(point, param, step)
        if any(x > 100 for x in usual_point):
            lam = [inf] * amount_of_lyap
            return lam

        for g in range(amount_of_lyap):
            for gg in range(dim_sys):
                start_point_for_perturb[g][gg] = point[gg] + perturb[g][gg]

        for j in range(amount_of_lyap):
            perturbed_point[j] = runge_kutta(start_point_for_perturb[j], param, step)
            if any(x > 100 for x in perturbed_point[j]):
                lam = [inf] * amount_of_lyap
                return lam

        for f in range(amount_of_lyap):
            for ff in range(dim_sys):
                deltas[f][ff] = perturbed_point[f][ff] - usual_point[ff]

        # Ортогонализация Грамма-Шмидта
        deltas_ort = ort_gram_shmidt([deltas[f] for f in range(amount_of_lyap)])

        for k in range(len(deltas_ort)):
            norms[k] = norma(deltas_ort[k])
            if any(x > 1000 for x in norms):
                lam = [inf] * amount_of_lyap
                return lam
            summa[k] += log(norms[k] / eps)

        # Перемасштабирование возмущений на базе ОРТОГОНАЛЬНЫХ векторов (deltas_ort)
        for jj in range(amount_of_lyap):
            for jjj in range(dim_sys):
                perturb[jj][jjj] = eps * (deltas_ort[jj][jjj] / norms[jj])

        point = usual_point

    for i in range(amount_of_lyap):
        lam[i] = summa[i] / (iterations * step)

    return lam
