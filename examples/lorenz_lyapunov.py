# examples/lorenz_lyapunov.py
from methods.benettin import count_of_lyap_var2

def main():
    print("Запуск алгоритма расчета показателей Ляпунова для Аттрактора Лоренца...")
    
    # Входные параметры
    init_point = [1.0, 1.0, 1.0]
    lorenz_params = [10.0, 28.0, 8.0/3.0]  # sigma, rho, beta
    
    # Конфигурация спектра
    amount_of_lyap = 3
    dt = 0.005
    iterations = 2000  # Можно увеличить для большей точности
    eps = 1e-8

    print(f"Вычисление запущено. Количество шагов: {iterations}...")
    spectrum = count_of_lyap_var2(
        point=init_point,
        param=lorenz_params,
        amount_of_lyap=amount_of_lyap,
        step=dt,
        iterations=iterations,
        eps=eps
    )
    
    print("\nРассчитанный спектр показателей Ляпунова:")
    for i, val in enumerate(spectrum):
        print(f"λ_{i+1} = {val:.5f}")


if __name__ == "__main__":
    main()
