#include <glut.h>
#include <SDL.h>
#undef main
#include <SDL_opengl.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <cstdio>

static SDL_Window* window;

// Параметры системы Рабиновича-Фабриканта
/*double alpha = 0.05;
double gamma = 0.1;
double state[3] = { 0.1, -0.1, 0.1 };*/

double alpha = 1.1;
double gamma = 0.87;
double state[3] = { -1.0, -0.0, 0.5 };

int screen_size[2] = { 800, 600 };
double min_coords[3] = { 1e9, 1e9, 1e9 };
double max_coords[3] = { -1e9, -1e9, -1e9 };
std::vector<std::vector<double>> trajectory;

// Система Рабиновича-Фабриканта
void rabinovich_fabrikant(double* x, double* fx, void* params) {
    if (isnan(x[0]) || isnan(x[1]) || isnan(x[2])) {
        fprintf(stderr, "Ошибка: Некорректные входные значения\n");
        return;
    }
    fx[0] = x[1] * (x[2] - 1 + x[0] * x[0]) + gamma * x[0];
    fx[1] = x[0] * (3 * x[2] + 1 - x[0] * x[0]) + gamma * x[1];
    fx[2] = -2 * x[2] * (alpha + x[0] * x[1]);
}

// Метод Эйлера
void euler_step(int n, double* x0, double* xh, double h,
    void (*f)(double* x, double* fx, void* params), void* params) {
    if (!x0 || !xh) {
        fprintf(stderr, "Ошибка: Некорректные входные параметры\n");
        return;
    }
    double* k1 = (double*)malloc(sizeof(double) * n);
    if (!k1) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        free(k1);
        return;
    }
    f(x0, k1, params);
    for (int i = 0; i < n; i++) {
        xh[i] = x0[i] + h * k1[i];
    }
    free(k1);
}

// Метод предиктор-корректор Адамса-Башфорта-Мултона
void predictor_corrector_step(int n, double* x0, double* xh, double h,
    void (*f)(double* x, double* fx, void* params), void* params,
    double* x_prev, double* x_prev2) {
    if (!x0 || !xh || !x_prev || !x_prev2) {
        fprintf(stderr, "Ошибка: Некорректные входные параметры\n");
        return;
    }
    double* f0 = (double*)malloc(sizeof(double) * n);
    double* f_prev = (double*)malloc(sizeof(double) * n);
    double* f_prev2 = (double*)malloc(sizeof(double) * n);
    double* x_predictor = (double*)malloc(sizeof(double) * n);
    double* f_predictor = (double*)malloc(sizeof(double) * n);
    if (!f0 || !f_prev || !f_prev2 || !x_predictor || !f_predictor) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        free(f0);
        free(f_prev);
        free(f_prev2);
        free(x_predictor);
        free(f_predictor);
        return;
    }

    f(x0, f0, params);
    f(x_prev, f_prev, params);
    f(x_prev2, f_prev2, params);

    for (int i = 0; i < n; ++i) {
        x_predictor[i] = x0[i] + h * (1.5 * f0[i] - 0.5 * f_prev[i]);
    }

    f(x_predictor, f_predictor, params);

    for (int i = 0; i < n; ++i) {
        xh[i] = x0[i] + h * 0.5 * (f0[i] + f_predictor[i]);
    }
    free(f0);
    free(f_prev);
    free(f_prev2);
    free(x_predictor);
    free(f_predictor);
}

// Метод Дормана-Принса 8-го порядка
void dormand_prince_step(int n, double* x0, double* xh, double h,
    void (*f)(double* x, double* fx, void* params), void* params) {
    if (!x0 || !xh) {
        fprintf(stderr, "Ошибка: Некорректные входные параметры\n");
        return;
    }
    double* k1 = (double*)malloc(sizeof(double) * n);
    double* k2 = (double*)malloc(sizeof(double) * n);
    double* k3 = (double*)malloc(sizeof(double) * n);
    double* k4 = (double*)malloc(sizeof(double) * n);
    double* k5 = (double*)malloc(sizeof(double) * n);
    double* k6 = (double*)malloc(sizeof(double) * n);
    double* k7 = (double*)malloc(sizeof(double) * n);
    double* tmp = (double*)malloc(sizeof(double) * n);
    if (!k1 || !k2 || !k3 || !k4 || !k5 || !k6 || !k7 || !tmp) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        free(k1);
        free(k2);
        free(k3);
        free(k4);
        free(k5);
        free(k6);
        free(k7);
        free(tmp);
        return;
    }
    f(x0, k1, params);
    for (int i = 0; i < n; i++)
        tmp[i] = x0[i] + h * 0.2 * k1[i];
    f(tmp, k2, params);
    for (int i = 0; i < n; i++)
        tmp[i] = x0[i] + h * (0.075 * k1[i] + 0.225 * k2[i]);
    f(tmp, k3, params);

    for (int i = 0; i < n; i++)
        tmp[i] = x0[i] + h * (0.3 * k1[i] - 0.9 * k2[i] + 1.2 * k3[i]);
    f(tmp, k4, params);

    for (int i = 0; i < n; i++)
        tmp[i] = x0[i] + h * (-11.0 / 54.0 * k1[i] + 2.5 * k2[i]
        - 70.0 / 27.0 * k3[i] + 35.0 / 27.0 * k4[i]);
    f(tmp, k5, params);

    for (int i = 0; i < n; i++)
        tmp[i] = x0[i] + h * (1631.0 / 55296.0 * k1[i] + 175.0 / 512.0 * k2[i]
        + 575.0 / 13824.0 * k3[i] + 44275.0 / 110592.0 * k4[i]
            + 253.0 / 4096.0 * k5[i]);
    f(tmp, k6, params);

    for (int i = 0; i < n; i++) {
        xh[i] = x0[i] + h * (37.0 / 378.0 * k1[i] + 0.0 * k2[i]
            + 250.0 / 621.0 * k3[i] + 125.0 / 594.0 * k4[i]
                + 0.0 * k5[i] + 512.0 / 1771.0 * k6[i]);
    }
    f(tmp, k7, params);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
    free(k7);
    free(tmp);
}

// Метод Рунге-Кутты 4-го порядка
void runge_kutta_step(int n, double* x0, double* xh, double h,
    void (*f)(double* x, double* fx, void* params), void* params) {
    if (!x0 || !xh) {
        fprintf(stderr, "Ошибка: Некорректные входные параметры\n");
        return;
    }
    double* k1 = (double*)malloc(sizeof(double) * n);
    double* k2 = (double*)malloc(sizeof(double) * n);
    double* k3 = (double*)malloc(sizeof(double) * n);
    double* k4 = (double*)malloc(sizeof(double) * n);
    double* tmp = (double*)malloc(sizeof(double) * n);
    if (!k1 || !k2 || !k3 || !k4 || !tmp) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        free(k1);
        free(k2);
        free(k3);
        free(k4);
        free(tmp);
        return;
    }
    f(x0, k1, params);
    for (int i = 0; i < n; i++)
        tmp[i] = x0[i] + h * k1[i] / 2.0;
    f(tmp, k2, params);
    for (int i = 0; i < n; i++)
        tmp[i] = x0[i] + h * k2[i] / 2.0;
    f(tmp, k3, params);
    for (int i = 0; i < n; i++)
        tmp[i] = x0[i] + h * k3[i];
    f(tmp, k4, params);
    for (int i = 0; i < n; i++) {
        xh[i] = x0[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(tmp);
}

// Метод Рунге-Кутты 4-ого порядка с адаптивным шагом
void adaptive_runge_kutta_step(int n, double* x0, double* xh, double& h,
    void (*f)(double* x, double* fx, void* params), void* params, double tolerance) {
    if (!x0 || !xh) {
        fprintf(stderr, "Ошибка: Некорректные входные параметры\n");
        return;
    }
    double* k1 = (double*)malloc(sizeof(double) * n);
    double* k2 = (double*)malloc(sizeof(double) * n);
    double* k3 = (double*)malloc(sizeof(double) * n);
    double* k4 = (double*)malloc(sizeof(double) * n);
    double* tmp = (double*)malloc(sizeof(double) * n);
    double* x_half = (double*)malloc(sizeof(double) * n);

    if (!k1 || !k2 || !k3 || !k4 || !tmp || !x_half) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        free(k1); free(k2); free(k3); free(k4); free(tmp); free(x_half);
        return;
    }
    f(x0, k1, params);
    for (int i = 0; i < n; i++) tmp[i] = x0[i] + h * k1[i] / 2.0;
    f(tmp, k2, params);
    for (int i = 0; i < n; i++) tmp[i] = x0[i] + h * k2[i] / 2.0;
    f(tmp, k3, params);
    for (int i = 0; i < n; i++) tmp[i] = x0[i] + h * k3[i];
    f(tmp, k4, params);
    for (int i = 0; i < n; i++) {
        xh[i] = x0[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }

    double half_h = h / 2.0;
    f(x0, k1, params);
    for (int i = 0; i < n; i++) tmp[i] = x0[i] + half_h * k1[i] / 2.0;
    f(tmp, k2, params);
    for (int i = 0; i < n; i++) tmp[i] = x0[i] + half_h * k2[i] / 2.0;
    f(tmp, k3, params);
    for (int i = 0; i < n; i++) tmp[i] = x0[i] + half_h * k3[i];
    f(tmp, k4, params);
    for (int i = 0; i < n; i++) {
        x_half[i] = x0[i] + half_h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    f(x_half, k1, params);
    for (int i = 0; i < n; i++) tmp[i] = x_half[i] + half_h * k1[i] / 2.0;
    f(tmp, k2, params);
    for (int i = 0; i < n; i++) tmp[i] = x_half[i] + half_h * k2[i] / 2.0;
    f(tmp, k3, params);
    for (int i = 0; i < n; i++) tmp[i] = x_half[i] + half_h * k3[i];
    f(tmp, k4, params);
    for (int i = 0; i < n; i++) {
        x_half[i] += half_h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }

    double error = 0.0;
    for (int i = 0; i < n; i++) {
        error = fmax(error, fabs(xh[i] - x_half[i]));
    }
    if (error > tolerance) {
        h *= 0.9 * pow(tolerance / error, 0.25);
    }
    else {
        h *= 0.9 * pow(tolerance / error, 0.2);
    }
    memcpy(xh, x_half, sizeof(double) * n);
    free(k1); free(k2); free(k3); free(k4); free(tmp); free(x_half);
}

// Метод Адамса-Мултона 2-го порядка (неявный метод)
void adams_moulton_step(int n, double* x0, double* xh, double h,
    void (*f)(double* x, double* fx, void* params), void* params,
    double* x_prev, double* x_prev2) {
    if (!x0 || !xh || !x_prev || !x_prev2) {
        fprintf(stderr, "Ошибка: Некорректные входные параметры\n");
        return;
    }
    double* f0 = (double*)malloc(sizeof(double) * n);
    double* f_prev = (double*)malloc(sizeof(double) * n);
    double* f_prev2 = (double*)malloc(sizeof(double) * n);
    double* x_predictor = (double*)malloc(sizeof(double) * n);
    double* f_predictor = (double*)malloc(sizeof(double) * n);
    double* x_corrector = (double*)malloc(sizeof(double) * n);
    double* f_corrector = (double*)malloc(sizeof(double) * n);
    if (!f0 || !f_prev || !f_prev2 || !x_predictor || !f_predictor || !x_corrector || !f_corrector) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        free(f0); free(f_prev); free(f_prev2); free(x_predictor);
        free(f_predictor); free(x_corrector); free(f_corrector);
        return;
    }
    f(x0, f0, params);
    f(x_prev, f_prev, params);
    f(x_prev2, f_prev2, params);
    for (int i = 0; i < n; ++i) {
        x_predictor[i] = x0[i] + h * (1.5 * f0[i] - 0.5 * f_prev[i]);
    }
    double tolerance = 1e-6;
    int max_iterations = 100;
    double max_delta;
    for (int i = 0; i < n; ++i) {
        x_corrector[i] = x_predictor[i];
    }
    for (int iter = 0; iter < max_iterations; ++iter) {
        f(x_corrector, f_corrector, params);
        for (int i = 0; i < n; ++i) {
            xh[i] = x0[i] + h * 0.5 * (f0[i] + f_corrector[i]);
        }
        // Проверка сходимости
        max_delta = 0.0;
        for (int i = 0; i < n; ++i) {
            double delta = fabs(xh[i] - x_corrector[i]);
            if (delta > max_delta) {
                max_delta = delta;
            }
        }
        if (max_delta < tolerance) {
            break;
        }
        for (int i = 0; i < n; ++i) {
            x_corrector[i] = xh[i];
        }
    }
    free(f0); free(f_prev); free(f_prev2); free(x_predictor);
    free(f_predictor); free(x_corrector); free(f_corrector);
}


// Вычисление Якобиана
void numerical_jacobian(double* x, double** J,
    void (*f)(double*, double*, void*), void* params, int n, double h_numerical) {
    double* fx = (double*)malloc(sizeof(double) * n);
    double* x_perturbed = (double*)malloc(sizeof(double) * n);
    double* fx_perturbed = (double*)malloc(sizeof(double) * n);
    if (!fx || !x_perturbed || !fx_perturbed) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        exit(1);
    }
    f(x, fx, params);
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
            x_perturbed[k] = x[k];
        }
        x_perturbed[j] += h_numerical;
        f(x_perturbed, fx_perturbed, params);
        for (int i = 0; i < n; ++i) {
            J[i][j] = (fx_perturbed[i] - fx[i]) / h_numerical;
        }
    }
    free(fx);
    free(x_perturbed);
    free(fx_perturbed);
}

// Функция решения системы линейных уравнений методом Гаусса
void gauss_elimination(double** A, double* b, double* x, int n) {
    double** A_copy = (double**)malloc(sizeof(double*) * n);
    double* b_copy = (double*)malloc(sizeof(double) * n);
    if (!A_copy || !b_copy) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        exit(1);
    }
    for (int i = 0; i < n; ++i) {
        A_copy[i] = (double*)malloc(sizeof(double) * n);
        if (!A_copy[i]) {
            fprintf(stderr, "Ошибка: Не удалось выделить память\n");
            exit(1);
        }
        for (int j = 0; j < n; ++j) {
            A_copy[i][j] = A[i][j];
        }
        b_copy[i] = b[i];
    }
    for (int i = 0; i < n; ++i) {
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(A_copy[k][i]) > fabs(A_copy[max_row][i])) {
                max_row = k;
            }
        }
        if (max_row != i) {
            for (int k = 0; k < n; ++k) {
                double temp = A_copy[i][k];
                A_copy[i][k] = A_copy[max_row][k];
                A_copy[max_row][k] = temp;
            }
            double temp = b_copy[i];
            b_copy[i] = b_copy[max_row];
            b_copy[max_row] = temp;
        }
        for (int k = i + 1; k < n; ++k) {
            double factor = A_copy[k][i] / A_copy[i][i];
            for (int j = i; j < n; ++j) {
                A_copy[k][j] -= factor * A_copy[i][j];
            }
            b_copy[k] -= factor * b_copy[i];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b_copy[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A_copy[i][j] * x[j];
        }
        x[i] /= A_copy[i][i];
    }
    for (int i = 0; i < n; ++i) {
        free(A_copy[i]);
    }
    free(A_copy);
    free(b_copy);
}

// Функция решения уравнения
void solve_newton(double* x_next, double* x_guess, double h,
    void (*f)(double*, double*, void*), void* params, int n,
    double tolerance, int max_iterations, double h_numerical) {
    double* fx = (double*)malloc(sizeof(double) * n);
    double** J = (double**)malloc(sizeof(double*) * n);
    double* delta_x = (double*)malloc(sizeof(double) * n);
    double* b = (double*)malloc(sizeof(double) * n);
    if (!fx || !J || !delta_x || !b) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        exit(1);
    }
    for (int i = 0; i < n; ++i) {
        J[i] = (double*)malloc(sizeof(double) * n);
        if (!J[i]) {
            fprintf(stderr, "Ошибка: Не удалось выделить память\n");
            exit(1);
        }
    }
    for (int iter = 0; iter < max_iterations; ++iter) {
        f(x_guess, fx, params);
        for (int i = 0; i < n; ++i) {
            b[i] = x_guess[i] - x_guess[i] - h * fx[i]; // F(x_guess) = x_guess - x0 - h * f(x_guess), x0 = x_guess
        }
        numerical_jacobian(x_guess, J, f, params, n, h_numerical);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                J[i][j] = (i == j ? 1.0 : 0.0) - h * J[i][j];
            }
        }
        // J * delta_x = -F(x_guess)
        gauss_elimination(J, b, delta_x, n);
        for (int i = 0; i < n; ++i) {
            x_next[i] = x_guess[i] - delta_x[i];
        }
        // Проверка на сходимость
        double max_delta = 0.0;
        for (int i = 0; i < n; ++i) {
            double delta = fabs(x_next[i] - x_guess[i]);
            if (delta > max_delta) {
                max_delta = delta;
            }
        }
        if (max_delta < tolerance) {
            break;
        }
        for (int i = 0; i < n; ++i) {
            x_guess[i] = x_next[i];
        }
    }
    free(fx);
    for (int i = 0; i < n; ++i) {
        free(J[i]);
    }
    free(J);
    free(delta_x);
    free(b);
}

// Неявный метод Эйлера
void implicit_euler_step(double* x0, double* xh, double h,
    void (*f)(double*, double*, void*), void* params, int n,
    double tolerance, int max_iterations, double h_numerical) {
    double* x_guess = (double*)malloc(sizeof(double) * n);
    if (!x0 || !xh || !x_guess) {
        fprintf(stderr, "Ошибка: Не удалось выделить память\n");
        exit(1);
    }
    for (int i = 0; i < n; ++i) {
        x_guess[i] = x0[i];
    }
    solve_newton(xh, x_guess, h, f, params, n, tolerance, max_iterations, h_numerical);
    free(x_guess);
}


// Функции для вызовов методов отдельно
// Метод Рунге-Кутты 4
void runge_kutta_callback(double dt, double* state) {
    double newstate[3];
    runge_kutta_step(3, state, newstate, dt, rabinovich_fabrikant, NULL);
    memcpy(state, newstate, sizeof(newstate));
}

// Метод Дормана-Принса 8
void dormand_prince_callback(double dt, double* state) {
    double newstate[3];
    dormand_prince_step(3, state, newstate, dt, rabinovich_fabrikant, NULL);
    memcpy(state, newstate, sizeof(newstate));
}

// Метод Эйлера
void euler_callback(double dt, double* state) {
    double newstate[3];
    euler_step(3, state, newstate, dt, rabinovich_fabrikant, NULL);
    memcpy(state, newstate, sizeof(newstate));
}

// Метод Рунге-Кутты 4 с адаптивным шагом
void adaptive_runge_kutta_callback(double dt, double* state) {
    double newstate[3];
    double tolerance = 1e-6;
    adaptive_runge_kutta_step(3, state, newstate, dt, rabinovich_fabrikant, NULL, tolerance);
    memcpy(state, newstate, sizeof(newstate));
}

// Неявный метод Эйлера
void implicit_euler_step_callback(double dt, double* state) {
    double newstate[3];
    int max_iterations = 100;
    double h_numerical = 1e-6;
    implicit_euler_step(state, newstate, dt, rabinovich_fabrikant, NULL, 3, 1e-6, max_iterations, h_numerical);
    memcpy(state, newstate, sizeof(newstate));
}

// Метод предиктор - корректор Адамса - Башфорта - Мултона
void predictor_corrector_callback(double dt, double* state, double* prev_state, double* prev_state2, int step_count) {
    double newstate[3];
    predictor_corrector_step(3, state, newstate, dt, rabinovich_fabrikant, NULL, prev_state, prev_state2);
    memcpy(prev_state2, prev_state, sizeof(prev_state));
    memcpy(prev_state, state, sizeof(state));
    memcpy(state, newstate, sizeof(newstate));
}

// Неявный метод Адамса - Мултона
void adams_multon_callback(double dt, double* state, double* prev_state, double* prev_state2, int step_count) {
    double newstate[3];
    adams_moulton_step(3, state, newstate, dt, rabinovich_fabrikant, NULL, prev_state, prev_state2);
    memcpy(prev_state2, prev_state, sizeof(prev_state));
    memcpy(prev_state, state, sizeof(state));
    memcpy(state, newstate, sizeof(newstate));
}

// Запись в файл
void write_to_file(const char* filename, double t, int step_count, double x, double y, double z) {
    FILE* file;
    errno_t err = fopen_s(&file, filename, "a");
    if (err != 0 || !file) {
        fprintf(stderr, "Ошибка: Не удалось открыть файл '%s' для записи\n", filename);
        return;
    }
    fprintf(file, "%d %.15le %.15le %.15le %.15le\n", step_count, t, x, y, z);
    fclose(file);
}

Uint32 callback(Uint32 interval, void* name) {
    static double dt = 0.01;
    static int step_count = 0;
    static int max_steps = 20000;
    double newstate[3];
    static double prev_state[3];
    static double prev_state2[3];
    static double t = 0.0;
    static double t_max = 60.0;

    //runge_kutta_callback(dt, state); //0.01
    //dormand_prince_callback(dt, state); // 0.01
    //euler_callback(dt, state); // 0.003
    //adaptive_runge_kutta_callback(dt, state); // 0.01
    //predictor_corrector_callback(dt, state, prev_state, prev_state2, step_count);
    //implicit_euler_step_callback(dt, state); // только 0.001
    //adams_multon_callback(dt, state, prev_state, prev_state2, step_count); // 0.003 либо 0.005
    
    //write_to_file("dp8_01.txt", t, step_count, state[0], state[1], state[2]);
    
    t += dt;
    step_count++;
    for (int i = 0; i < 3; ++i) {
        if (isnan(state[i])) {
            fprintf(stderr, "Ошибка: Значения стали NaN на шаге %d\n", step_count);
            return 0;
        }
    }
    for (int i = 0; i < 3; ++i) {
        if (state[i] < min_coords[i]) min_coords[i] = state[i];
        if (state[i] > max_coords[i]) max_coords[i] = state[i];
    }

    printf("Step: %d, State: %.3f, %.3f, %.3f\n", step_count, state[0], state[1], state[2]);
    trajectory.push_back({ state[0], state[1], state[2] });

    if (step_count >= max_steps) {
        return 0;
    }
    /*if (step_count >= max_steps || t >= t_max) {
        return 0;
    }*/
    return interval;
}


static int get_input(void) {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
        case SDL_QUIT:
            return 0;
        case SDL_KEYDOWN:
            if (event.key.keysym.sym == SDLK_ESCAPE) {
                return 0;
            }
            break;
        }
    }
    return 1;
}


void draw_trajectory() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -6.0f);

    glBegin(GL_POINTS);
    //glBegin(GL_LINE_STRIP);
    glColor3f(2.0f, 2.0f, 2.0f);

    if (!trajectory.empty()) {
        for (const auto& point : trajectory) {
            double x = (point[0] - min_coords[0]) / (max_coords[0] - min_coords[0]) * 2.0 - 1.0;
            double y = (point[1] - min_coords[1]) / (max_coords[1] - min_coords[1]) * 2.0 - 1.0;
            double z = (point[2] - min_coords[2]) / (max_coords[2] - min_coords[2]) * 2.0 - 1.0;
            glVertex3f(x, y, z);
        }
    }
    glEnd();

    SDL_GL_SwapWindow(window);
}


void main_loop() {
    while (get_input()) {
        draw_trajectory();
        SDL_Delay(16);
    }
}


void setup_opengl() {
    glViewport(0, 0, screen_size[0], screen_size[1]);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (GLfloat)screen_size[0] / (GLfloat)screen_size[1], 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glEnable(GL_DEPTH_TEST);
}


int main(int argc, char* argv[]) {
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) < 0) {
        fprintf(stderr, "Ошибка инициализации SDL: %s\n", SDL_GetError());
        return EXIT_FAILURE;
    }
    window = SDL_CreateWindow("Rabinovich-Fabrikant",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        screen_size[0], screen_size[1],
        SDL_WINDOW_OPENGL);
    if (!window) {
        fprintf(stderr, "Ошибка создания окна: %s\n", SDL_GetError());
        SDL_Quit();
        return EXIT_FAILURE;
    }
    SDL_GL_CreateContext(window);
    setup_opengl();
    SDL_ShowWindow(window);
    SDL_RaiseWindow(window);
    SDL_AddTimer(16, callback, nullptr);
    main_loop();
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}