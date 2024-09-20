#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

const double M_PI = acos(-1);

class FFT {
public:
    // Основной алгоритм FFT для длины, кратной 2
    static void fft(vector<complex<double>>& x) {
        int n = x.size();
        if (n <= 1) return;

        // Разделяем массив на четные и нечетные элементы
        vector<complex<double>> even(n / 2);
        vector<complex<double>> odd(n / 2);
        for (int i = 0; i < n / 2; ++i) {
            even[i] = x[i * 2];
            odd[i] = x[i * 2 + 1];
        }

        fft(even); // Рекурсивно применяем FFT к четным элементам
        fft(odd);  // Рекурсивно применяем FFT к нечетным элементам

        for (int k = 0; k < n / 2; ++k) {
            complex<double> t = polar(1.0, -2 * M_PI * k / n) * odd[k];
            x[k] = even[k] + t;
            x[k + n / 2] = even[k] - t;
        }
    }

    static void inverse(vector<complex<double>>& x) {
        int n = x.size();
        for (int i = 0; i < n; ++i) {
            x[i] = conj(x[i]); // Конъюгирование
        }
        fft(x);
        for (int i = 0; i < n; ++i) {
            x[i] = conj(x[i]) / static_cast<double>(n);
        }
    }
};

// FFT для длины, кратной 3
class FFT3 {
private:
    static void fft3(vector<complex<double>>& x) {
        int n = x.size();
        if (n <= 1) return;

        vector<complex<double>> y(n);
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                y[k] += x[j] * polar(1.0, -2 * M_PI * k * j / n);
            }
        }

        x = y;
    }

public:
    static void forward(vector<complex<double>>& x) {
        fft3(x);
    }

    static void inverse(vector<complex<double>>& x) {
        for (int i = 0; i < x.size(); ++i) {
            x[i] = conj(x[i]);
        }
        fft3(x);
        for (int i = 0; i < x.size(); ++i) {
            x[i] = conj(x[i]) / static_cast<double>(x.size());
        }
    }
};

// FFT для длины, кратной 5
class FFT5 {
public:
    static void forward(vector<complex<double>>& x) {
        int n = x.size();
        if (n % 5 != 0) throw invalid_argument("Size of input must be a multiple of 5.");

        vector<complex<double>> y(n);
        for (int k = 0; k < n; ++k) {
            for (int j = 0; j < n; ++j) {
                y[k] += x[j] * polar(1.0, -2 * M_PI * k * j / n);
            }
        }
        x = y;
    }

    static void inverse(vector<complex<double>>& x) {
        int n = x.size();
        for (int i = 0; i < n; ++i) {
            x[i] = conj(x[i]);
        }
        forward(x);
        for (int i = 0; i < n; ++i) {
            x[i] = conj(x[i]) / static_cast<double>(n);
        }
    }
};

// Основная функция
int main() {
    srand(static_cast<unsigned int>(time(0)));

    int N;
    cout << "Enter the array size (must be a multiple of 2, 3 or 5):";
    cin >> N;

    if (N <= 0) {
        cout << "The array size must be a positive number." << endl;
        return 1;
    }
    vector<complex<double>> input(N);

    // Заполнение случайными комплексными числами
    for (int i = 0; i < N; ++i) {
        input[i] = complex<double>(rand() % 100, rand() % 100);
    }

    // Вывод входных данных
    cout << "Input data:\n";
    for (auto& val : input) {
        cout << val << "\n";
    }

    // Выполнение прямого преобразования
    vector<complex<double>> output = input;

    if (N % 5 == 0) {
        FFT5::forward(output);
    }
    else if (N % 3 == 0) {
        FFT3::forward(output);
    }
    else if (N % 2 == 0) {
        FFT::fft(output);
    }
    else {
        cout << "Invalid size for FFT. Size must be a multiple of 2, 3, or 5.\n";
        return 1;  // Завершение программы из-за неверного размера
    }

    // Вывод результата прямого преобразования
    cout << "\nFFT output:\n";
    for (auto& val : output) {
        cout << val << "\n";
    }

    // Выполнение обратного преобразования
    if (N % 5 == 0) {
        FFT5::inverse(output);
    }
    else if (N % 3 == 0) {
        FFT3::inverse(output);
    }
    else if (N % 2 == 0) {
        FFT::inverse(output);
    }

    // Вывод результата обратного преобразования
    cout << "\nInverse FFT output (should be close to input):\n";
    for (auto& val : output) {
        cout << val << "\n";
    }

    // Сравнение ошибки между входными и выходными данными
    double error = 0.0;
    for (int i = 0; i < N; ++i) {
        error += abs(input[i] - output[i]);
    }
    cout << "\nTotal absolute error: " << error << "\n";

    return 0;
}