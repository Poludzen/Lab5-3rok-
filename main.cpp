#include<iostream>
#include<math.h>
#include<random>
#include<time.h>
#include<sys/wait.h>
#include<unistd.h>
// function in integral
double f(double x) {
	return std::pow(4,x)-std::pow(6,x)+5;
}
// calculate integral by trapezoid method
double calculate_trapezoid(double* pts, int count_of_pts) {
	double result = 0;
	for (int i = 0; i < count_of_pts - 1; i++) {
		result += (f(pts[i + 1]) + f(pts[i])) * (pts[i + 1] - pts[i]) / 2;
	}
	return result;
}


// for gaussian quadrature we need Legandre polynomial roots (and derivative)
double legandre_polynomial(double x, int n) {
	if (n == 0) {
		return 1;
	}
	if (n == 1) {
		return x;
	}
	return (2 * n - 1) * x * legandre_polynomial(x, n - 1) / n - (n - 1) * legandre_polynomial(x, n - 2) / n;
}
double legandre_polynomial_derivative(double x, int n) {
	return n * (legandre_polynomial(x, n - 1) - x * legandre_polynomial(x, n)) / (1 - x * x);
}
void legandre_roots_and_weights(int n, double* arr_of_roots, double* arr_of_weights) {
	// cos(pi) = -1 => arccos(-1) = pi
	double pi = std::acos(-1);
	double eps = 0.001;
	for (int i = 1; i <= n; i++) {
		// finding n roots with Euler method
		double prev = std::cos(pi * (4 * i - 1) / (4 * n + 2));
		double curr = prev - legandre_polynomial(prev, n) / legandre_polynomial_derivative(prev, n);
		while (std::abs(curr - prev) > eps) {
			prev = curr;
			curr = prev - legandre_polynomial(prev, n) / legandre_polynomial_derivative(prev, n);
		}
		arr_of_roots[i - 1] = curr;
		arr_of_weights[i - 1] = 2 / ((1 - curr * curr) * std::pow(legandre_polynomial_derivative(curr, n), 2));
	}
}

// calculate integral by gaussian method
double calculate_gaussian(double* pts, int count_of_pts) {
	double* legandre_points = new double[count_of_pts];
	double* weights = new double[count_of_pts];
	double a = pts[0];
	double b = pts[count_of_pts - 1];
	legandre_roots_and_weights(count_of_pts, legandre_points, weights);
	double result = 0;
	for (int i = 0; i < count_of_pts; i++) {
		result += weights[i] * f(0.5 * (a + b + (b - a) * legandre_points[i]));
	}
	return result * (b - a) / 2;
}

// calculate pi by Lebniz formula
double calculate_pi(int n) {
	double result = 0;
	for (int i = 1; i <= n; i++) {
		result += std::pow(-1, i - 1) / (2 * i - 1);
	}
	return 4 * result;
}
void integral_calculation_task() {
	int P = 3;

	for (int i = 0; i < P; i++) {
		int fork_status = fork();
		// go to child procces
		if (fork_status == 0) {
			// to get new random values at each fork
			srand(time(0) + getpid());
			// let it be -8<=a<8 and a+1<=b<16 but you can choose others 
			int a = rand() % 16 - 8;
			int b = rand() % 8 + a + 1;
			int n_gauss = rand() % 10 + 2;
			//int n_trapez = rand() % 100+2;
			int n_trapez = n_gauss;
			double* gauss_points = new double[n_gauss];
			// [a,b] on n points, so a and b also 2 points
			double gauss_h = double((b - a)) / (n_gauss - 1);
			for (int i = 0; i < n_gauss; i++) {
				gauss_points[i] = a + i * gauss_h;
			}
			double* trapez_points = new double[n_trapez];
			// [a,b] on n points, so a and b also 2 points
			double trapez_h = double((b - a)) / (n_trapez - 1);
			for (int i = 0; i < n_trapez; i++) {
				trapez_points[i] = a + i * trapez_h;
			}
			std::cout << "Calculate with gauss on " << n_gauss << " points" << " on " << a << " and " << b << ":\n";
			std::cout << calculate_gaussian(gauss_points, n_gauss) << "\n";
			std::cout << "Calculate with trapezoid on " << n_trapez << " points" << " on " << a << " and " << b << ":\n";
			std::cout << calculate_trapezoid(trapez_points, n_trapez) << "\n";
			exit(0);
		}
		else if (fork_status < 0) {
			std::cout << "\nERROR!\n";
		}
	}
	for (int i = 0; i < P; i++) {
		wait(NULL);
	}
}
void pi_calculation_task() {
	int P = 3;
	for (int i = 0; i < P; i++) {
		int fork_status = fork();
		// go to child procces
		if (fork_status == 0) {
			// to get new random values at each fork
			srand(time(0) + getpid());
			// to get values from 100 to 5000
			int n = rand()%4900 + 100;
			std::cout << "Calculated pi with Leibniz with n = " << n<<" is : "<<calculate_pi(n)<<"\n";
			exit(0);
		}
		else if (fork_status < 0) {
			std::cout << "\nERROR!\n";
		}
	}
	for (int i = 0; i < P; i++) {
		wait(NULL);
	}
}

int main() {
	integral_calculation_task();
	pi_calculation_task();
}
