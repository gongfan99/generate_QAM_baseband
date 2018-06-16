#include <iostream>
#include <random>
#include <thread>

#include "ThreadPool.h"

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_int_distribution<int> dist(-500, 500);
auto rnd = std::bind(dist, mt);


void simulate_hard_computation() {
  std::this_thread::sleep_for(std::chrono::milliseconds(500 + rnd()));
}

// Simple function that adds multiplies two numbers and prints the result
void multiply(const int a, const int b) {
  simulate_hard_computation();
  const int res = a * b;
  std::cout << a << " * " << b << " = " << res << std::endl;
}

// Same as before but now we have an output parameter
void multiply_output(int & out, const int a, const int b) {
  simulate_hard_computation();
  out = a * b;
  std::cout << a << " * " << b << " = " << out << std::endl;
}

// Same as before but now we have an output parameter
int multiply_return(const int a, const int b) {
  simulate_hard_computation();
  const int res = a * b;
  std::cout << a << " * " << b << " = " << res << std::endl;
  return res;
}


int main() {
  // Create pool with 3 threads
  ThreadPool pool(3);

  // Initialize pool
  pool.start();
  std::cout << "pool started" << std::endl;

  // Submit (partial) multiplication table
  for (int i = 1; i < 3; ++i) {
    for (int j = 1; j < 4; ++j) {
      pool.submit(multiply, i, j);
      std::cout << i << " * " << j << " submitted" << std::endl;
    }
  }

  // Submit function with output parameter passed by ref
  int output_ref;
  pool.submit(multiply_output, std::ref(output_ref), 5, 6);
  std::cout << 5 << " * " << 6 << " submitted" << std::endl;

  // Submit (partial) multiplication table
  for (int i = 1; i < 3; ++i) {
    for (int j = 1; j < 4; ++j) {
      pool.submitFuture(multiply, i, j);
      std::cout << i << " * " << j << " submitted" << std::endl;
    }
  }

  // Submit function with output parameter passed by ref
  auto ft = pool.submitFuture(multiply_return, 5, 6);
  std::cout << 5 << " * " << 6 << " submitted" << std::endl;
  std::cout << 5 << " * " << 6 << " = " << ft.get() << std::endl;

  auto ft2 = pool.submitFuture(multiply_output, std::ref(output_ref), 7, 8);
  std::cout << 7 << " * " << 8 << " submitted" << std::endl;
  ft2.wait();
  std::cout << 7 << " * " << 8 << " = " << output_ref << std::endl;

  std::this_thread::sleep_for(std::chrono::milliseconds(2000));
  std::cout << "pool shutting down" << std::endl;
  pool.shutdown();
  std::cout << "pool shut down" << std::endl;
}