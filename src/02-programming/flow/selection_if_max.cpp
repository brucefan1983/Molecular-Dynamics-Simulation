#include <iostream> // std::cout, std::endl;

int main()
{
  double a = 1.0;
  double b = 2.0;

  // method 1: using if-else
  double max1 = 0.0;
  if (a > b) {
    max1 = a;
  } else {
    max1 = b;
  }
  std::cout << "max1 = " << max1 << std::endl; // max1 = 2

  // method 2: using conditional expression
  auto max2 = (a > b) ? a : b;
  std::cout << "max2 = " << max2 << std::endl; // max2 = 2

  return 0;
}