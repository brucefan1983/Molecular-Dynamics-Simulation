#include <iomanip>  // std::setprecision
#include <iostream> // std::cout, std::fixed, std::endl;

int main()
{
  // make sure to output 20 digits after the decimal point
  std::cout << std::fixed << std::setprecision(20);

  float pi1{3.14159265358979323846f};
  double pi2{3.14159265358979323846};
  long double pi3{3.14159265358979323846L};
  auto pi4 = 3.14159265358979323846f;
  auto pi5 = 3.14159265358979323846;
  auto pi6 = 3.14159265358979323846L;

  std::cout << pi1 << std::endl; // 3.14159274101257324219
  std::cout << pi2 << std::endl; // 3.14159265358979311600
  std::cout << pi3 << std::endl; // 3.14159265358979323851
  std::cout << pi4 << std::endl; // 3.14159274101257324219
  std::cout << pi5 << std::endl; // 3.14159265358979311600
  std::cout << pi6 << std::endl; // 3.14159265358979323851

  return 0;
}