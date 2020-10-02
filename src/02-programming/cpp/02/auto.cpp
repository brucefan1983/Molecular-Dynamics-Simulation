#include <iomanip>  // std::setprecision
#include <iostream> // std::cout, std::fixed, std::endl;

int main()
{
  // make sure to output 20 digits after the decimal point
  std::cout << std::fixed << std::setprecision(20);

  auto float_variable = 3.14159265358979323846f;       // suffix f or F
  auto double_variable = 3.14159265358979323846;       // no suffix
  auto long_double_variable = 3.14159265358979323846L; // suffix L or l

  std::cout << float_variable << std::endl;       // 3.14159274101257324219
  std::cout << double_variable << std::endl;      // 3.14159265358979311600
  std::cout << long_double_variable << std::endl; // 3.14159265358979323851

  return 0;
}