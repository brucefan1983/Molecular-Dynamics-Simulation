#include <iomanip>  // std::setprecision
#include <iostream> // std::cout, std::fixed, std::endl;

int main()
{
  // make sure to output 15 digits after the decimal point
  std::cout << std::fixed << std::setprecision(15);

  // original value
  double double_variable = 3.141592653589793; // double
  std::cout << double_variable << std::endl;  // 3.141592653589793

  // using =
  float float_variable = double_variable;   // double->float
  std::cout << float_variable << std::endl; // 3.141592741012573
  int int_variable = double_variable;       // doubel->int
  std::cout << int_variable << std::endl;   // 3

  // using {}
  float float_variable_2{double_variable}; // warning or error
  int int_variable_2{double_variable};     // warning or error
  return 0;
}