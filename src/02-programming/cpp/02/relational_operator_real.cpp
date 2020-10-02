#include <iomanip>  // std::setprecision
#include <iostream> // std::cout, std::fixed, std::endl;

int main()
{
  double x = 3.0;
  double y = 1.0 / 3.0;
  double z = 1.0;

  // make sure to output 16 digits after the decimal point
  std::cout << std::fixed << std::setprecision(16);

  bool relation = (z == x * y);        // false
  std::cout << relation << std::endl;  // 0
  std::cout << z - x * y << std::endl; // 0.0000000000000001

  return 0;
}