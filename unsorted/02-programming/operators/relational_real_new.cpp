#include <iomanip>  // std::setprecision
#include <iostream> // std::cout, std::fixed, std::endl;

int main()
{
  double x = 3.0;
  double y = 1.0 / 3.0;
  double z = 1.0;
  double x_times_y = x * y;
  double epsilon = 1.0e-15;

  bool relation_1 = x_times_y < z + epsilon; // true
  bool relation_2 = x_times_y > z - epsilon; // true
  bool relation = relation_1 && relation_2;  // true
  std::cout << relation << std::endl;        // 0

  return 0;
}