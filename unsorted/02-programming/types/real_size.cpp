#include <iostream>
#include <limits>

int main()
{
  std::cout << "float: " << sizeof(float) << " B, from "
            << std::numeric_limits<float>::min() << " to "
            << std::numeric_limits<float>::max() << std::endl;
  std::cout << "double: " << sizeof(double) << " B, from "
            << std::numeric_limits<double>::min() << " to "
            << std::numeric_limits<double>::max() << std::endl;
  std::cout << "long double: " << sizeof(long double) << " B, from "
            << std::numeric_limits<long double>::min() << " to "
            << std::numeric_limits<long double>::max() << std::endl;

  return 0;
}