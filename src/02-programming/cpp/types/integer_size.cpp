#include <iostream>
#include <limits>

int main()
{
  std::cout << "short: " << sizeof(short) << " B, from "
            << std::numeric_limits<short>::min() << " to "
            << std::numeric_limits<short>::max() << std::endl;
  std::cout << "int: " << sizeof(int) << " B, from "
            << std::numeric_limits<int>::min() << " to "
            << std::numeric_limits<int>::max() << std::endl;
  std::cout << "long: " << sizeof(long) << " B, from "
            << std::numeric_limits<long>::min() << " to "
            << std::numeric_limits<long>::max() << std::endl;
  std::cout << "long long: " << sizeof(long long) << " B, from "
            << std::numeric_limits<long long>::min() << " to "
            << std::numeric_limits<long long>::max() << std::endl;
  std::cout << "unsigned short: " << sizeof(unsigned short)
            << " B, from "
            << std::numeric_limits<unsigned short>::min() << " to "
            << std::numeric_limits<unsigned short>::max()
            << std::endl;
  std::cout << "unsigned: " << sizeof(unsigned) << " B, from "
            << std::numeric_limits<unsigned>::min() << " to "
            << std::numeric_limits<unsigned>::max() << std::endl;
  std::cout << "unsigned long: " << sizeof(unsigned long)
            << " B, from "
            << std::numeric_limits<unsigned long>::min() << " to "
            << std::numeric_limits<unsigned long>::max() << std::endl;
  std::cout << "unsigned long long: " << sizeof(unsigned long long)
            << " B, from "
            << std::numeric_limits<unsigned long long>::min()
            << " to "
            << std::numeric_limits<unsigned long long>::max()
            << std::endl;
  return 0;
}