#include <iostream> // std::cout, std::endl;

int main()
{
  float x = -1.0f;
  if (x < 0.0f)
    x = -x;
  std::cout << "x = " << x << std::endl; // x = 1
  return 0;
}