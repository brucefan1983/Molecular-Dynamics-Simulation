#include <iostream> // std::cout, std::endl;

int main()
{
  float x[] = {1.0f, 2.0f, 3.0f, 4.0f};
  for (int n = 0; n < 4; ++n) {
    std::cout << "x[" << n << "] = " << x[n] << std::endl;
    std::cout << "*(x + " << n << ") = " << *(x + n) << std::endl;
  }
  return 0;
}
