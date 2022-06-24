#include <iostream> // std::cout, std::endl;

int main()
{
  int n = 6;
  for (int i = -n; i <= n; ++i) {
    if (i == 0)
      continue;
    std::cout << 1.0 / static_cast<double>(i) << std::endl;
  }

  return 0;
}