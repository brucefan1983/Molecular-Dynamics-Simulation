#include <iostream> // std::cout, std::endl;

void my_swap1(float x, float y) // wrong implementation
{
  float z = x;
  x = y;
  y = z;
}

void my_swap2(float* x, float* y) // using pointer
{
  float z = *x;
  *x = *y;
  *y = z;
}

void my_swap3(float& x, float& y) // using reference
{
  float z = x;
  x = y;
  y = z;
}

int main()
{
  float x = 1.0f;
  float y = 2.0f;
  my_swap1(x, y);
  std::cout << "x = " << x << ", y = " << y << std::endl;
  my_swap2(&x, &y);
  std::cout << "x = " << x << ", y = " << y << std::endl;
  my_swap3(x, y);
  std::cout << "x = " << x << ", y = " << y << std::endl;
  return 0;
}
