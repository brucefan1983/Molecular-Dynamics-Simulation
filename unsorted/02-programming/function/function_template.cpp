#include <iostream> // std::cout, std::endl;

template <typename T>
void my_swap(T* x, T* y)
{
  T z = *x;
  *x = *y;
  *y = z;
}

int main()
{
  float xf = 1.0f;
  float yf = 2.0f;
  my_swap<float>(&xf, &yf);
  std::cout << "xf = " << xf << ", yf = " << yf << std::endl;

  double xd = 1.0;
  double yd = 2.0;
  my_swap<double>(&xd, &yd);
  std::cout << "xd = " << xd << ", yd = " << yd << std::endl;

  int xi = 1;
  int yi = 2;
  my_swap<int>(&xi, &yi);
  std::cout << "xi = " << xi << ", yi = " << yi << std::endl;

  return 0;
}
