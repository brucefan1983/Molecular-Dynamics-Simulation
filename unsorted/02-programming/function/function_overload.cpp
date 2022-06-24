#include <iostream> // std::cout, std::endl;

void my_swap(float* x, float* y) // float type
{
  float z = *x;
  *x = *y;
  *y = z;
}

void my_swap(double* x, double* y) // double type
{
  double z = *x;
  *x = *y;
  *y = z;
}

void my_swap(int* x, int* y) // int type
{
  int z = *x;
  *x = *y;
  *y = z;
}

int main()
{
  float xf = 1.0f;
  float yf = 2.0f;
  my_swap(&xf, &yf);
  std::cout << "xf = " << xf << ", yf = " << yf << std::endl;

  double xd = 1.0;
  double yd = 2.0;
  my_swap(&xd, &yd);
  std::cout << "xd = " << xd << ", yd = " << yd << std::endl;

  int xi = 1;
  int yi = 2;
  my_swap(&xi, &yi);
  std::cout << "xi = " << xi << ", yi = " << yi << std::endl;

  return 0;
}
