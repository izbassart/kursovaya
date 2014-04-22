#include "magnit.h"

short Magnit :: PA[] = { 0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32, 39, 48,
                         56, 67, 80, 94, 111, 132, 154, 179, 207, 236, 300, 400 };

double Magnit :: ap2kp(double ap)
{
  int j = 0;
  while (ap < PA[j + 1])
      j++;
  return (j - 1 + (ap - PA[j]) / (PA[j + 1] - PA[j])) / 3;
}
